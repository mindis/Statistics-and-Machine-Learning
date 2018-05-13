#####Factorized Latent Dirichlet Allocation#####
options(warn=0)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####データの発生####
##データの設定
k <- 15   #トピック数
hh <- 5000   #ユーザー数
d <- 2500   #アイテム数
v <- 1000   #語彙数 
w <- rpois(d, rgamma(d, 50, 0.5))   #1文書あたりの語彙数
f <- sum(w)   #総語彙数

#IDの設定
d_id <- rep(1:d, w)   #文書ID
user_id0 <- rep(1:hh, rep(d, hh))
item_id0 <- rep(1:d, hh)


####トピックモデルのデータを生成####
##パラメータを設定
#ディリクレ事前分布のパラメータ
alpha1 <- rep(0.15, k)
alpha2 <- rep(0.1, v)

#全種類の単語が出現するまで繰り返す
for(rp in 1:1000){
  print(rp)
  
  #ディリクレ分布からパラメータを生成
  theta <- thetat <- extraDistr::rdirichlet(d, alpha1)
  phi <- phit <- extraDistr::rdirichlet(k, alpha2)
  
  ##トピックと単語を生成
  word_list <- Z_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  
  for(i in 1:d){
    #トピックを生成
    z <- rmnom(w[i], 1, theta[i, ])
    z_vec <- as.numeric(z %*% 1:k)
    
    #単語を生成
    word <- rmnom(w[i], 1, phi[z_vec, ])
    WX[i, ] <- colSums(word)
    
    #データを格納
    word_list[[i]] <- as.numeric(word %*% 1:v)
    Z_list[[i]] <- z
  }
  if(min(colSums(WX)) > 0) break
}

#リストを変換
wd <- unlist(word_list)
Z <- do.call(rbind, Z_list)

#インデックスを設定
doc_list <- doc_vec <- list()
wd_list <- wd_vec <- list()
for(i in 1:d){
  doc_list[[i]] <- which(d_id==i)
  doc_vec[[i]] <- rep(1, length(doc_list[[i]]))
}
for(j in 1:v){
  wd_list[[j]] <- which(wd==j)
  wd_vec[[j]] <- rep(1, length(wd_list[[j]]))
}

#トピック分布の期待値を設定
Z_score <- matrix(0, nrow=d, ncol=k) 
for(i in 1:d){
  Z_score[i, ] <- colMeans(Z[doc_list[[i]], ])
}


####行列分解のデータ生成####
##パラメータを設定
N <- hh*d
mu <- mut <- 5.5
tau <- taut <- diag(3.5, k)
sigma <- sigmat <- 0.5   #モデルの標準偏差
U <- Ut <- mvrnorm(hh, rep(mu, k), tau)   #パラメータ

##正規分布からスコアを生成
Mu <- as.numeric(t(U %*% t(Z_score)))
y0 <- rnorm(N, Mu, sigma)

#生成したスコアを評価データに変換
y0_censor <- ifelse(y0 < 1, 1, ifelse(y0 > 10, 10, y0)) 
y_full <- y0 #round(y0_censor, 0)   #スコアを丸める


##欠損ベクトルを生成
#欠損有無のベータ分布のパラメータを設定
beta1 <- rbeta(hh, 5.5, 10.0)   #ユーザ-購買確率
beta2 <- rbeta(d, 5.0, 8.0)   #アイテム購買確率

#欠損がある購買データを生成
Zt <- matrix(0, nrow=hh, ncol=d)
for(j in 1:d){
  deficit <- rbinom(hh, 1, beta1 * beta2[j])
  Zt[, j] <- deficit   #欠損を代入
}

#欠損インデックス
zt_vec <- as.numeric(t(Zt))
index_zt1 <- which(zt_vec==1)
index_zt0 <- which(zt_vec==0)
N <- length(index_zt1)

#欠損ベクトルに応じてデータを抽出
user_id <- user_id0[index_zt1]
item_id <- item_id0[index_zt1]
y <- y_full[index_zt1]
n1 <- plyr::count(user_id)$freq
n2 <- plyr::count(item_id)$freq

#生成した応答変数のヒストグラム
hist(y0, col="grey", xlab="スコア", main="ユーザー×アイテムのスコア分布")   #元データ
hist(y_full, col="grey", xlab="スコア", main="ユーザー×アイテムのスコア分布")   #完全データのスコア分布
hist(y, col="grey", xlab="スコア", main="ユーザー×アイテムのスコア分布")   #購買データのスコア分布

####モンテカルロEMアルゴリズムでfLDAを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k){
  #負担係数を計算
  Bur <- theta[w, ] * t(phi)[wd, ]   #尤度
  Br <- Bur / rowSums(Bur)   #負担率
  r <- colSums(Br) / sum(Br)   #混合率
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
tol <- 1
iter <- 1
dl <- 100
L <- 500   #モンテカルロサンプリング数

##事前分布の設定
alpha1 <- 0.1
alpha2 <- 0.1
beta1 <- 1
beta2 <- 1

##パラメータの真値
theta <- thetat
phi <- phit
sigma <- sigmat
mu <- rep(mut, k)
tau <- taut
inv_tau <- solve(tau)
u <- Ut
Zi <- Z
Zi_T <- t(Zi)
wsum0 <- matrix(0, nrow=d, ncol=k)
for(i in 1:d){
  wsum0[i, ] <- Zi_T[, doc_list[[i]]] %*% doc_vec[[i]]
}
Zi_score <- wsum0 / w


##初期値の設定
theta <- extraDistr::rdirichlet(d, rep(1.0, k))
phi <- extraDistr::rdirichlet(k, rep(1.0, v))
sigma <- 0.5
mu <- rep(mean(y), k)
tau <- diag(0.3, k)
inv_tau <- solve(tau)
u <- mvrnorm(hh, mu, tau)
Zi <- rmnom(f, 1, burden_fr(theta, phi, wd, d_id, k)$Br)
Zi_T <- t(Zi)
wsum0 <- matrix(0, nrow=d, ncol=k)
for(i in 1:d){
  wsum0[i, ] <- Zi_T[, doc_list[[i]]] %*% doc_vec[[i]]
}
Zi_score <- wsum0 / w


##インデックスを作成
user_list <- item_list <- list()
ui_list <- iu_list <- y_list <- y_vec <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
  iu_list[[i]] <- item_id[user_list[[i]]]
}
for(j in 1:d){
  item_list[[j]] <- which(item_id==j)
  ui_list[[j]] <- user_id[item_list[[j]]]
  y_list[[j]] <- y[item_list[[j]]]
  y_vec[[j]] <- rep(1, length(y_list[[j]]))
}
w_vec <- w[d_id]



####モンテカルロEMアルゴリズムでパラメータを推定####
while(abs(dl) > tol){   #dlがtol以上なら繰り返す

  ###モンテカルロEステップで潜在変数をサンプリング
  ##トピック尤度からトピック割当確率を推定
  #完全データの対数尤度
  z_par <- theta[d_id, ] * t(phi)[wd, ]

  #トピックの割当からトピックをサンプリング
  z_rate <- z_par / rowSums(z_par)   #トピックの割当確率
  Zi <- rmnom(f, L, z_rate) / L
  Zi_T <- t(Zi)
  
  ##ユーザー特徴行列をサンプリング
  #トピックスコアを設定
  wsum0 <- matrix(0, nrow=d, ncol=k)
  for(i in 1:d){
    wsum0[i, ] <- Zi_T[, doc_list[[i]]] %*% doc_vec[[i]]
  }
  Zi_score <- (wsum0 / w)
  
  #ユーザーごとに特徴ベクトルをサンプリング
  for(i in 1:hh){
    
    #特徴ベクトルの事後分布のパラメータ
    index <- iu_list[[i]]   #アイテムインデックス
    Xy <- sigma^-2 * (t(Zi_score[index, ]) %*% y[user_list[[i]]])
    XXV <- sigma^-2 * (t(Zi_score[index, ]) %*% Zi_score[index, ]) + inv_tau
    inv_XXV <- solve(XXV)
    x <- inv_XXV %*% (Xy + inv_tau %*% mu)   #事後分布の平均
    
    #多変量正規分布からユーザー特徴ベクトルをサンプリング
    u[i, ] <- colMeans(mvrnorm(L, x, inv_XXV))   #モンテカルロ平均
  }

  ###Mステップで完全データの尤度を最大化
  ##観測モデルの誤差パラメータを更新
  uz_mu <- as.numeric(t(u %*% t(Zi_score)))[index_zt1]
  er <- y - uz_mu
  sigma <- sd(er)

  
  ##トピックモデルのパラメータを推定
  #トピック分布のパラメータを更新
  wsum <- wsum0 + alpha1
  theta <- wsum / (w + alpha1*k)
  
  #単語分布のパラメータを更新
  vsum0 <- matrix(0, nrow=k, ncol=v)
  for(j in 1:v){
    vsum0[, j] <- Zi_T[, wd_list[[j]], drop=FALSE] %*% wd_vec[[j]]
  }
  vsum <- vsum0 + alpha2
  phi <- vsum / as.numeric(vsum %*% rep(1, v))
  
  
  ##ユーザー特徴行列の階層モデルのパラメータを推定
  mu <- colMeans(u)   #平均ベクトル
  tau <- diag(diag(var(u) + alpha1))   #分散共分散行列
  inv_tau <- solve(tau)
  
  ##アルゴリズムの収束判定
  LLs1 <- sum(dnorm(y, uz_mu, sigma, log=TRUE))   #完全データの対数尤度
  LLs2 <- sum(log(rowSums(theta[d_id, ] * t(phi)[wd, ])))
  LL <- LLs1 + LLs2
  iter <- iter + 1
  dl <- LL - LL1
  LL1 <- LL
  print(c(LL, LLs1, LLs2))
}



