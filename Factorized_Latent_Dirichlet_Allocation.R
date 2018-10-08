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
item <- 2500   #アイテム数
v <- 1000   #語彙数 
w <- rpois(item, rgamma(item, 60, 0.4))   #1文書あたりの語彙数
pt <- rtpois(hh, rgamma(hh, 25.0, 0.225), a=1, b=Inf)   #購買接触数
f <- sum(w)   #総語彙数
hhpt <- sum(pt)   #総スコア数
vec_k <- rep(1, k)

#IDの設定
d_id <- rep(1:item, w)   #文書ID
no_id <- as.numeric(unlist(tapply(1:f, d_id, rank)))
user_id <- rep(1:hh, pt)   #ユーザーID
t_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}

##アイテムの割当を生成
#セグメント割当を生成
topic <- 25
phi <- extraDistr::rdirichlet(topic, rep(0.5, item))
z <- as.numeric(rmnom(hh, 1,  extraDistr::rdirichlet(hh, rep(2.5, topic))) %*% 1:topic)

#多項分布からアイテムを生成
item_id_list <- list()
for(i in 1:hh){
  if(i%%100==0){
    print(i)
  }
  item_id_list[[i]] <- as.numeric(rmnom(pt[i], 1, phi[z[user_id[user_list[[i]]]], ]) %*% 1:item)
}
item_id <- unlist(item_id_list)
item_list <- list()
for(j in 1:item){
  item_list[[j]] <- which(item_id==j)
}

#スパース行列を作成
user_data <- sparseMatrix(1:hhpt, user_id, x=rep(1, hhpt), dims=c(hhpt, hh))
user_data_T <- t(user_data)
item_data <- sparseMatrix(1:hhpt, item_id, x=rep(1, hhpt), dims=c(hhpt, item))
item_data_T <- t(item_data)

#生成したデータを可視化
freq_item <- plyr::count(item_id); freq_item$x <- as.character(freq_item$x)
hist(freq_item$freq, breaks=25, col="grey", xlab="アイテムの購買頻度", main="アイテムの購買頻度分布")
gc(); gc()


##素性ベクトルを生成
k1 <- 3; k2 <- 5; k3 <- 5
x1 <- matrix(runif(hhpt*k1, 0, 1), nrow=hhpt, ncol=k1)
x2 <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #データを結合
col_x <- ncol(x)

##階層モデルの説明変数を生成
#ユーザーの説明変数
k1 <- 2; k2 <- 4; k3 <- 5
u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #データを結合
col_u <- ncol(u)

#アイテムの説明変数
k1 <- 2; k2 <- 4; k3 <- 4
g1 <- matrix(runif(item*k1, 0, 1), nrow=item, ncol=k1)
g2 <- matrix(0, nrow=item, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  g2[, j] <- rbinom(item, 1, pr)
}
g3 <- rmnom(item, 1, runif(k3, 0.2, 1.25)); g3 <- g3[, -which.min(colSums(g3))]
g <- cbind(1, g1, g2, g3)   #データを結合
col_g <- ncol(g)


####トピックモデルのデータを生成####
##パラメータを設定
#ディリクレ事前分布のパラメータ
alpha1 <- rep(0.15, k)
alpha2 <- rep(0.05, v)

#全種類の単語が出現するまで繰り返す
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  
  #ディリクレ分布からパラメータを生成
  theta <- thetat <- extraDistr::rdirichlet(item, alpha1)
  phi <- extraDistr::rdirichlet(k, alpha2)
  
  #単語出現確率が低いトピックを入れ替える
  index <- which(colMaxs(phi) < (k*10)/f)
  for(j in 1:length(index)){
    phi[as.numeric(rmnom(1, 1, extraDistr::rdirichlet(1, rep(2.0, k))) %*% 1:k), index[j]] <- (k*10)/f
  }
  phit <- phi

  ##トピックと単語を生成
  word_list <- Z_list <- list()
  WX <- matrix(0, nrow=item, ncol=v)
  
  for(i in 1:item){
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

#スパース行列を作成
d_data <- sparseMatrix(1:f, d_id, x=rep(1, f), dims=c(f, item))
d_data_T <- t(d_data)
word_data <- sparseMatrix(1:f, wd, x=rep(1, f), dims=c(f, v))
word_data_T <- t(word_data)

#インデックスを設定
doc_list <- wd_list <- list()
for(i in 1:item){
  doc_list[[i]] <- which(d_id==i)
}
for(j in 1:v){
  wd_list[[j]] <- which(wd==j)
}

#トピック分布の期待値を設定
Z_score <- matrix(0, nrow=item, ncol=k) 
for(i in 1:d){
  Z_score[i, ] <- colMeans(Z[doc_list[[i]], ])
}


####評価ベクトルを生成####
rp <- 0
repeat { 
  rp <- rp + 1
  
  ##パラメータを設定
  #素性ベクトルのパラメータ
  sigma <- sigmat <- 0.5
  beta <- betat <- c(5.5, rnorm(col_x-1, 0, 0.75))
  
  #階層モデルの分散パラメータ
  Cov_u <- Cov_ut <- runif(1, 0.1, 0.4)   #ユーザー-アイテムの階層モデルの標準偏差
  Cov_v <- Cov_vt <- runif(1, 0.1, 0.4)   #アイテムの階層モデルの標準偏差
  Cov_z <- Cov_zt <- diag(runif(k, 0.01, 0.1), k)   #トピックのユーザー特徴ベクトルの階層モデルの分散
  
  #階層モデルの回帰係数を設定
  alpha_u <- alpha_ut <- rnorm(col_u, 0, 0.35)
  alpha_v <- alpha_vt <- rnorm(col_g, 0, 0.35)
  alpha_z <- alpha_zt <- mvrnorm(col_u, rep(0, k), runif(k, 0.1, 0.25) * diag(k))
  
  #変量効果と特徴ベクトルのパラメータを生成
  theta_u <- theta_ut <- u %*% alpha_u + rnorm(hh, 0, Cov_u)
  theta_v <- theta_vt <- g %*% alpha_v + rnorm(item, 0, Cov_v)
  theta_z <- theta_zt <- u %*% alpha_z + mvrnorm(hh, rep(0, k), Cov_z)
  
  
  #正規分布からスコアを生成
  mu <- as.numeric(x %*% beta + theta_u[user_id] + theta_v[item_id] + (Z_score[item_id, ] * theta_z[item_id, ]) %*% vec_k)
  y0 <- rnorm(hhpt, mu, sigma)
  
  #break条件
  print(round(c(max(y0), min(y0)), 3))
  if(max(y0) < 15.0 & min(y0) > -4.0 & max(y0) > 11.0 & min(y0) < -1.0){
    break
  }
}

#生成したスコアを評価データに変換
y0_censor <- ifelse(y0 < 1, 1, ifelse(y0 > 10, 10, y0)) 
y <- round(y0_censor, 0)   #スコアを丸める

#スコア分布と単語分布
hist(y0, col="grey", breaks=25, xlab="スコア", main="完全データのスコア分布")
hist(y, col="grey", breaks=25, xlab="スコア", main="観測されたスコア分布")


####マルコフ連鎖モンテカルロ法でfLDAを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k, vec_k){
  #負担係数を計算
  Bur <- theta[w, ] * t(phi)[wd, ]   #尤度
  Br <- Bur / as.numeric(Bur %*% vec_k)   #負担率
  bval <- list(Br=Br, Bur=Bur)
  return(bval)
}

##アルゴリズムの設定
R <- 3000   #サンプリング回数
keep <- 2   #2回に1回の割合でサンプリング結果を格納
disp <- 10
iter <- 0
burnin <- 1000/keep


##事前分布の設定
#トピックモデルの事前分布
alpha01 <- 0.1
alpha02 <- 0.1

#行列分解の事前分布
beta01 <- 0
beta02 <- 0
s0 <- 0.1
v0 <- 0.1
Deltabar <- matrix(0, nrow=ncol(u), ncol=k)   #階層モデルの回帰係数の事前分布の平均
ADelta <- 0.01 * diag(1, ncol(u))   #階層モデルの回帰係数の事前分布の分散
nu <- k + 1   #逆ウィシャート分布の自由度
V <- nu1 * diag(rep(1, k)) #逆ウィシャート分布のパラメータ


##パラメータの真値
#トピックモデルのパラメータ
theta <- thetat
phi <- phit
wsum <- as.matrix(d_data_T %*% Z)
Z_score <- wsum / w

#素性ベクトルのパラメータ
sigma <- sigmat
beta <- betat

#階層モデルの分散パラメータ
Cov_u <- Cov_ut
Cov_v <- Cov_vt   
Cov_z <- Cov_zt  

#階層モデルの回帰係数を設定
alpha_u <- alpha_ut
alpha_v <- alpha_vt
alpha_z <- alpha_zt

#変量効果と特徴ベクトルのパラメータ
theta_u <- theta_ut
theta_v <- theta_vt
theta_z <- theta_zt


##初期値の設定
#トピックモデルの初期値
theta <- extraDistr::rdirichlet(item, rep(1.0, k))
phi <- extraDistr::rdirichlet(k, rep(1.0, v))
Zi <- rmnom(f, 1, rep(1/k, k))
wsum <- as.matrix(d_data_T %*% Z)
Z_score <- wsum / w

#素性ベクトルのパラメータ
sigma <- 1.0
beta <- as.numeric(solve(t(x) %*% x) %*% t(x) %*% y)

#階層モデルの分散パラメータ
Cov_u <- 0.2
Cov_v <- 0.2
Cov_z <- 0.01 * diag(k)

#階層モデルの回帰係数を設定
alpha_u <- rnorm(col_u, 0, 0.1)
alpha_v <- rnorm(col_g, 0, 0.1)
alpha_z <- mvrnorm(col_u, rep(0, k), 0.01 * diag(k))

#変量効果と特徴ベクトルのパラメータを生成
theta_u <- u %*% alpha_u + rnorm(hh, 0, Cov_u)
theta_v <- g %*% alpha_v + rnorm(item, 0, Cov_v)
theta_z <- u %*% alpha_z + mvrnorm(hh, rep(0, k), Cov_z)



####ギブスサンプリングムでパラメータをサンプリング####
for(rp in 1:R){
  
  
  ##トピックごとの評価スコアの尤度を設定
  #データの設定
  mu_uv <- as.numeric(x %*% beta + theta_u[user_id] + theta_v[item_id])   #トピック因子を除いた期待値
  wsum_z <- wsum[d_id, ] - Z   #単語ベクトルでのトピック割当
  w_vec <- w[d_id]
  
  wsum_topic <- wsum_z
  wsum_topic[, j]  <- wsum_topic[, j] + 1
  Z_score <- wsum_topic / w_vec

  
  #トピック尤度からトピック割当確率を推定
  theta_z
  
  
  
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
  LLs <- c(LLs, LL1)
  print(c(LL, LLs1, LLs2))
}



