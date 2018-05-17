#####Switching Multinomial LDA model#####
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
r <- 5   #評価スコア数
s <- 3   #極性値数
a <- 3   #分岐数
k11 <- 10   #ユーザーのテキストのトピック数
k12 <- 15   #アイテムのテキストのトピック数
hh <- 1000   #レビュアー数
item <- 200   #アイテム数
v1 <- 300   #評価スコアの語彙数
v2 <- 350   #ユーザートピックの語彙数
v3 <- 350   #アイテムトピックの語彙数
v <- v1 + v2 + v3   #総語彙数
spl <- matrix(1:v1, nrow=s, ncol=v1/s, byrow=T)
v1_index <- 1:v1
v2_index <- (v1+1):v2
v3_index <- (v2+1):v

##IDと欠損ベクトルの作成
#IDを仮設定
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)

#欠損ベクトルを作成
for(rp in 1:100){
  m_vec <- rep(0, hh*item)
  for(i in 1:item){
    prob <- runif(1, 0.025, 0.16)
    m_vec[item_id0==i] <- rbinom(hh, 1, prob)
  }
  m_index <- which(m_vec==1)
  
  #完全なIDを設定
  user_id <- user_id0[m_index]
  item_id <- item_id0[m_index]
  d <- length(user_id)   #総レビュー数
  
  #すべてのパターンが生成されればbreak
  if(length(unique(user_id))==hh & length(unique(item_id))==item) break
}

#単語数を設定
w <- rpois(d, rgamma(d, 25, 0.5))   #文書あたりの単語数
f <- sum(w)   #総単語数
n_user <- plyr::count(user_id)$freq
n_item <- plyr::count(item_id)$freq

#単語IDを設定
u_id <- rep(user_id, w)
i_id <- rep(item_id, w)
d_id <- rep(1:d, w)

#インデックスを設定
user_index <- list()
item_index <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
}

##パラメータの設定
#ディリクレ分布の事前分布の設定
alpha11 <- rep(0.15, k11)
alpha12 <- rep(0.15, k12)
alpha21 <- c(rep(0.5, v1/s), rep(0.025, v1/s), rep(0.0025, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha22 <- c(rep(0.3, v1/s), rep(0.1, v1/s), rep(0.025, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha23 <- c(rep(0.2, v1/s), rep(1.0, v1/s), rep(0.2, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha24 <- c(rep(0.025, v1/s), rep(0.1, v1/s), rep(0.3, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha25 <- c(rep(0.0025, v1/s), rep(0.025, v1/s), rep(0.5, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha2 <- rbind(alpha21, alpha22, alpha23, alpha24, alpha25)
alpha31 <- c(rep(0.001, v1/s), rep(0.001, v1/s), rep(0.001, v1/s), rep(0.1, v2), rep(0.002, v3))
alpha32 <- c(rep(0.001, v1/s), rep(0.001, v1/s), rep(0.001, v1/s), rep(0.002, v2), rep(0.1, v3))
beta1 <- c(1.6, 4.8, 5.6)

##すべての単語が出現するまでデータの生成を続ける
for(rp in 1:1000){
  print(rp)
  
  #事前分布からパラメータを生成
  theta11 <- thetat11 <- extraDistr::rdirichlet(hh, alpha11)
  theta12 <- thetat12 <- extraDistr::rdirichlet(item, alpha12)
  omega <- omegat <- extraDistr::rdirichlet(r, alpha2)
  phi <- phit <- extraDistr::rdirichlet(k11, alpha31)
  gamma <- gammat <- extraDistr::rdirichlet(k12, alpha32)
  lambda <- lambdat <- extraDistr::rdirichlet(hh, beta1)
  
  ##モデルに基づきデータを生成
  WX <- matrix(0, nrow=d, ncol=v)
  y <- rep(0, d)
  Z1_list <- Z21_list <- Z22_list <- wd_list <- list()
  
  for(i in 1:d){
    #ユーザーとアイテムを抽出
    u_index <- user_id[i]
    i_index <- item_id[i]
    
    #評価スコアを生成
    y[i] <- as.numeric(rmnom(1, 1, c(0.1, 0.225, 0.3, 0.25, 0.125)) %*% 1:r)

    #多項分布からスイッチング変数を生成
    z1 <- rmnom(w[i], 1, lambda[u_index, ])
    z1_vec <- as.numeric(z1 %*% 1:a)
    index_z11 <- which(z1[, 1]==1)
    
    #ユーザートピックを生成
    z21 <- matrix(0, nrow=w[i], ncol=k11)
    index_z21 <- which(z1[, 2]==1)
    if(sum(z1[, 2]) > 0){
      z21[index_z21, ] <- rmnom(sum(z1[, 2]), 1, theta11[u_index, ])
    }
    z21_vec <- as.numeric(z21 %*% 1:k11)
    
    #アイテムトピックを生成
    z22 <- matrix(0, nrow=w[i], ncol=k12)
    index_z22 <- which(z1[, 3]==1)
    if(sum(z1[, 3]) > 0){
      z22[index_z22, ] <- rmnom(sum(z1[, 3]), 1, theta12[i_index, ])
    }
    z22_vec <- as.numeric(z22 %*% 1:k12)
    
    #トピックから単語を生成
    words <- matrix(0, nrow=w[i], ncol=v)
    if(sum(z1[, 1]) > 0){
      words[index_z11, ] <- rmnom(sum(z1[, 1]), 1, omega[y[i], ])
    }
    if(sum(z1[, 2]) > 0){
      words[index_z21, ] <- rmnom(sum(z1[, 2]), 1, phi[z21_vec[index_z21], ])
    }
    if(sum(z1[, 3]) > 0){
      words[index_z22, ] <- rmnom(sum(z1[, 3]), 1, gamma[z22_vec[index_z22], ])
    }
    word_vec <- as.numeric(words %*% 1:v)
    WX[i, ] <- colSums(words)
    
    #データを格納
    wd_list[[i]] <- word_vec
    Z1_list[[i]] <- z1
    Z21_list[[i]] <- z21
    Z22_list[[i]] <- z22
  }
  if(min(colSums(WX)) > 0) break
}

#リストを変換
wd <- unlist(wd_list)
Z1 <- do.call(rbind, Z1_list)
Z21 <- do.call(rbind, Z21_list)
Z22 <- do.call(rbind, Z22_list)
storage.mode(Z1) <- "integer"
storage.mode(Z21) <- "integer"
storage.mode(Z22) <- "integer"
storage.mode(WX) <- "integer"


####マルコフ連鎖モンテカルロ法でSwitching Binomial LDAを推定####
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
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##事前分布の設定
alpha11 <- 0.1; alpha12 <- 0.1
alpha21 <- 0.1; alpha22 <- 0.1; alpha23 <- 0.1
beta <- 0.5

##パラメータの真値
theta11 <- thetat11
theta12 <- thetat12
phi <- phit
gamma <- gammat
omega <- omegat
lambda <- lambdat

##パラメータの初期値を設定
#トピック分布の初期値
theta11 <- extraDistr::rdirichlet(hh, rep(1.0, k11))
theta12 <- extraDistr::rdirichlet(item, rep(1.0, k12))

#単語分布の初期値
phi <- extraDistr::rdirichlet(k11, rep(1.0, v))   #ユーザーの単語分布の初期値
gamma <- extraDistr::rdirichlet(k12, rep(1.0, v))   #アイテムの単語分布の初期値
omega <- extraDistr::rdirichlet(r, rep(5.0, v))   #評価スコアの単語分布の初期値

#スイッチング変数の初期値
lambda <- matrix(1/s, nrow=hh, ncol=s)


##パラメータの格納用配列
THETA11 <- array(0, dim=c(hh, k11, R/keep))
THETA12 <- array(0, dim=c(item, k12, R/keep))
PHI <- array(0, dim=c(k11, v, R/keep))
GAMMA <- array(0, dim=c(k12, v, R/keep))
OMEGA <- array(0, dim=c(r, v, R/keep))
LAMBDA <- array(0, dim=c(hh, s, R/keep))
SEG1 <- matrix(0, nrow=f, ncol=a)
SEG21 <- matrix(0, nrow=f, ncol=k11)
SEG22 <- matrix(0, nrow=f, ncol=k12)
storage.mode(SEG21) <- "integer"
storage.mode(SEG22) <- "integer"

##データとインデックスの設定
#インデックスの設定
user_list <- user_vec <- list()
item_list <- item_vec <- list()
wd_list <- wd_vec <- list()
user_n <- rep(0, hh)
item_n <- rep(0, item)
for(i in 1:hh){
  user_list[[i]] <- which(u_id==i)
  user_vec[[i]] <- rep(1, length(user_list[[i]]))
  user_n[i] <- sum(user_vec[[i]])
}
for(i in 1:item){
  item_list[[i]] <- which(i_id==i)
  item_vec[[i]] <- rep(1, length(item_list[[i]]))
  item_n[i] <- sum(item_vec[[i]])
}
for(j in 1:v){
  wd_list[[j]] <- which(wd==j)
  wd_vec[[j]] <- rep(1, length(wd_list[[j]]))
}

#データの設定
y_vec <- y[d_id]
y_data <- matrix(as.numeric(table(1:f, y_vec)), nrow=f, ncol=r)
storage.mode(y_data) <- "integer"
r_vec <- rep(1, r)
a_vec <- rep(1, a)
vec11 <- rep(1, k11)
vec12 <- rep(1, k12)

##対数尤度の基準値
par <- colSums(WX) / f
LLst <- sum(WX %*% log(par))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){

  ##多項分布よりスイッチング変数をサンプリング
  #評価スコア、ユーザーおよびアイテムの期待尤度を設定
  Li_score <- as.numeric((t(omega)[wd, ] * y_data) %*% r_vec)   #スコア尤度
  Li_user <- theta11[u_id, ] * t(phi)[wd, ]   #ユーザー尤度
  par_user <- as.numeric(Li_user %*% vec11)   #ユーザーの期待尤度
  Li_item <- theta12[i_id, ] * t(gamma)[wd, ]   #アイテム尤度
  par_item <- as.numeric(Li_item %*% vec12)   #アイテムの期待尤度
  par <- cbind(Li_score, par_user, par_item)
  
  #潜在確率からスイッチング変数を生成
  lambda_r <- lambda[u_id, ]   #スイッチング変数の事前分布
  par_r <- lambda_r * par
  s_prob <- par_r / as.numeric(par_r %*% a_vec)   #スイッチング変数の割当確率
  Zi1 <- rmnom(f, 1, s_prob)   #多項分布からスイッチング変数を生成
  Zi1_T <- t(Zi1)
  index_z21 <- which(Zi1[, 2]==1)
  index_z22 <- which(Zi1[, 3]==1)
  
  #ディリクレ分布から混合率をサンプリング
  rsum0 <- matrix(0, nrow=hh, ncol=a)
  for(i in 1:hh){
    rsum0[i, ] <- Zi1_T[, user_list[[i]]] %*% user_vec[[i]]
  }
  rsum <- rsum0 + beta   #ディリクレ分布のパラメータ
  lambda <- extraDistr::rdirichlet(hh, rsum)   #ディリクレ分布からlambdaをサンプリング
  
  
  ##ユーザーおよびアイテムのトピックをサンプリング
  #トピックの割当確率を推定
  z_rate1 <- Li_user[index_z21, ] / par_user[index_z21]   #ユーザーのトピック割当確率
  z_rate2 <- Li_item[index_z22, ] / par_item[index_z22]   #アイテムのトピック割当確率
  
  #多項分布からトピックを生成
  Zi21 <- matrix(0, nrow=f, ncol=k11)
  Zi22 <- matrix(0, nrow=f, ncol=k12)
  Zi21[index_z21, ] <- rmnom(nrow(z_rate1), 1, z_rate1)
  Zi22[index_z22, ] <- rmnom(nrow(z_rate2), 1, z_rate2)
  Zi21_T <- t(Zi21)
  Zi22_T <- t(Zi22)
  
  
  ##トピックモデルのパラメータをサンプリング
  #ユーザーのトピック分布をサンプリング
  wusum0 <- matrix(0, nrow=hh, ncol=k11)
  for(i in 1:hh){
    wusum0[i, ] <- Zi21_T[, user_list[[i]]] %*% user_vec[[i]]
  }
  wusum <- wusum0 + alpha21   #ディリクレ分布のパラメータ
  theta11 <- extraDistr::rdirichlet(hh, wusum)   #ディリクレ分布からtheta21をサンプリング
  
  #アイテムのトピック分布をサンプリング
  wisum0 <- matrix(0, nrow=item, ncol=k12)
  for(i in 1:item){
    wisum0[i, ] <- Zi22_T[, item_list[[i]]] %*% item_vec[[i]]
  }
  wisum <- wisum0 + alpha22   #ディリクレ分布のパラメータ
  theta12 <- extraDistr::rdirichlet(item, wisum)   #ディリクレ分布からtheta22をサンプリング
  
  
  ##評価スコア、ユーザーおよびアイテムの単語分布をサンプリング
  y_data_t <- t(y_data * Zi1[, 1])
  vssum0 <- matrix(0, nrow=r, ncol=v)
  vusum0 <- matrix(0, nrow=k11, ncol=v)
  visum0 <- matrix(0, nrow=k12, ncol=v)
  for(j in 1:v){
    vssum0[, j] <- y_data_t[, wd_list[[j]], drop=FALSE] %*% wd_vec[[j]]
    vusum0[, j] <- Zi21_T[, wd_list[[j]], drop=FALSE] %*% wd_vec[[j]]
    visum0[, j] <- Zi22_T[, wd_list[[j]], drop=FALSE] %*% wd_vec[[j]]
  }
  vssum <- vssum0 + alpha21; vusum <- vusum0 + alpha22; visum <- visum0 + alpha23
  omega <- extraDistr::rdirichlet(r, vssum)
  phi <- extraDistr::rdirichlet(k11, vusum)
  gamma <- extraDistr::rdirichlet(k12, visum)
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA11[, , mkeep] <- theta11
    THETA12[, , mkeep] <- theta12
    PHI[, , mkeep] <- phi
    GAMMA[, , mkeep] <- gamma
    OMEGA[, , mkeep] <- omega
    LAMBDA[, , mkeep] <- lambda
  }  
  
  #トピック割当はバーンイン期間を超えたら格納する
  if(rp%%keep==0 & rp >= burnin){
    SEG1 <- SEG1 + Zi1
    SEG21 <- SEG21 + Zi21
    SEG22 <- SEG22 + Zi22
  }
  
  if(rp%%disp==0){
    #対数尤度を計算
    LL <- sum(log(rowSums(Zi1 * par)))
    
    #サンプリング結果を確認
    print(rp)
    print(c(LL, LLst))
    print(round(c(colMeans(Zi1), colMeans(Z1)), 3))
    print(round(cbind(phi[, (v1-5):(v1+4)], phit[, (v1-5):(v1+4)]), 3))
  }
}

####サンプリング結果の可視化と要約####
burnin <- 1000/keep
RS <- R/keep

##サンプリング結果のプロット
#スイッチング変数の混合率の可視化
matplot(t(LAMBDA[, 1, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(LAMBDA[, 2, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(LAMBDA[, 3, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")

#トピック分布のサンプリング結果をプロット
matplot(t(THETA11[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA11[100, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA11[250, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA11[500, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA11[1000, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA12[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA12[50, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA12[100, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA12[150, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA12[200, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")

#単語分布のサンプリング結果の可視化
matplot(t(PHI[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI[7, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI[9, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(GAMMA[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(GAMMA[4, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(GAMMA[8, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(GAMMA[12, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(GAMMA[15, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(OMEGA[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(OMEGA[2, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(OMEGA[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(OMEGA[4, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(OMEGA[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")


##サンプリング結果の要約
#サンプリング結果の事後平均
round(cbind(apply(LAMBDA[, , burnin:RS], c(1, 2), mean), lambdat), 3)   #スイッチング変数の混合率の事後平均
round(cbind(apply(THETA11[, , burnin:RS], c(1, 2), mean), thetat11), 3)   #ユーザーのトピック割当の事後平均
round(cbind(apply(THETA12[, , burnin:RS], c(1, 2), mean), thetat12), 3)   #アイテムのトピック割当の事後平均
round(cbind(t(apply(PHI[, , burnin:RS], c(1, 2), mean)), t(phit)), 3)   #ユーザーの単語分布の事後平均
round(cbind(t(apply(GAMMA[, , burnin:RS], c(1, 2), mean)), t(gammat)), 3)   #アイテムの単語分布の事後平均
round(cbind(t(apply(OMEGA[, , burnin:RS], c(1, 2), mean)), t(omegat)), 3)   #評価スコアの単語分布の事後平均


#サンプリング結果の事後信用区間
round(apply(LAMBDA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.025)), 3)
round(apply(LAMBDA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.975)), 3)
round(apply(THETA1[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.025)), 3)
round(apply(THETA2[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.975)), 3)
round(t(apply(PHI[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.025))), 3)
round(t(apply(PHI[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.975))), 3)
round(t(apply(GAMMA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.025))), 3)
round(t(apply(GAMMA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.975))), 3)
round(t(apply(OMEGA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.025))), 3)
round(t(apply(OMEGA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.975))), 3)


##サンプリングされた潜在変数の要約
n <- max(SEG1)
round(cbind(SEG1/n, Z1), 3)
round(cbind(rowSums(SEG21), SEG21/max(rowSums(SEG21)), Z21 %*% 1:k11), 3)
round(cbind(rowSums(SEG22), SEG22/max(rowSums(SEG22)), Z22 %*% 1:k12), 3)

