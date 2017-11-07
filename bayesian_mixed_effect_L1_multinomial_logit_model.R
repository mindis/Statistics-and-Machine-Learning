#####変量効果正則化多項ロジットモデル#####
library(MASS)
library(mclust)
library(reshape2)
library(bayesm)
library(ExtDist)
library(extraDistr)
library(matrixStats)
library(glmnet)
library(monomvn)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)
gc(); gc()

#set.seed(534798)

####データの発生####
##データの設定
hh <- 2000
pt <- ceiling(rgamma(hh, 3.7, 0.9))
hhpt <- sum(pt)
select <- 10

#IDの設定
id <- rep(1:hh, pt)
time <- c()
for(i in 1:hh){
  time <- c(time, 1:pt[i])
}
ID <- data.frame(no=1:hhpt, id, time)

#IDのベクトル化
id_v <- as.numeric(t(matrix(id, nrow=hhpt, ncol=select)))
time_v <- as.numeric(t(matrix(time, nrow=hhpt, ncol=select)))
ID_v <- data.frame(no=1:length(id_v), id=id_v, time=time_v)


####説明変数の発生####
k <- 40   #特徴量数
data1 <- extraDistr::rdirichlet(hhpt, rep(0.2, k/2))[, -k/2]
data2 <- extraDistr::rdirichlet(hhpt, rep(0.2, k/2))[, -k/2]
Data <- cbind(1, data1, data2)

####応答変数の発生#### 
##切片をベクトル変換
X <- matrix(diag(select), nrow=hhpt*select, ncol=select, byrow=T)[, -select]

for(i in 1:1000){
  print(i) 
  
  ##パラメータの設定
  #lassoモデルの回帰パラメータを設定
  b01 <- runif(select-1, -0.8, 0.8)
  b02 <- matrix(runif((k-2)*(select-1), -3.25, 4.0), nrow=k-2, ncol=select-1) * 
              matrix(rbinom((k-2)*(select-1), 1, 0.3), nrow=k-2, ncol=select-1)
  b0 <- rbind(b01, b02)
  
  #分散を設定
  sigma0 <- 0.3
  
  ##変量効果のパラメータを設定
  cov0 <- diag(0.25, select-1)
  theta0 <- mvrnorm(hh, rep(0, select-1), diag(0.25, select-1))  
  
  ##応答変数の発生
  logit <- matrix(0, nrow=hhpt, ncol=select)
  tau0 <- matrix(0, nrow=hhpt, ncol=select-1)
  
  #ロジットの定義
  for(j in 1:(select-1)){
    tau0[, j] <- rnorm(hhpt, 0, sigma0)
    logit[, j] <- (Data %*% b0[, j] + tau0[, j]) + theta0[ID$id, j]
  }
  
  #多項分布から応答変数を発生
  Pr0 <- exp(logit) / rowSums(exp(logit))
  y <- t(apply(Pr0, 1, function(x) rmultinom(1, 1, x)))
  
  if(max(colMeans(y)) < 0.45 & min(colMeans(y)) > 0.05) break
}

#発生させたデータを確認
colSums(y); round(colMeans(Pr0), 3)
rownames(b0) <- NULL
sum(b0!=0)

####マルコフ連鎖モンテカルロ法で階層ベイズ正則化多項ロジットモデルを推定####
##多項ロジットモデルの対数尤度関数を設定
fr <- function(y, X, beta, gamma, select, hhpt){
  
  #ロジットと確率の計算
  logit <- cbind(beta + gamma, 0)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hhpt, ncol=select)
  
  #対数尤度を定義
  LLi <- rowSums(y * log(Pr))
  return(LLi)
}


##アルゴリズムの設定
R <- 40000
keep <- 4
sbeta <- 1.5
iter <- 0

##事前分布の設定
nu <- select
V <- solve(nu * diag(select-1))


##初期値の設定
#ベイジアンlasso回帰の初期値
beta0 <- scale(colSums(y))
oldgamma <- mvrnorm(hhpt, beta0[-select], diag(0.2, select-1))
oldtheta <- matrix(0, nrow=k-2+1, ncol=select-1)
oldsigma <- diag(0.1, select-1)
inv_sigma <- solve(oldsigma)

#変量効果の初期値
oldcov <- diag(0.1, select-1)
inv_cov <- solve(oldcov)
oldbeta <- mvrnorm(hh, rep(0, select-1), oldcov)

##サンプリング結果の保存用配列
THETA <- array(0, dim=c(k-2+1, select-1, R/keep))
BETA <- array(0, dim=c(hh, select-1, R/keep))
COV <- array(0, dim=c(select-1, select-1, R/keep))
SIGMA <- matrix(0, nrow=R/keep, ncol=select-1)

##インデックスの設定
index_id <- list()
index_y <- list()
for(i in 1:hh){
  index_id[[i]] <- which(ID_v$id==i)
  index_y[[i]] <- which(ID$id==i)
}

lognew1 <- rep(0, hh)
logold1 <- rep(0, hh)
logpnew1 <- rep(0, hh)
logpold1 <- rep(0, hh)
lambda <- rep(1, select-1)
er_new <- matrix(0, nrow=hh, select-1)
er_old <- matrix(0, nrow=hh, select-1)

#IDの設定
id_vec <- c()
for(i in 1:hh){
  id_vec <- c(id_vec, rep(i, sum(ID$id==i)*select))
}


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){

  ##変量効果をサンプリング
  #新しいパラメータをサンプリング
  betad <- oldbeta 
  betan <- betad + mvrnorm(hh, rep(0, select-1), diag(0.025, select-1))
  
  #対数尤度と対数事前分布を計算
  lognew0 <- fr(y, X, betan[ID$id, ], oldgamma, select, hhpt)
  logold0 <- fr(y, X, betad[ID$id, ], oldgamma, select, hhpt)
  logpnew1 <- -0.5 * rowSums(betan %*% inv_cov * betan)
  logpold1 <- -0.5 * rowSums(betad %*% inv_cov * betad)
  
  for(i in 1:hh){
    lognew1[i] <- sum(lognew0[index_y[[i]]])
    logold1[i] <- sum(logold0[index_y[[i]]])
  }
  
  #メトロポリスヘイスティング法でパラメータの採択を決定
  rand <- runif(hh)   #一様分布から乱数を発生
  LLind_diff <- exp(lognew1 + logpnew1 - logold1 - logpold1)   #採択率を計算
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=hh, ncol=select-1)
  oldbeta <- flag*betan + (1-flag)*betad   #alphaがrandを上回っていたら採択
  
  ##逆ウィシャート分布から階層モデルの分散共分散行列をサンプリング
  #逆ウィシャート分布のパラメータ
  R_par <- V + t(oldbeta) %*% oldbeta
  Sn <- nu + hh
  
  #逆ウィシャート分布から分散共分散行列をサンプリング
  oldcov <- rwishart(Sn, solve(R_par))$IW
  inv_cov <- solve(oldcov)
  
  
  ##多項ロジットモデルのサンプルごとのパラメータをサンプリング
  #新しいパラメータをサンプリング
  gammad <- oldgamma
  gamman <- gammad + mvrnorm(hhpt, rep(0, select-1), diag(0.05 ,select-1))
  
  #誤差を設定
  mu <- Data %*% oldtheta
  er_new <- gamman - mu
  er_old <- gammad - mu
  
  #lassoの分散共分散行列を推定
  oldsigma <- diag(diag(var(er_old)))
  inv_sigma <- solve(oldsigma)
  
  
  #対数尤度と対数事前分布を計算
  lognew2 <- fr(y, X, oldbeta[ID$id, ], gamman, select, hhpt)
  logold2 <- fr(y, X, oldbeta[ID$id, ], gammad, select, hhpt)
  logpnew2 <- -0.5 * rowSums(er_new %*% inv_sigma * er_new)
  logpold2 <- -0.5 * rowSums(er_old %*% inv_sigma * er_old)
  
  #メトロポリスヘイスティング法でパラメータの採択を決定
  rand <- runif(hhpt)   #一様分布から乱数を発生
  LLind_diff <- exp(lognew2 + logpnew2 - logold2 - logpold2)   #採択率を計算
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff

  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=hhpt, ncol=select-1)
  oldgamma <- flag*gamman + (1-flag)*gammad   #alphaがrandを上回っていたら採択
  
  
  ##ベイジアンlassoで階層モデルの回帰パラメータをサンプリング
  for(j in 1:(select-1)){
    res <- blasso(X=Data[, -1], y=oldgamma[, j], beta=oldtheta[-1, j], lambda2=lambda[j], s2=diag(oldsigma)[j], T=2)
    oldtheta[, j] <- c(res$mu[2], res$beta[2, ])
    lambda[j] <- res$lambda2[2]
  }

  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- oldtheta
    BETA[, , mkeep] <- oldbeta
    COV[, , mkeep] <- oldcov
    SIGMA[mkeep, ] <- diag(oldsigma)
    print(rp)
    print(sum(lognew1))
    print(round(lambda, 3))
    print(round(t(cbind(oldtheta, b0)[1:15, ]), 3))
    print(round(rbind(diag(oldcov), diag(cov0)), 3))
  }
}

matplot(t(THETA[4, , ]), type="l")
matplot(t(BETA[10, , ]), type="l")



