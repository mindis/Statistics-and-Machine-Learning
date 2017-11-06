#####階層ベイズ正則化多項ロジットモデル#####
library(MASS)
library(mclust)
library(reshape2)
library(bayesm)
detach("package:bayesm", unload=TRUE)
library(ExtDist)
library(extraDistr)
library(matrixStats)
library(glmnet)
library(monomvn)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(534798)

####データの発生####
##データの設定
hh <- 2000
pt0 <- rpois(hh, 4.4)
pt <- ifelse(pt0==0, 1, pt0)
hhpt <- sum(pt)
select <- 8

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
k <- 100   #特徴量数
cont <- 70   #連続変数数
bin <- 30   #二値変数数

p <- runif(1, 0.4, 0.6)
x.cont <- matrix(rnorm(hh*cont, 0, 1), nrow=hh, ncol=cont)
x.bin <- matrix(rbinom(hh*bin, 1, p), nrow=hh, ncol=bin)
Data <- cbind(1, x.cont, x.bin)


####応答変数の発生#### 
##切片をベクトル変換
X <- matrix(diag(select), nrow=hhpt*select, ncol=select, byrow=T)[, -select]

for(i in 1:1000){
  print(i) 
  
  ##パラメータの設定
  #階層モデルの回帰パラメータを設定
  b01 <- runif(select-1, -0.8, 0.8)
  b02 <- matrix(runif(cont*(select-1), 0, 1.1), nrow=cont, ncol=select-1) * 
    matrix(rbinom(cont*(select-1), 1, 0.3), nrow=cont, ncol=select-1)
  b03 <- matrix(runif(bin*(select-1), -1.1, 1.3), nrow=bin, ncol=select-1) * 
    matrix(rbinom(bin*(select-1), 1, 0.3), nrow=bin, ncol=select-1)
  b0 <- rbind(b01, b02, b03)
  
  #階層モデルの分散を設定
  cov0 <- 0.35
  
  ##応答変数の発生
  logit <- matrix(0, nrow=hhpt, ncol=select)
  tau0 <- matrix(0, nrow=hh, ncol=select-1)
  
  #ロジットの定義
  for(j in 1:(select-1)){
    tau0[, j] <- rnorm(hh, 0, cov0)
    logit[, j] <- (Data %*% b0[, j] + tau0[, j])[ID$id, ]
  }
  
  #多項分布から応答変数を発生
  Pr0 <- exp(logit) / rowSums(exp(logit))
  y <- t(apply(Pr0, 1, function(x) rmultinom(1, 1, x)))
  
  if(max(colMeans(y)) < 0.4 & min(colMeans(y)) > 0.05) break
}

#発生させたデータを確認
colSums(y); round(colMeans(Pr0), 3)
rownames(b0) <- NULL

####マルコフ連鎖モンテカルロ法で階層ベイズ正則化多項ロジットモデルを推定####
##多項ロジットモデルの対数尤度関数を設定
fr <- function(y, X, beta, select, hhpt){
  
  #ロジットと確率の計算
  logit <- matrix(rowSums(X * beta), nrow=hhpt, ncol=select, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hhpt, ncol=select)
  
  #対数尤度を定義
  LLi <- rowSums(y * log(Pr))
  return(LLi)
}


##アルゴリズムの設定
R <- 20000
keep <- 4
sbeta <- 1.5
iter <- 0

##初期値の設定
beta0 <- scale(colSums(y))
oldbeta <- mvrnorm(hh, beta0[-select], diag(0.2, select-1))
oldtheta <- matrix(0, nrow=k+1, ncol=select-1)
oldcov <- diag(0.1, select-1)
inv_cov <- solve(oldcov)

##サンプリング結果の保存用配列
THETA <- array(0, dim=c(k+1, select-1, R/keep))
BETA <- array(0, dim=c(hh, select-1, R/keep))
TAU <- array(0, dim=c(select-1, select-1, R/keep))

##インデックスの設定
index_id <- list()
index_y <- list()
for(i in 1:hh){
  index_id[[i]] <- which(ID_v$id==i)
  index_y[[i]] <- which(ID$id==i)
}
lognew <- rep(0, hh)
logold <- rep(0, hh)
logpnew <- rep(0, hh)
logpold <- rep(0, hh)
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
  
  ##多項ロジットモデルのID別パラメータをサンプリング
  #新しいパラメータをサンプリング
  betad <- oldbeta
  betan <- betad + mvrnorm(hh, rep(0, select-1), diag(0.025 ,select-1))
  
  #誤差を設定
  mu <- Data %*% oldtheta
  er_new <- betan - mu
  er_old <- betad - mu
  
  #階層モデルの分散共分散行列を推定
  oldcov <- var(er_old)
  inv_cov <- solve(oldcov)
  
  #対数尤度と対数事前分布を計算
  lognew0 <- fr(y, X, betan[id_vec, ], select, hhpt)
  logold0 <- fr(y, X, betad[id_vec, ], select, hhpt)
  
  for(i in 1:hh){
    lognew[i] <- sum(lognew0[index_y[[i]]])
    logold[i] <- sum(logold0[index_y[[i]]])
    logpnew[i] <- -0.5 * er_new[i, ] %*% inv_cov %*% er_new[i, ]
    logpold[i] <- -0.5 * er_old[i, ] %*% inv_cov %*% er_old[i, ]
  }
  
  #MHサンプリングでパラメータの採択を決定
  #メトロポリスヘイスティング法でパラメータの採択を決定
  rand <- runif(hh)   #一様分布から乱数を発生
  LLind_diff <- exp(lognew + logpnew - logold - logpold)   #採択率を計算
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=hh, ncol=select-1)
  oldbeta <- flag*betan + (1-flag)*betad   #alphaがrandを上回っていたら採択
  
  
  ##ベイジアンlassoで階層モデルの回帰パラメータをサンプリング
  for(j in 1:(select-1)){
    res <- blasso(X=Data[, -1], y=oldbeta[, j], beta=oldtheta[-1, j], lambda2=lambda[j], s2=diag(oldcov)[j], T=2)
    oldtheta[, j] <- c(res$mu[2], res$beta[2, ])
    lambda[j] <- res$lambda2[2]
  }
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- oldtheta
    BETA[, , mkeep] <- oldbeta
    TAU[, , mkeep] <- oldcov
    print(rp)
    print(sum(lognew))
    print(round(lambda, 3))
    print(round(t(cbind(oldtheta, b0)[1:20, ]), 3))
  }
}

matplot(t(THETA[7, , ]), type="l")



