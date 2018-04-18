#####mixed effect poisson block model#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(flexmix)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(506832)

####データの発生####
##データの設定
d <- 150   #アイテム数
k <- 7   #潜在変数数
N <- d*(d-1)/2   #総サンプル数
vec <- rep(1, k)

##IDと潜在変数の設定
#潜在変数を生成
Z <- rmnom(d, 1, runif(k, 1, 3))
z <- as.numeric(Z %*% 1:k)

#IDを設定
id1 <- id2 <- c()
for(i in 1:(d-1)){
  id1 <- c(id1, rep(i, length((i+1):d)))
  id2 <- c(id2, (i+1):d)
}

##応答変数の生成
#パラメータを生成
cov <- covt <- 0.75
beta <- betat <- 0.8
alphat <- alpha <- rnorm(d, beta, cov)   #変量効果のパラメータ
theta0 <- matrix(rnorm(k*k, 0, 0.85), nrow=k, ncol=k)   #潜在変数のパラメータ
theta0[upper.tri(theta0)] <- 0
theta <- theta0 + t(theta0)
diag(theta) <- diag(theta0)
thetat <- theta

#ポアソン分布の平均構造
mu <- alpha[id1] + alpha[id2] + (theta[z[id1], ] * Z[id2, ]) %*% vec
lambda <- exp(mu)


#ポアソン分布から応答変数を生成
y <- rpois(N, lambda)
sum(y); mean(y)
hist(y, xlab="頻度", main="アイテム間の出現頻度", col="grey", breaks=25)


####マルコフ連鎖モンテカルロ法でmixed effect poisson block modelを推定####
##変量効果ポアソン回帰モデルの対数尤度
loglike <- function(alpha, theta, y, y_factorial, z, Z, vec, id1, id2){
  
  #尤度を定義する
  lambda <- exp(alpha[id1] + alpha[id2] + (theta[z[id1], ] * Z[id2, ]) %*% vec)   #平均構造
  LLi <- as.numeric(y*log(lambda)-lambda - lfactorial(y))   #対数尤度
  LL <- sum(LLi)   #対数尤度の和
  
  #結果を返す
  LL_value <- list(LLi=LLi, LL=LL)
  return(LL_value)
}


##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##事前分布の設定
#回帰係数の事前分布
alpha01 <- 0
beta01 <- 0
tau01 <- 0.01

#変量効果の事前分布
alpha02 <- 0
s02 <- 0.01
v02 <- 0.01

##パラメータの真値
theta <- thetat
alpha <- alphat
beta <- betat
cov <- covt

##初期値の設定
#変量効果の初期値
cov <- 0.5   #階層モデルの標準偏差の初期値
beta <- 0.8   #階層モデルの平均の初期値
mu <- rep(0, d)
for(i in 1:d){
  mu[i] <- mean(y[which(id1==i | id2==i)])
}
x <- rnorm(d, 0, 0.5)
rank_mu <- ceiling(rank(mu))
alpha <- sort(x, decreasing=TRUE)[rank_mu]   #変量効果の初期値

#潜在変数のパラメータ
theta0 <- matrix(rnorm(k*k, 0, 0.3), nrow=k, ncol=k)   #潜在変数のパラメータ
theta0[upper.tri(theta0)] <- 0
theta <- theta0 + t(theta0)

#潜在変数の初期値
Zi <- rmnom(d, 1, rep(1/k, k))
z_vec <- as.numeric(Zi %*% 1:k)


##定数を計算
y_factorial <- lfactorial(y)
Y_factorial <- matrix(y_factorial, nrow=N, ncol=k*k)
upper_tri <- matrix(as.logical(upper.tri(theta) + diag(1, k)), nrow=k, ncol=k)
vec <- rep(1, k)

##インデックスを設定
item_list <- list()
for(i in 1:d){
  item_list[[i]] <- which(id1==i | id2==i)
}


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
##メトロポリスヘイスティング法で変量効果をサンプリング
#新しいパラメータをサンプリング
alphad <- alpha
alphan <- alphad + rnorm(d, 0, 0.2)

#事前分布の誤差
er_new <- alphan - beta 
er_old <- alphad - beta

#対数尤度と対数事前分布を設定
lognew0 <- loglike(alpha, theta, y, y_factorial, z_vec, Zi, vec, id1, id2)$LLi
logold0 <- loglike(alpha, theta, y, y_factorial, z_vec, Zi, vec, id1, id2)$LLi
logpnew <- -0.5 * (er_new^2 / cov)
logpold <- -0.5 * (er_old^2 / cov)
 
#アイテムごとに対数尤度の和を取る
lognew <- logold <- rep(0, d)
for(i in 1:d){
  lognew[i] <- sum(lognew0[item_list[[i]]])
  logold[i] <- sum(logold0[item_list[[i]]])
}

#MHサンプリング
rand <- runif(d)   #一様分布から乱数を発生
LLind_diff <- exp(lognew + logpnew - logold - logpold)   #採択率を計算
alpha2 <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff

#alphaの値に基づき新しいbetaを採択するかどうかを決定
flag <- (alpha2 >= rand)*1 + (alpha2 < rand)*0
alpha <- flag*alphan + (1-flag)*alphad   #alphaがrandを上回っていたら採択


##ギブスサンプリングで潜在変数をサンプリング

item_list[[1]]
theta
item_list[[1]] 


alpha[id1] + alpha[id2]
item_list[[1]]
id2


Z_mu <- matrix(theta[(upper.tri(theta) + diag(k))==1], nrow=N, ncol=k*(k-1)/2 + k, byrow=T)
lambda <- exp(matrix(alpha[id1] + alpha[id2], nrow=N, ncol=k*(k-1)/2 + k) + Z_mu)
LLi <- y*log(lambda)-lambda - y_factorial   #対数尤度
LLi
i <- 31
z_par <- exp(LLi - rowMaxs(LLi))
round(z_par / rowSums(z_par), 3)[i, ]
y[i]
theta[(upper.tri(theta) + diag(k))==1]  


rmnom(1, 1, z_par / sum(z_par)) %*% z_index
(Z %*% 1:k)[i, ]


z_index <- rep(1:k, rep(k, k))

