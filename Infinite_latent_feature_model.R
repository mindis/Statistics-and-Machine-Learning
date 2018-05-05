#####無限因子分析モデル#####
options(warn=0)
library(MASS)
library(NMF)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(529867)

####データの発生####
##データの設定
k <- 10   #基底数
hh <- 3000   #ユーザー数
v <- 100   #項目数

##パラメータの設定
beta1 <- 17.5
beta2 <- 25.0
mu_vec <- rep(0, k)
cov <- diag(k)
sigma <- 1

##応答変数の生成
Z <- ZT <- matrix(rbinom(hh*k, 1, rbeta(hh*k, beta1, beta2)), nrow=hh, ncol=k)   #潜在変数行列
X <- XT<- t(mvrnorm(v, mu_vec, cov))   #特徴行列
Mu <- Z %*% X   #平均構造を生成
Y <- Mu + matrix(rnorm(hh*v, 0, sigma), nrow=hh, ncol=v)   #応答変数


####マルコフ連鎖モンテカルロ法で無限因子分析を推定####
##アルゴリズムの設定
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10
g <- 5

##事前分布の設定
sigma <- 1
mu_vec <- rep(0, k)
cov <- diag(k)
alpha <- 0.0001   #IBPの事前分布

##初期値の設定
k0 <- 1
Z <- matrix(rbinom(hh, 1, 0.4), nrow=hh, ncol=)
X <- t(rnorm(v, 0, 1))
Mu <- Z %*% X

##データの設定
max_seg <- 20
vec <- 1:v
z_vec <- rep(1, hh)


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##ベルヌーイ分布より既存潜在変数zを生成
  Z_new <- matrix(0, nrow=hh, ncol=k0+1)
  Z0 <- Z1 <- Z
  
  for(j in 1:k0){
    #対数尤度を設定
    Z1[, j] <- 1; Z0[, j] <- 0   #潜在変数を置き換え
    LLi0 <- as.numeric(dnorm(Y, Z0 %*% X, sigma, log=TRUE) %*% vec)
    LLi1 <- as.numeric(dnorm(Y, Z1 %*% X, sigma, log=TRUE) %*% vec)
    LLi <- cbind(LLi0, LLi1)
    Li <- exp(LLi - rowMaxs(LLi)) + 10^-100
    
    #潜在変数の割当確率からzを生成
    gamma0 <- sum(Z[, j]) / hh   #事前分布
    gamma_par <- cbind((1-gamma0) * Li[, 1], gamma0 * Li[, 2])
    gamma <- gamma_par[, 2] / rowSums(gamma_par)   #潜在変数の割当確率
    Z_new[, j] <- rbinom(hh, 1, gamma)   #ベルヌーイ分布から潜在変数を生成
  }
  Z <- Z_new
  
  ##IBPから新しい潜在変数zを生成
  #対数尤度を設定
  X_new <- t(rnorm(v, 0, sigma))
  lambda <- alpha/hh
  LLi0 <- as.numeric(dnorm(Y, Z[, 1:k0] %*% X, sigma, log=TRUE) %*% vec)
  LLi1 <- as.numeric(dnorm(Y, Z[, 1:k0] %*% X + z_vec %*% X_new, sigma, log=TRUE) %*% vec) + dpois(k0+1, lambda, log=TRUE)
  LLi <- cbind(LLi0, LLi1)
  Li <- exp(LLi - rowMaxs(LLi)) 

  #潜在変数の割当確率からzを生成
  gamma_par <- cbind((1-alpha/hh) * Li[, 1], alpha/hh * Li[, 2])
  gamma <- gamma_par[, 2] / rowSums(gamma_par)   #潜在変数の割当確率
  Z[, (k0+1)] <- rbinom(hh, 1, gamma)   #ベルヌーイ分布から潜在変数を生成
  
  #新しい潜在変数が一定数以下なら潜在変数を削除
  if(sum(Z[, (k0+1)]) <= g){
    Z <- Z[, -(k0+1)]
  }
  k0 <- ncol(Z)   #基底数を更新
  
  lambda <- alpha/hh
  dpois(10, 10-lambda, log=TRUE)
  
  
  ##特徴行列を更新
  #多変量回帰モデルの平均ベクトル
  ZZ <- t(Z) %*% Z + diag(k0)
  inv_ZZ <- solve(ZZ)
  beta_mu <- inv_ZZ %*% t(Z) %*% Y
  
  #多変量正規分布より特徴行列をサンプリング
  X <- matrix(0, nrow=k0, ncol=v)
  for(j in 1:v){
    X[, j] <- mvrnorm(1, beta_mu[, j], inv_ZZ)
  }
  
  LL <- sum(dnorm(Y, Z %*% X, sigma, log=TRUE))
  print(LL)
  print(k0)
  print(colSums(Z))
}


sum(dnorm(Y, ZT %*% XT, sigma, log=TRUE))



