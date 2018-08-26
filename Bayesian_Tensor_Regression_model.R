#####Bayesian Tensor Regression#####
options(warn=0)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(HMM)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####データの発生####
##データの設定
r <- 5   #基底数
d <- 3   #テンソル数
hh <- 5000   #ユーザー数
w <- rpois(hh, (rgamma(hh, 10, 0.15)))   #ユーザーごとのサンプル数
f <- sum(w)   #総サンプル数

##素性ベクトルを生成
m1 <- 2; m2 <- 3; m3 <- 4
z1 <- matrix(runif(f*m1, 0, 1), nrow=f, ncol=m1)
z2 <- matrix(0, nrow=f, ncol=m2)
for(j in 1:m2){
  pr <- runif(1, 0.25, 0.55)
  z2[, j] <- rbinom(f, 1, pr)
}
z3 <- rmnom(f, 1, runif(m3, 0.2, 1.25)); z3 <- z3[, -which.min(colSums(z3))]
z <- cbind(z1, z2, z3)   #データを結合
m <- ncol(z)

##階層モデルの説明変数
m1 <- 1; m2 <- 3; m3 <- 5
u1 <- matrix(runif(hh*m1, 0, 1), nrow=hh, ncol=m1)
u2 <- matrix(0, nrow=hh, ncol=m2)
for(j in 1:m2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #データを結合

##テンソルを生成
k1 <- 20; k2 <- 7; k3 <- 5
X_list <- id_list <- no_list <- list()
N <- 0

for(i in 1:hh){
  #パラメータを生成
  lambda <- runif(k1, 0.25, 3.0)
  par1 <- extraDistr::rdirichlet(k2, rep(0.5, k2))
  par2 <- extraDistr::rdirichlet(k3, rep(0.5, k3))
  
  #データを生成
  n <- w[i]*k1
  x <- array(rpois(n, rep(lambda, w[i])) * rmnom(n, 1, par1), dim=c(w[i]*k1, k2, k3))   #1、2階層目のデータを生成
  z <- rmnom(n, 1, par2)   #3階層目のデータを生成
  for(j in 1:k3){
    x[, , j] <- x[, , j] * z[, j]   #テンソルを生成
  }
  
  #IDを生成
  id <- rep(1:w[i], rep(k1, w[i]))
  no <- rep(1:k1, w[i])
  
  #データを格納
  storage.mode(x) <- "integer"
  X_list[[i]] <- x
  id_list[[i]] <- id
  no_list[[i]] <- no
  N <- N + dim(x)[1]
}


##パラメータを生成
#素性ベクトルのパラメータ
beta <- c(5.0, rnorm(m-1, 0, 0.5))


##階層モデルとテンソルのパラメータを生成
#1階層目のモデルのパラメータ
alpha1 <- array(0, dim=c(ncol(u), r, k1))
Cov1 <- array(0, dim=c(r, r, k1))
theta1 <- array(0, dim=c(k1, r, hh))
for(j in 1:k1){
  alpha1[, , j] <- mvrnorm(ncol(u), rep(0, r), diag(runif(r, 0.005, 0.05)))
  Cov1[, , j] <- diag(runif(r, 0.01, 0.1))
  theta1[j, , ] <- t(u %*% alpha1[, , j] + mvrnorm(hh, rep(0, r), Cov1[, , j]))
}

#2階層目のモデルのパラメータ
alpha2 <- array(0, dim=c(ncol(u), r, k2))
Cov2 <- array(0, dim=c(r, r, k2))
theta2 <- array(0, dim=c(k2, r, hh))
for(j in 2:k2){
  alpha2[, , j] <- mvrnorm(ncol(u), rep(0, r), diag(runif(r, 0.005, 0.05)))
  Cov2[, , j] <- diag(runif(r, 0.01, 0.1))
  theta2[j, , ] <- t(u %*% alpha2[, , j] + mvrnorm(hh, rep(0, r), Cov2[, , j]))
}

#3階層目のモデルのパラメータ
alpha3 <- array(0, dim=c(ncol(u), r, k3))
Cov3 <- array(0, dim=c(r, r, k3))
theta3 <- array(0, dim=c(k3, r, hh))
for(j in 2:k3){
  alpha3[, , j] <- mvrnorm(ncol(u), rep(0, r), diag(runif(r, 0.005, 0.05)))
  Cov3[, , j] <- diag(runif(r, 0.01, 0.1))
  theta3[j, , ] <- t(u %*% alpha3[, , j] + mvrnorm(hh, rep(0, r), Cov3[, , j]))
}


