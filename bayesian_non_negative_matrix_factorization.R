#####ベイジアン非負値行列因子分解#####
library(MASS)
library(matrixStats)
library(FAdist)
library(NMF)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

####データの発生####
#データの設定
hh <- 5000   #ユーザー数
category <- 200   #カテゴリー数
k <- 10   #潜在変数数

##非負値行列因子分解の仮定に従いデータを生成
#ガンマ分布よりパラメータを設定
alpha01 <- 0.2; beta01 <- 1.0
alpha02 <- 0.15; beta02 <- 0.8
W0 <- matrix(rgamma(hh*k, alpha01, beta01), nrow=hh, ncol=k)
H0 <- matrix(rgamma(category*k, alpha02, beta02), nrow=k, ncol=category)
WH <- W0 %*% H0

#ポアソン分布よりデータを生成
Data <- matrix(0, nrow=hh, ncol=category)
for(j in 1:category){
  Data[, j] <- rpois(hh, WH[, j])
}
colSums(Data)
rowSums(Data)
LL <- sum(dpois(Data, W0 %*% H0, log=TRUE))


####マルコフ連鎖モンテカルロ法でNMFを推定####
##アルゴリズムの設定
R <- 10000
keep <- 4
iter <- 0

##事前分布の設定
alpha1 <- 0.01; beta1 <- 0.01
alpha2 <- 0.01; beta2 <- 0.01

##初期値の設定
W <- matrix(rgamma(hh*k, 0.1, 0.1), nrow=hh, ncol=k)
H <- matrix(rgamma(category*k, 0.1, 0.1), nrow=k, ncol=category)

##サンプリング結果の保存用配列
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, category, R/keep))
lambda <- array(0, dim=c(hh, category, k))



####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##補助変数lambdaを更新
  lambda <- array(0, dim=c(hh, category, k))
  WH <- W %*% H
  for(j in 1:k){
    lambda[, , j] <- W[, j] %*% t(H[j, ]) / WH
  }
  
  ##ガンマ分布よりWをサンプリング
  for(j in 1:k){
    w1 <- alpha1 + rowSums(lambda[, , j] * Data)
    w2 <- beta1 + sum(H[j, ])
    W[, j] <- rgamma(hh, w1, w2)   
  }

  ##ガンマ分布よりHをサンプリング
  for(j in 1:k){
    h1 <- alpha2 + colSums(lambda[, , j] * Data)
    h2 <- beta2 + sum(W[, j])
    H[j, ] <- rgamma(category, h1, h2)  
  }


  ##サンプリング結果の保存と表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    W_array[, , mkeep] <- W
    H_array[, , mkeep ] <- H
    print(rp)
    print(c(sum(dpois(Data, W %*% H, log=T)), LL))
    print(round(cbind(W[1:10, 1:7], W0[1:10, 1:7]), 3))
    print(round(cbind(H[, 1:7], H0[, 1:7]), 3))
    
  }
}


      