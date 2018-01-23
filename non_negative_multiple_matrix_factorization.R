#####非負値多重行列因子分解#####
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
hh <- 3000   #ユーザー数
item <- 300   #アイテム数
group <- 20   #グループ数
category <- 15   #カテゴリー数
k <- 10   #潜在変数数

#対応行列の設定
pr1 <- runif(group, 0.25, 1.0)
pr2 <- runif(category, 0.25, 1.0)
V <- rmnom(hh, 1, pr1)   #ユーザーグループ対応行列
W <- rmnom(item, 1, pr2)   #商品カテゴリ対応行列


##非負値行列因子分解の仮定に従いデータを生成
#ガンマ分布よりパラメータを設定
A0 <- matrix(0, nrow=hh, ncol=k)
B0 <- matrix(0, nrow=k, ncol=item)
alpha01 <- matrix(0, nrow=k, ncol=group)
beta01 <- matrix(0, nrow=k, ncol=group)
alpha02 <- matrix(0, nrow=k, ncol=category)
beta02 <- matrix(0, nrow=k, ncol=category)

for(i in 1:k){
  #ユーザー特徴行列を生成
  for(j1 in 1:group){
    alpha01[i, j1] <- runif(1, 0.01, 0.4)
    beta01[i, j1] <- runif(1, 0.6, 1.8)
    A0[V[, j1]==1, i] <- rgamma(sum(V[, j1]), alpha01[i, j1], beta01[i, j1])
  }
  #アイテム特徴行列を生成
  for(j2 in 1:category){
    alpha02[i, j2] <- runif(1, 0.01, 0.3)
    beta02[i, j2] <- runif(1, 0.6, 1.5)
    B0[i, W[, j2]==1] <- rgamma(sum(W[, j2]), alpha02[i, j2], beta02[i, j2])
  }
}

#ポアソン分布の期待値を設定
AB <- A0 %*% B0

#グループ特徴行列およびカテゴリー特徴行列を生成
C0 <- t(V) %*% A0 
D0 <- t(W) %*% t(B0)


##カルバックライブラーダイバージェンス距離に基づき応答変数を生成
#ユーザー商品購買行列を生成
X <- matrix(0, nrow=hh, ncol=item)
for(j in 1:item){
  X[, j] <- rpois(hh, AB[, j])
}
rowSums(X)
sum(X)

#ユーザーカテゴリ購買行列を設定
Y <- matrix(0, nrow=hh, ncol=category)
for(j in 1:category){
  Y[, j] <- rowSums(X[, W[, j]==1])
}
colSums(Y)

#グループ商品購買行列を設定
Z <- matrix(0, nrow=group, ncol=item)
for(j in 1:group){
  Z[j, ] <- colSums(X[V[, j]==1, ])
}


####マルコフ連鎖モンテカルロ法でMMNFを推定####
##アルゴリズムの設定
R <- 10000
keep <- 2
iter <- 0
disp <- 10

##事前分布の設定
alpha1 <- 0.01; beta1 <- 0.01
alpha2 <- 0.01; beta2 <- 0.01

##初期値の設定
A <- matrix(rgamma(hh*k, 0.1, 0.1), nrow=hh, ncol=k)
B <- matrix(rgamma(category*k, 0.1, 0.1), nrow=k, ncol=category)

##サンプリング結果の保存用配列
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, category, R/keep))
lambda <- array(0, dim=c(hh, category, k))





mean(rgamma(1000000, 0.5, 8))
mean(rgamma(1000000, 1, 16))
