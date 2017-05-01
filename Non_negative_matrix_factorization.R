#####非負値行列因子分解#####
library(MASS)
library(NMF)
library(reshape2)
library(plyr)
library(lattice)
library(ggplot2)

####データの発生####
seg <- 30   #セグメント数
n <- 50   #セグメントごとのサンプル数
N <- seg*n   #全サンプル数
category <- 100   #カテゴリー数
lambda <- 50   #購買数のポアソン平均
y <- rpois(N, lambda)   #購買発生数
Y <- data.frame(seg=rep(1:seg, rep(n, seg)), y)

#購買頻度行列を発生させる
X.freq <- matrix(0, 0, ncol=category)
for(i in 1:seg){
  p <- runif(100, 0, 5)
  freq <- t(apply(Y[Y$seg==i, ], 1, function(x) rmultinom(1, x[2], p)))
  X.freq <- rbind(X.freq, freq)
}

####非負値行列因子分解を推定####
#パラメータの初期値を設定
H1 <- matrix(runif(N*seg, 0, 1), nrow=N, ncol=seg)
U1 <- matrix(runif(N*seg, 0, 1), nrow=seg, ncol=category) 

##パラメータを更新(二乗誤差基準)
X <- H1 %*% U1
for(i in 1:100){
  U <- t(t(U1) * (t(X.freq) %*% H1 / (t(X) %*% H1)))
  H <- H1 * (X.freq %*% t(U1) / X %*% t(U1))
  U1 <- U
  H1 <- H
}

t(U1) * t(X.freq) %*% H1 / (t(X) %*% H1)
H1 * (X.freq %*% t(U1) / X %*% t(U1))

t(X.freq) %*% rep(1, 1500)
X
H1
