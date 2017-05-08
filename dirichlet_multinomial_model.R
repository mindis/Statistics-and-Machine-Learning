#####ディクレリ多項分布モデル#####
library(MASS)
library(gtools)
library(reshape2)
library(plyr)
library(ggplot2)

####データの発生####
#set.seed(3479)
#データの設定
N1 <- 1500   #クラス1のサンプル数
N2 <- 1500   #クラス2のサンプル数
N <- N1 + N2   #総サンプル数
col <- 50   #変数数
n <- round(rpois(N, 200), 0)   #サンプルあたりの頻度

##ディクレリ多項分布からデータを発生させる
##クラス1のデータ発生
#ディクレリ分布からの乱数発生
alpha1 <- runif(col, 0.1, 1.5)   #パラメータ
theta1 <- rdirichlet(N1, alpha1)   #ディレクリ乱数を発生
round(theta1[, 1:15], 3)   #データを確認

#多項分布からの乱数発生
Z1 <- cbind(n[1:N1], theta1)
W1 <- t(apply(Z1, 1, function(x) rmultinom(1, x[1], x[-1])))
W1[, 1:20]

##クラス2のデータ発生
#ディクレリ分布からの乱数発生
alpha2 <- runif(col, 0.5, 2.0)   #パラメータ
theta2 <- rdirichlet(N2, alpha2)   #ディレクリ乱数を発生
round(theta2[, 1:15], 3)   #データを確認

#多項分布からの乱数発生
Z2 <- cbind(n[(N1+1):N], theta2)
W2 <- t(apply(Z2, 1, function(x) rmultinom(1, x[1], x[-1])))
W2[, 1:20]
