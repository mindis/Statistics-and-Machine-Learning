#####LASSO#####
library(MASS)
library(kernlab)
library(quadprog)
library(glmnet)
library(lars)
library(reshape2)
library(plyr)
####データの発生####
set.seed(2234)
n <- 1000   #サンプル数
p <- 300   #説明変数の数
b <- c(3, 1.2, 0.8, 2.0, -1.4, -1.5, 2.2, -0.7, 1.8, 1.0, rep(0, p-9))   #真の回帰係数
X <- cbind(1, matrix(rnorm(n*p), nrow=n, ncol=p, byrow=T))   #説明変数の発生
Z <- X %*% b   #真の平均構造
Y <- Z + rnorm(n, 0, 2)   #応答変数の発生

####L1正則化lassoをL推定####
##最小化する関数を設定
Q <- t(X) %*% X
c <- -2*(t(Y) %*% X)

##制約条件
A <- t(rep(1, p+1))
b <- rep(-0.5, p+1)

##凸二次計画問題を解く
sv <- solve.QP(H, c, A, b)

http://qiita.com/hogefugabar/items/71916560f7efc6eededf
