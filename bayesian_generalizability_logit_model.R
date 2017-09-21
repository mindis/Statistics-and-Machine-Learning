#####変量効果ロジットモデルによる一般化可能性理論#####
library(MASS)
library(nlme)
library(glmm)
library(survival)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(57389)

####データの発生####
n1 <- 500   #ユーザー数
n2 <- 150   #対象ゲーム数
N <- n1*n2   #総サンプル数
genre <- 6   #ゲームジャンル数

##IDの設定
u.id <- rep(1:n1, rep(n2, n1))
g.id <- rep(1:n2, n1)
ID <- data.frame(no=1:N, u.id, g.id)

####説明変数の発生####
##ゲームの特性変数行列を作成
Disc <- rbinom(n2, 1, 0.6)
Genre0 <- t(rmultinom(n2, 1, runif(genre, 0.5, 2.0)))
Genre <- Genre0[, -which.min(colSums(Genre0))]
X1 <- cbind(1, Disc, Genre)   #データの結合

##変量効果のデザイン行列を設定
Z1 <- matrix(0, nrow=N, ncol=n1*ncol(X1))
Z2 <- matrix(0, nrow=N, ncol=n2)
for(i in 1:n1){
  print(i)
  r <- ((i-1)*ncol(X1)+1):((i-1)*ncol(X1)+ncol(X1))
  Z1[ID$u.id==i, r] <- X1
}

for(i in 1:n2){
  print(i)
  Z2[ID$g.id==i, i] <- 1 
}

####応答変数の発生####
##パラメータの設定
#平均構造を設定
theta0 <- runif(1, -1.0, -0.6)

#個人嗜好の変量効果の設定
Cov01 <- diag(c(0.75, 0.35, runif(genre-1, 0.2, 0.4)))
alpha0 <- as.numeric(t(mvrnorm(n1, rep(0, ncol(X1)), Cov01)))

#ゲーム価値の変量効果の設定
Cov02 <- runif(1, 0.75, 0.8)
beta0 <- rnorm(n2, 0, Cov02)

##ロジットと確率の計算
logit <- theta0 + Z1 %*% alpha0 + Z2 %*% beta0
Pr0 <- exp(logit)/(1+exp(logit))

##応答変数の発生
Y <- rbinom(N, 1, Pr0)

#発生させた変数の可視化と集計
mean(Y); sum(Y); mean(Pr0); summary(Pr0)
hist(Pr0, col="grey", xlab="確率", main="購買確率の分布")


