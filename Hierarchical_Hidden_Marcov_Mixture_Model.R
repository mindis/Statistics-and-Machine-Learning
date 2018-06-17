#####Hierarchical Hidden Marcov Mixture Model#####
options(warn=2)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(data.table)
library(ggplot2)

#set.seed(5723)

####データの発生####
##データの設定
k1 <- 4   #文脈数
k2 <- 7   #状態数
k3 <- 10   #トピック数
d <- 3000   #文書数
s <- rep(300, k1)   #文脈ごとの語彙数
v <- sum(s)   #総語彙数
v1 <- c(1, (cumsum(s)+1)[-k1])
v2 <- cumsum(s)
w <- rpois(d, rgamma(d, 70, 0.5))   #文書ごとの単語数
f <- sum(w)   #総語彙数

##IDの設定
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){t_id <- c(t_id, 1:w[i])}

##パラメータの事前分布の設定
#マルコフ推移確率のディレクリ分布の事前分布を設定
alpha11 <- rep(5.0, k1)
alpha12 <- matrix(1.0, nrow=k1, ncol=k1)
diag(alpha12) <- 2.5
alpha21 <- rep(1.5, k2)
alpha22 <- matrix(0.2, nrow=k2, ncol=k2)
diag(alpha22) <- 1.0

#切り替え確率のベータ分布の事前分布を設定
s0 <- 7.5
v0 <- 50.0

#トピック分布のディリクレ分布の事前分布を設定
alpha3 <- array(0.15, dim=c(k2, k3, k1))

#トピックごとの単語分布のディリクレ分布の事前分布を設定
gamma <- matrix(0.0025, nrow=k1, ncol=v)
for(j in 1:k1){
  gamma[j, v1[j]:v2[j]] <- 0.1
}

#トピックごとの単語分布の事前分布を設定


##事前分布からパラメータを生成
#ディリクレ分布からマルコフ推移確率を生成
theta11 <- thetat11 <- <- as.numeric(extraDistr::rdirichlet(1, alpha11))
theta12 <- thetat12 <- extraDistr::rdirichlet(k1, alpha12)
theta21 <- thetat21 <- extraDistr::rdirichlet(k1, alpha21)
theta22 <- array(0, dim=c(k2, k2, k1))
for(j in 1:k1){
  theta22[, , j] <- extraDistr::rdirichlet(k2, alpha22)
}
thetat22 <- theta22

#ベータ分布から切り替え変数を生成
beta <- betat <- rbeta(d, s0, v0)


#ディリクレ分布からトピック分布を生成
theta3 <- array(0, dim=c(k2, k3, k1))
for(j in 1:k1){
  theta3[, , j] <- extraDistr::rdirichlet(k2, alpha3[, , j])
}
thetat <- theta

#ディリクレ分布から単語分布を生成
phi <- array(0, dim=c(k3, v, k1))
for(j in 1:k1){
  phi[, , j] <- extraDistr::rdirichlet(k3, gamma[j, ])
}
phit <- phi
phit[, 1:10, ]














