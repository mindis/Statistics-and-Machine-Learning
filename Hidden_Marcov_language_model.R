#####Hidden Marcov language model#####
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
library(ggplot2)

#set.seed(92483)

####データの発生####
##データの設定
k1 <- 7   #機能語の品詞混合数
k2 <- 5   #一般語の文法混合数
d <- 2000   #文書数
v1 <- 50   #機能語数
v2 <- 1000   #一般語数
s <- rpois(d, 15)   #文章数
s[s < 5] <- ceiling(runif(sum(s < 5), 5, 10))
a <- sum(s)   #総文章数
w <- rpois(a, 13)   #文章あたりの単語数
w[w < 5] <- ceiling(runif(sum(w < 5), 5, 10))
f <- sum(w)   #総単語数

#文書IDの設定
u_id <- rep(1:d, s)
t_id <- c()
for(i in 1:d){t_id <- c(t_id, 1:s[i])}
words <- as.numeric(tapply(w, u_id, sum))

#文章区切りのベクトルを作成
ID_d <- rep(1:d, words)
td_d <- c()
for(i in 1:d){
  td_d <- c(td_d, rep(1:s[i], w[u_id==i]))
}
nd_d <- rep(1:a, w)
x_vec <- rep(0, f)
x_vec[c(1, cumsum(w[-a])+1)] <- 1

#インデックスを設定
s_list <- list()
vec_list <- list()
for(i in 1:a){
  if(i%%1000==0){
    print(i)
  }
  s_list[[i]] <- which(nd_d==i)
  vec_list[[i]] <- rep(1, length(s_list[[i]]))
}


##パラメータの設定
#マルコフ推移行列のパラメータ
alpha01 <- rep(0.1, k2)
alpha02 <- rep(0.2, k1)
alpha03 <- rep(0.2, k2)
theta1 <- thetat1 <- as.numeric(extraDistr::rdirichlet(1, alpha01))
theta2 <- thetat2 <- extraDistr::rdirichlet(k2, alpha02)
theta3 <- thetat3 <- extraDistr::rdirichlet(k1, alpha03)

#単語出現率のパラメータ
alpha11 <- rep(0.1, v1)
alpha12 <- rep(0.1, v2)
phi1 <- phit1 <- extraDistr::rdirichlet(k1, alpha11)
phi2 <- phit2 <- extraDistr::rdirichlet(k2, alpha12)


##潜在推移と観測データの生成
Z1 <- 

Z1 <- rmnom(1, 1, theta1)









