#####Hierarchical Probablistic Automaton model#####
options(warn=0)
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
k1 <- 6   #上位階層
k2 <- 6   #中間階層
k3 <- 10   #下位階層
d <- 3000   #文書数
v1 <- 700   #一般語
v2 <- 200   #独立語
v3 <- 200   #機能語
v4 <- 25   #区切り記号
v <- v1 + v2 + v3 + v4   #総語彙数
w <- rpois(d, rgamma(d, 65, 0.4))   #語彙数
f <- sum(w)   #総語彙数

#IDの設定
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}

##パラメータの事前分布の設定
#初期分布の事前分布
pi01 <- rep(5.0, k1-1)

#マルコフ推移行列の事前分布
alpha01 <- c(rep(50.0, k1-1), 30)   #最上位階層の事前分布
















