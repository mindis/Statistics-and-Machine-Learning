#####Tree Structured Mixture Multinomial Model#####
options(warn=0)
library(stringr)
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
#set.seed(2506787)
####データの発生####
##データの設定
#木構造のトピックの設定
m <- 3   #木の最大深さ
k1 <- 3   #1階層目のトピック数
k2 <- rtpois(k1, a=1, b=Inf, 2.5)   #2階層目のトピック数
k3 <- c()
for(j in 1:k1){
  k3 <- c(k3, rtpois(k2[j], a=0, b=Inf, 2.0))   #3階層目のトピック数
}
index_k3 <- cbind(c(1, cumsum(k2)[-length(k2)]+1), cumsum(k2))
k <- 1 + sum(c(k1, k2, k3))   #総トピック数

#文書の設定
d <- 3000   #文書数
w <- rpois(d, rgamma(d, 75, 0.5))   #単語数
f <- sum(w)   #総単語数
v0 <- 250 
v1 <- rep(300, k1)
v <- sum(c(v0, v1))   #総語彙数
index_v0 <- 1:v0
index_v1 <- cbind(c(v0+1, (v0+cumsum(v1)[-k1])+1), (v0+cumsum(v1)))

#IDを設定
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}


##パラメータを生成
#


