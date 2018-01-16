#####階層pitman-yor過程言語モデル#####
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
library(ggplot2)

#set.seed(5723)
#データの設定
d <- 1000   #文書数
v <- 1000   #語彙数
w <- rpois(d, rgamma(d, 70, 0.40))   #1文書あたりの単語数
f <- sum(w)

#bi-gramモデルのパラメータを設定
par1 <- c(5.0, 0.5)
par2 <- c(0.5, 0.8)
alpha0 <- matrix(0, nrow=v, ncol=v) 
for(i in 1:v){
  index <- rbinom(1, 1, 0.3) + 1
  alpha0[i, ] <- rep(rgamma(1, par1[index], par2[index]), v)
}
theta <- extraDistr::rdirichlet(v, alpha0) 
summary(as.numeric(theta))

