#####Multilayer Topic Model#####
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

#set.seed(93471)

####データの発生####
#データの設定
k0 <- 4   #上位階層のトピック数
k1 <- 4   #下位階層のトピック数
k2 <- 4
k3 <- 5
k <- c(k1, k2, k3)
k_sums <- sum(k)
d <- 2000   #文書数
v <- 700   #語彙数
s <- length(k) 
w <- rpois(d, rgamma(d, 100, 0.6))
f <- sum(w)

#文書IDの設定
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){t_id <- c(t_id, 1:w[i])}

#語彙トピックのインデックス
for(j in 1:100){
  v_index <- c(1, round(sort(runif(k-1, 100, v-100))), v)
  x <- v_index[2:length(v_index)] - v_index[1:(length(v_index)-1)]
  if(min(x) > 100) break
}

##パラメータの設定
#ディリクレ分布のパラメータ
alpha01 <- rep(0.25, k0)
alpha11 <- rep(0.2, k1)
alpha12 <- rep(0.2, k2)
alpha13 <- rep(0.2, k3)
beta01 <- rep(0.25, s)
beta11 <- beta12 <- beta13 <- rep(0.001, v)
beta11[v_index[1]:v_index[2]] <- 0.1
beta12[(v_index[2]+1):v_index[4]] <- 0.1
beta13[(v_index[4]+1):v_index[5]] <- 0.1

#パラメータを生成
theta0 <- thetat0 <- extraDistr::rdirichlet(d, alpha01)
theta1 <- thetat1 <- extraDistr::rdirichlet(d, alpha11)
theta2 <- thetat2 <- extraDistr::rdirichlet(d, alpha12)
theta3 <- thetat3 <- extraDistr::rdirichlet(d, alpha13)
phi0 <- phit0 <- extraDistr::rdirichlet(k0, beta01)
phi1 <- phit1 <- extraDistr::rdirichlet(k1, beta11)
phi2 <- phit2 <- extraDistr::rdirichlet(k2, beta12)
phi3 <- phit3 <- extraDistr::rdirichlet(k3, beta13)


##階層トピックモデルからデータを生成
Z1_list <- list()
Z2

for(i in 1:d){
  #上位階層のトピックを生成
   <- rmnom(w[i], 1, theta0[i, ])
}



