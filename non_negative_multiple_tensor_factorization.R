#####非負値テンソル因子分解#####
library(MASS)
library(matrixStats)
library(FAdist)
library(NMF)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

##データの設定
hh <- 3000   #ユーザー数
item <- 500   #アイテム数
time <- 12   #観測期間数
k <- 10   #基底数

##非負値テンソル因子分解の仮定に従いデータを生成
#ガンマ分布よりパラメータを設定
alpha01 <- 0.2; beta01 <- 0.9
alpha02 <- 0.25; beta02 <- 0.85
alpha03 <- 0.15; beta03 <- 1.5

W0 <- matrix(rgamma(hh*k, alpha01, beta01), nrow=hh, ncol=k)
H0 <- matrix(rgamma(item*k, alpha02, beta02), nrow=k, ncol=item)
C0 <- matrix(rgamma(time*k, alpha03, beta03), nrow=k, ncol=time)

#ポアソン分布の期待値を計算
WHC <- array(0, dim=c(hh, item, time))
for(j in 1:k){
  WHC <- WHC + W0[, j] %o% H0[j, ] %o% C0[j, ]
}

#ポアソン分布から購買量の生成と対数尤度の計算
Data <- array(0, dim=c(hh, item, time))
LLbest <- 0
for(rp in 1:time){
  for(j in 1:item){
    Data[, j, rp] <- rpois(hh, WHC0[, j, rp])
  }
  LLbest <- LLbest + sum(dpois(Data[, , rp], WHC0[, , rp], log=TRUE))
}
sum(Data)


####マルコフ連鎖モンテカルロ法で非負値テンソル分解を推定####
##アルゴリズムの設定
R <- 5000
keep <- 2
disp <- 10
iter <- 0

##事前分布の設定
alpha1 <- 0.01; beta1 <- 0.01
alpha2 <- 0.01; beta2 <- 0.01
alpha3 <- 0.01; beta3 <- 0.01

##初期値の設定
W <- matrix(rgamma(hh*k, 0.1, 1.0), nrow=hh, ncol=k)
H <- matrix(rgamma(item*k, 0.1, 1.0), nrow=k, ncol=item)
C <- matrix(rgamma(time*k, 0.1, 1.0), nrow=k, ncol=time)

##サンプリング結果の保存用配列
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, item, R/keep))
C_array <- array(0, dim=c(k, time, R/keep))
lambda <- array(0, dim=c(hh, item, k))


####ギブスサンプリングでパラメータをサンプリング####

