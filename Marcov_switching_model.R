#####マルコフ切り替えモデル#####
library(MASS)
library(MSwM) 
library(reshape2)
library(gtools)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
n <- 500   #サンプル数
k1 <- 4   #切り替え数
k2 <- 6   #観測確率のパラメータ数


##初期確率の定義
Pf <- c(0.4, 0.1, 0.3, 0.2)

##推移行列の定義
pr1 <- c(0.2, 0.1, 0.4, 0.3)
pr2 <- c(0.1, 0.5, 0.2, 0.2)
pr3 <- c(0.3, 0.1, 0.3, 0.3)
pr4 <- c(0.1, 0.2, 0.3, 0.4)
Pr <- rbind(pr1, pr2, pr3, pr4)

##観測確率の定義
P <- matrix(0, nrow=k1, ncol=k2)
for(i in 1:k1){
  alpha <- runif(1, 0.4, 1)
  P[i, ] <- rdirichlet(1, rep(alpha, k2))
}


##応答変数の発生
Z <- matrix(0, nrow=n, ncol=k1)
Y <- matrix(0, nrow=n, ncol=k2)

#潜在変数の初期値
Z[1, ] <- rmultinom(1, 1, Pf)
Y[1, ] <- rmultinom(1, 1, P[which.max(Z[1, ]), ])

#2回目以降の応答変数を逐次的に発生させる
for(i in 2:n){
  Z[i, ] <- rmultinom(1, 1, Pr[which.max(Z[i-1, ]), ])
  Y[i, ] <- rmultinom(1, 1, P[which.max(Z[i, ]), ])
}


####EMアルゴリズムでマルコフ切り替えモデルを推定####

