#####入れ子型マルチレベルモデル#####
library(MASS)
library(nlme)
library(glmm)
library(reshape2)
library(plyr)
library(lattice)
library(ggplot2)

####データの発生####
#set.seed(94327)
n <- 1000   #評価者数
g1 <- 50   #対象アニメ数
g2 <- round(runif(g1, 2.3, 5.8))   #登場キャラ数
g2s <- sum(g2)

####対象アニメを視聴しているかどうかを発生####
##説明変数の発生
cont <- 3; bin <- 4; multi <- 4
X.cont <- matrix(rnorm(n*cont), nrow=n, ncol=cont)
X.bin <- matrix(0, nrow=n, ncol=bin)
X.multi <- matrix(0, nrow=n, ncol=multi)

#二値説明変数を設定
for(i in 1:bin){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(n, 1, p)
}

#多値説明変数を設定
p <- runif(multi)
X.multi <- t(rmultinom(n, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))] #冗長な変数は削除

#データを結合
X <- cbind(1, X.cont, X.bin, X.multi)

##アニメの割当を発生
#パラメータの設定
alpha0 <- runif(g1, -1.2, 1.0)
alpha1 <- matrix(runif(g1*cont, 0, 0.6), nrow=cont, ncol=g1)
alpha2 <- matrix(runif(g1*(bin+multi-1), -0.7, 0.6), nrow=bin+multi-1, ncol=g1)
alpha <- rbind(alpha0, alpha1, alpha2)

#ロジットと確率の計算
logit <- X %*% alpha
Pr <- exp(logit)/(1+exp(logit))

#二項分布から割当を発生
R <- apply(Pr, 2, function(x) rbinom(n, 1, x))
colMeans(R); mean(R)


####デザイン行列を定義####
#デザイン行列の格納用配列
Z1 <- matrix(0, nrow=n*g2s, ncol=n)
Z2 <- matrix(0, nrow=n*g2s, ncol=g1)
Z3 <- matrix(0, nrow=n*g2s, ncol=g2s)

#インデックスを作成
index_g21 <- c(1, cumsum(g2))
index_g21[2:length(index_g21)] <- index_g21[2:length(index_g21)] + 1
index_g22 <- cumsum(g2)

for(i in 1:n){
  print(i)
  #個人別の格納用配列
  z2 <- matrix(0, nrow=g2s, ncol=g1)
  z3 <- matrix(0, nrow=g2s, ncol=g2s)
  
  r <- ((i-1)*g2s+1):((i-1)*g2s+g2s)
  Z1[r, i] <- 1
  
  for(j in 1:g1){
    if(R[i, j]==1){
      z2[index_g21[j]:index_g22[j], j] <- 1
      z3[index_g21[j]:index_g22[j], index_g21[j]:index_g22[j]] <- diag(g2[j])
    }
  }
  Z2[r, ] <- z2
  Z3[r, ] <- z3
}

Z <- cbind(Z1, Z2, Z3)   #データを結合

#評価していないアニメは欠損させる
index_zeros <- subset(1:nrow(Z2), rowSums(Z2)==0)
Z <- Z[-index_zeros, ]
Z




