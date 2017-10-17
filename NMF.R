#####選択の曖昧性のためのマルチクラス分類モデル####
library(MASS)
library(matrixStats)
library(flexmix)
library(glmnet)
library(mlogit)
library(nnet)
library(FAdist)
library(NMF)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

####データの発生####
hh <- 3000
seg <- 2

##説明変数の発生
#もとの説明変数行列を発生
topic <- 10   #トピック数
k1 <- 300   #説明変数数
freq <- rpois(hh, 200)   #ポアソン分布から頻度を発生

#ディレクリ分布から出現確率を発生
#パラメータの設定
alpha0 <- runif(topic, 0.1, 1.0)   #文書のディレクリ事前分布のパラメータ
theta0 <- rdirichlet(hh, alpha0)   #文書のトピック分布をディレクリ乱数から発生

alpha1 <- matrix(0, nrow=topic, ncol=k1)
phi0 <- matrix(0, nrow=topic, ncol=k1)
for(i in 1:topic){
  alpha1[i, ] <- rgamma(k1, 0.4, 0.1)   #単語のディレクリ事前分布のパラメータ
  phi0[i, ] <- rdirichlet(1, alpha1[i, ])   #単語のトピック分布をディレクリ乱数から発生
}

#多項分布の乱数からデータを発生
X0 <- matrix(0, hh, k1)
Topic <- list()

for(i in 1:hh){
  z <- t(rmultinom(freq[i], 1, theta0[i, ]))   #文書のトピック分布を発生
  
  zn <- z %*% c(1:topic)   #0,1を数値に置き換える
  zdn <- cbind(zn, z)   #apply関数で使えるように行列にしておく
  
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi0[x[1], ])))   #文書のトピックから単語を生成
  wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
  X0[i, ] <- wdn  
  Topic[[i]] <- zdn[, 1]
  print(i)
}

##非負値行列因子分解で頻度行列を圧縮
k <- 10
X0_trance <- t(X0)
res1 <- nmf(X0, k, "brunet")
res2 <- nmf(X0_trance, k, "brunet")

H <- res2@fit@H
W <- res2@fit@W

test <- round(cbind(t(solve(t(W) %*% W) %*% t(W) %*% X0_trance), t(H)), 3)


#真のトピックの出現確率と推定されたトピック確率を比較
t_rate <- matrix(0, hh, topic) 
for(i in 1:hh){
  rate0 <- table(Topic[[i]])/freq[i]
  rate <- rep(0, topic)
  index <- subset(1:topic, 1:topic %in% names(rate0))
  rate[index] <- rate0
  t_rate[i, ] <- rate
}

#最適なトピック数は10
opt.topic <- 10
Topic_rate <- round(cbind(t_rate, t(res2[[opt.topic-1]]@fit@H)/matrix(rowSums(t(res2[[opt.topic-1]]@fit@H)), nrow=hh, ncol=topic)), 3)

#トピックの出現確率を説明変数とする
X <- (t(res2[[opt.topic-1]]@fit@H)/matrix(rowSums(t(res2[[opt.topic-1]]@fit@H)), nrow=hh, ncol=opt.topic))[, -opt.topic]
XM <- cbind(1, X)