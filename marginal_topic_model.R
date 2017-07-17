#####周辺化トピックモデル#####
library(MASS)
library(lda)
library(RMeCab)
library(gtools)
library(reshape2)
library(plyr)
library(ggplot2)


####データの発生####
#set.seed(423943)
#データの設定
k <- 5   #トピック数
d <- 200   #文書数
v <- 100   #語彙数
w <- 200   #1文書あたりの単語数 

#パラメータの設定
alpha0 <- runif(k, 0.1, 1.25)   #文書のディレクリ事前分布のパラメータ
alpha1 <- rep(0.5, v)   #単語のディレクリ事前分布のパラメータ

#ディレクリ乱数の発生
theta0 <- rdirichlet(d, alpha0)   #文書のトピック分布をディレクリ乱数から発生
phi0 <- rdirichlet(k, alpha1)   #単語のトピック分布をディレクリ乱数から発生

#多項分布の乱数からデータを発生
WX <- matrix(0, d, v)
ZS <- list()

for(i in 1:d){
  z <- t(rmultinom(w, 1, theta0[i, ]))   #文書のトピック分布を発生
  
  zn <- z %*% c(1:k)   #0,1を数値に置き換える
  zdn <- cbind(zn, z)   #apply関数で使えるように行列にしておく
  
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi0[x[1], ])))   #文書のトピックから単語を生成
  wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
  WX[i, ] <- wdn  
  ZS[[i]] <- zdn[, 1]
  print(i)
}
barplot(colSums(z), names.arg=c("seg1", "seg2", "seg3", "seg4", "seg5"))
round(colSums(WX)/sum(WX), 3)   #単語の出現頻度


####周辺化ギブスサンプリングでトピックモデルを推定####
