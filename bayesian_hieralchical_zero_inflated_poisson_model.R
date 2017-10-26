#####階層ベイズゼロ過剰ポアソン回帰モデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(extraDistr)
library(matrixStats)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####データの発生####
##データの設定
hh <- 2000   #サンプル数
pt <- ceiling(rgamma(hh, 4.8, 0.8))   #購買機会
hhpt <- sum(pt)

##IDの設定
id <- rep(1:hh, pt)
time <- c()
for(i in 1:hh){time <- c(time, 1:pt[i])}
ID <- data.frame(no=1:length(id), id, time)

####説明変数の発生####
##階層モデルの説明変数の発生
cont1 <- 3; bin1 <- 3; multi1 <- 4
X.cont <- matrix(rnorm(hh*cont1), nrow=hh, ncol=cont1)
X.bin <- matrix(0, nrow=hh, ncol=bin1)
X.multi <- matrix(0, nrow=hh, ncol=multi1)

#二値説明変数を設定
for(i in 1:bin1){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(hh, 1, p)
}

#多値説明変数を設定
p <- runif(multi1)
X.multi <- t(rmultinom(hh, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))] #冗長な変数は削除

#データを結合
ZX <- cbind(1, X.cont, X.bin, X.multi)


##ゼロ過剰ポアソン回帰モデルの説明変数
#発売ジャンル
quantity <- rpois(hhpt, 15)
status <- cbind(quantity, matrix(c(0.2, 0.3, 0.2, 0.15, 0.15), nrow=hhpt, ncol=5, byrow=T))
genre <- (t(apply(status, 1, function(x) rmultinom(1, x[1], x[-1])))/quantity)[, -5]

#プロモーション有無
promo <- rbinom(hhpt, 1, 0.4)

#イベント有無
event <- rbinom(hhpt, 1, 0.25)

#訪問回数
visit0 <- rep(0, hhpt)
for(i in 1:hh){
  par1 <- runif(1, 1.2, 4.7)
  par2 <- runif(1, 0.8, 2.5)
  visit0[ID$id==i] <- log(round(rgamma(sum(ID$id==i), par1, par2))) 
}
visit <- visit0 + 1
visit[is.infinite(visit)] <- 0
summary(visit)

##データの結合
X <- data.frame(bp=1, genre=genre, promo, event, visit)
XM <- as.matrix(X)


####応答変数の発生####
##購買潜在変数zの発生
z0 <- rep(0, hhpt)
for(i in 1:hh){
  z0[ID$id==i] <- rep(rbinom(1, 1, 0.7), pt[i])
}


