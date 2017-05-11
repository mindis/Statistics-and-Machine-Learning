#####変量効果ロジスティック回帰モデル#####
library(MASS)
library(glmm)
library(lme4)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
hh <- 1000   #消費者数
pt <- rpois(hh, 7)   #ひとりひとりの接触数
pt <- ifelse(pt==0, 1, pt)   #接触数が0なら1に置き換え
hhpt <- sum(pt)   #全サンプル数
col.fix <- 8   #固定効果の変数数
col.random <- 3   #ランダム効果数

##IDの記録と説明変数の発生
#idの記録
id <- rep(1:hh, rep(pt, 1))

#tの記録
t <- c()
for(i in 1:hh){
  ti <- 1:pt[i]
  t <- c(t, ti)
}
#データの結合
ID <- data.frame(no.=1:length(id), id=id, t=t)

#説明変数の発生
X1 <- matrix(rnorm(hhpt*(col.fix-3), 0, 1), nrow=hhpt, ncol=(col.fix-3))
X2 <- matrix(0, hhpt, col.fix-ncol(X1))
for(i in 1:(col.fix-ncol(X1))){
  bin <- rbinom(hhpt, 1, runif(1, 0.2, 0.7))
  X2[, i] <- bin
}
X <- data.frame(cont=X1, bin=X2)   #固定効果のデータ
Z <- data.frame(1, X$cont.1, X$bin.1)   #ランダム効果のデザイン行列

##回帰係数の発生
#固定効果の回帰係数の発生
beta1.fix <- c(runif(ncol(X1), 0, 1.2), runif(ncol(X2), -1.0, 1.0))
beta0.fix <- 0.5

#変量効果の発生
var.v <- c(0.4^2, 0.3^2, 0.3^2)
beta.M <- matrix(0, hhpt, col.random)
for(i in 1:hh){
  random <- rmvnorm(1, rep(0, col.random), diag(var.v))
  beta.M[ID$id==i, ] <- matrix(random, nrow=length(ID$id[ID$id==i]), ncol=col.random, byrow=T)
}
beta.random <- beta.M[ID$t==1, ]   #ユニークな変量効果行列

##反応変数の発生
#ロジットの計算
logit.fix <- beta0.fix + as.matrix(X) %*% beta1.fix 
logit.random <- rowSums(Z * beta.random)
logit <- logit.fix + logit.random

#確率の計算
P <- exp(logit)/(1+exp(logit))
hist(P, col="#0000ff40", border = "#0000ff",  breaks=20, 
     main="変量効果モデルの確率の分布", xlab="確率", ylab="頻度")

#ベルヌーイ乱数で応答変数を発生
y <- c()
for(i in 1:hhpt){
  y.bin <- rbinom(1, 1, P[i])
  y <- c(y, y.bin)
}
table(y)
round(YXZ <- data.frame(ID, y, P, random=beta.M), 3)

##ランダム効果の結果を可視化
#挙動を見る変数の範囲
val <- seq(-3.0, 3.0, length=200)

#変量効果と固定効果をプロット
mu.f <- beta0.fix + val * beta1.fix[1] 
p.f <- exp(mu.f)/(1+exp(mu.f))
plot(val, p.f, type="l", lwd=2, col=2, main="変量効果の可視化", ylab="確率", xlab="value", ylim=c(0, 1))
for(i in 1:(hh/2)){
  mu.r <- mu.f + beta.random[i, 1] + val * beta.random[i, 2]
  p.r <- exp(mu.r)/(1+exp(mu.r))
  lines(val, p.r, type="l")
}
lines(val, p.f, type="l", lwd=2, col=2)


####マルコフ連鎖モンテカルロ法で変量効果ロジスティック回帰モデルを推定####

