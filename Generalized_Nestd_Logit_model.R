#####generalized nested logit model#####
library(MASS)
library(mlogit)
library(nnet)
library(flexmix)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(8437)

####データの設定####
g <- 3   #ネスト数
hh <- 2000   #サンプル数
pt <- rpois(hh, 5); pt <- ifelse(pt==0, 1, pt)   #購買機会(購買機会数が0なら1に置き換え)
hhpt <- sum(pt)
member <- 9   #選択可能メンバー数
k <- 5   #説明変数の数

##IDの設定
id <- rep(1:hh, pt)
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}

#IDとセグメントを結合
ID <- data.frame(no=1:hhpt, id=id, t=t)   #データの結合


####説明変数の発生####
#衣装の設定
c.num <- member
CLOTH <- list()
for(i in 1:member){
  CLOTH[[i]] <- t(rmultinom(hhpt, 1, runif(c.num)))
  CLOTH[[i]] <- CLOTH[[i]][, -c.num]
}

#レベルの対数
lv.weib <- round(rweibull(hh*2, 1.8, 280), 0)
index.lv <- sample(subset(1:length(lv.weib), lv.weib > 80), hh)
lv <- log(lv.weib[index.lv])

#パネルに変更
LV <- c()
for(i in 1:hh){
  LV <- c(LV, rep(lv[i], pt[i]))
}

#スコアの対数
score.norm <- exp(rnorm(hhpt*2, 12.5, 0.5))
index.score <- sample(subset(1:length(score.norm), score.norm > 150000), hhpt)
score <- log(score.norm[index.score])
SCORE <- score

#どのメンバーの勧誘回だったか
prob <- 1/(member)
scout <- t(rmultinom(hhpt, 2, rep(prob, member)))

#メンバーで勧誘が重複しなくなるまで乱数を発生させ続ける
for(i in 1:10000){
  if(max(scout)==1) break
  index.scout <- subset(1:nrow(scout), apply(scout, 1, max) > 1)
  scout[index.scout, ] <- t(rmultinom(length(index.scout), 2, rep(prob, member)))
  print(i)
}
SCOUT <- scout


##説明変数をベクトル形式に変換
#切片の設定
p <- c(1, rep(0, member))
Pop <- matrix(p, nrow=hhpt*length(p), ncol=member, byrow=T)
Pop <- subset(Pop, rowSums(Pop) > 0)[, -member]

#多項型説明変数をベクトル形式に変更
LV.v <- matrix(0, nrow=hhpt*member, ncol=member-1)
SCORE.v <- matrix(0, nrow=hhpt*member, ncol=member-1)

for(i in 1:hhpt){
  index.v <- ((i-1)*member+1):((i-1)*member+member)
  LV.v[index.v, ] <- diag(LV[i], member)[, -member]
  SCORE.v[index.v, ] <- diag(SCORE[i], member)[, -member]
}

#条件付き説明変数をベクトル形式に変更
CLOTH.v <- matrix(0, nrow=hhpt*member, ncol=c.num-1)
for(i in 1:hhpt){
  print(i)
  for(j in 1:member){
    index.v <- (i-1)*member+j
    CLOTH.v[index.v, ] <- CLOTH[[j]][i, ]
  }
}

SCOUT.v <- as.numeric(t(SCOUT))

#データを結合
X <- data.frame(pop=Pop, lv=LV.v, score=SCORE.v, cloth=CLOTH.v, scout=SCOUT.v)
XM <- as.matrix(X)


####GNLモデルに基づき応答変数を発生####
##ネストを設定する
mus <- c("maki", "rin", "hanayo", "eri", "nozomi", "nico", "kotori", "umi", "honoka")

#学年ごと
first <- c(1, 1, 1, 0, 0, 0, 0, 0, 0)
second <- c(0, 0, 0, 1, 1, 1, 0, 0, 0)
third <- c(0, 0, 0, 0, 0, 0, 1, 1, 1)

#ミニユニット
Prim <- c(0, 0, 1, 0, 0, 0, 1, 0, 1)
BiBi <- c(1, 0, 0, 1, 0, 1, 0, 0, 0)
LW <- c(0, 1, 0, 0, 1, 0, 0, 1, 0)

#スクフェス属性
smile <- c(0, 1, 0, 0, 0, 1, 0, 0, 1)
pure <- c(0, 0, 1, 0, 1, 0, 1, 0, 0)
cool <- c(1, 0, 0, 1, 0, 0, 0, 1, 0)

#ネストを結合
nest <- rbind(first, second, third, Prim, BiBi, LW, smile, pure, cool)


##パラメータを設定
#回帰パラメータの設定
b0 <- c(1.9, -0.6, -0.9, 0.9, -1.2, 1.5, 1.7, 0.4)
b1 <- runif((member-1)*2, -0.7, 0.8)
b2 <- runif(NCOL(CLOTH.v)+NCOL(SCOUT.v), -1.0, 1.8)
b <- c(b0, b1, b2)
beta.t <- b


#ログサム変数のパラメータ
grade <- runif(g, 0.3, 0.6)
unit <- runif(g, 0.5, 0.8)
type <- runif(g, 0.3, 0.5)
logsum.par <- c(grade, unit, type)

#アロケーションパラメータの設定
gamma.k1 <- runif(member, 1.5, 4.5)
gamma.k2 <- runif(member, 0.5, 3.5)
gamma.k3 <- runif(member, 2.5, 4.5)

#アロケーションパラメータを正規化
gamma.grade <- gamma.k1/sum(gamma.k1)
gamma.unit <- gamma.k2/sum(gamma.k2)
gamma.type <- gamma.k3/sum(gamma.k3) 
gamma.par <- rbind(gamma.grade, gamma.unit, gamma.type)
gamma <- matrix(0, nrow=nrow(nest), ncol=3)

for(i in 1:g){
  for(j in 1:3){
    r <- (i-1)*3+j
    gamma[r, ] <- (gamma.par[i, ]*nest[r, ])[nest[r, ]==1]
  }
}


##GNLモデルに基づき確率を計算
#ロジットを計算
logit <- matrix(XM %*% b, nrow=hhpt, ncol=member, byrow=T)

##ネストの所属確率を計算
#ネストごとにログサム変数を計算
logsum <- matrix(0, nrow=hhpt, ncol=nrow(nest)) 
d2_2 <- matrix(0, nrow=hhpt, ncol=nrow(nest))
d2_1 <- array(0, dim=c(hhpt, member, nrow(nest)))

for(i in 1:nrow(nest)){
  #ネスト、ログサム、アロケーションパラメータをメンバーで行列に変更
  nest.k <- matrix(nest[i, ], nrow=hhpt, ncol=member, byrow=T)
  gamma.k <- matrix(gamma[i, ], nrow=hhpt, ncol=g, byrow=T)
  
  #ログサム変数を計算
  logsum[, i] <- logsum.par[i] * log(rowSums((gamma.k * exp((logit*nest.k)[, nest[i, ]==1]))^(1/logsum.par[i])))
  d2_2[, i] <- rowSums((gamma.k * exp((logit*nest.k)[, nest[i, ]==1]))^(1/logsum.par[i]))
  d2_1[, nest[i, ]==1, i] <- (gamma.k * exp((logit*nest.k)[, nest[i, ]==1]))^(1/logsum.par[i])
}


#ネストjの選択確率のパラメータを計算
U1_1 <- exp(logsum)
U1_2 <- matrix(rowSums(exp(logsum)), nrow=hhpt, ncol=nrow(nest))
U2_1 <- matrix(0, nrow=hhpt, ncol=member)
U2_2 <- matrix(0, nrow=hhpt, ncol=member)

#メンバーごとにネストに所属する確率の合計を計算
Pr1 <- matrix(0, nrow=hhpt, ncol=member)
for(i in 1:member){
  i <- 1
  Pr1[, i] <- rowSums((U1_1 / U1_2)[, nest[, i]==1])
  U2_2[, i] <- rowSums(d2_2[, nest[, i]==1])
  U2_1[, i] 
  d2_1
}

##ネストを条件つけた選択確率を計算




