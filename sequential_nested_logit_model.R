#####逐次型ネステッドロジットモデル#####
library(MASS)
library(mlogit)
library(nnet)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

#####データの発生#####
#set.seed(89361)
##データの設定
hh <- 1000   #プレイヤー数
pt <- 120   #観測期間
hhpt <- hh*pt   #全サンプル数

#IDの設定
ID <- matrix(0, hhpt, 3)
ID[, 1] <- 1:nrow(ID)
ID[, 2] <- rep(1:hh, rep(pt, hh))
ID[, 3] <- rep(1:pt, hh)


##パラメータの設定
#ログイン有無のパラメータの設定
login.para <- 7   #ログインモデルのパラメータ数
alpha.s0 <- runif(1, -1.6, -1.3)   #切片
alpha.s1 <- runif(1, 0.9, 1.5)   #過去60日のログイン率
alpha.s2 <- runif(1, -0.4, -0.3)   #前回からのログイン経過時間(単位=日)
alpha.s3 <- runif(1, 0.5, 1.2)   #前日のログイン有無
alpha.s4 <- runif(1, 0.13, 0.24)   #1日あたりのプレイ数平均の対数
alpha.s5 <- runif(1, 0.5, 1.0)   #イベント有無
alpha.s6 <- runif(1, 0.3, 0.6)   #プレイヤーのスキルレベル
alpha.s7 <- runif(1, 0.4, 0.9)   #ログサム変数のパラメータ
alpha.s <- c(alpha.s1, alpha.s2, alpha.s3, alpha.s4, alpha.s4, alpha.s6, alpha.s7)

#ガチャ有無のパラメータの設定
buy.para <- 14   #課金有無のパラメータ数
beta.s0 <- runif(1, -4.2, -3.6)   #切片
beta.s1 <- runif(1, 0.7, 1.2)   #過去60日のガチャ率
beta.s2 <- runif(1, -0.7, -0.5)   #過去7日のガチャ有無
beta.s3 <- runif(1, 0.15, 0.3)   #1回あたりのガチャ回数
beta.s4 <- runif(1, 0.3, 0.6)   #前回のガチャからの経過日数の対数
beta.s5 <- runif(1, 0.1, 0.22)   #1日あたりのプレイ数平均の対数
beta.s6 <- runif(9, 0.8, 1.6)   #URの確率アップキャラ
beta.s <- c(beta.s1, beta.s2, beta.s3, beta.s4, beta.s5, beta.s6)

##初期値の設定
X.login <- matrix(0, nrow=hhpt, ncol=login.para)
X.buy <- matrix(0, nrow=hhpt, ncol=buy.para)

#1日目を取り出す
index.f <- subset(1:nrow(X.login), ID[, 3]==1)

##初期データの発生
##ログインデータの初期値
#過去60日のログイン履歴
login.hist <- matrix(0, nrow=hh, ncol=60)
for(i in 1:hh){
 p <- runif(1, 0.2, 1.0) 
 login.hist[i, ] <- rbinom(60, 1, p)
}
login.ratef <- rowMeans(login.hist)
X.login[index.f, 1] <- login.ratef   #過去60日のログイン率

#前回からのログイン経過時間
login.last <- apply(login.hist, 1, function(x) tail(subset(1:length(x), x==1), 1))
X.login[index.f, 2] <- 61-login.last   #前回からのログイン経過時間   

#前日ログインの有無
login.pre <- ifelse(61-login.last==1, 1, 0)
X.login[index.f, 3] <- login.pre

#1日当たりのプレイ数平均
play.M <- matrix(0, nrow=hh, ncol=60) 

for(i in 1:hh){
  #平均プレイ回数の記録
  pois <- runif(1, 5, 18)
  play <- rpois(sum(login.hist[i, ]), pois)
  X.login[index.f[i], 4] <- log(sum(play)/sum(login.hist[i, ]))
  
  #プレイ履歴の記録
  index.play <- subset(1:length(login.hist[i, ]), login.hist[i, ]==1)
  play.M[i, index.play] <- play
}  
hist(X.login[index.f, 4], breaks=20, col="grey", xlab="1日あたりプレイ数平均",main="1日あたりプレイ数平均の分布")

#イベントの有無
X.login[, 5] <- rep(c(rep(0, 5), rep(1, 10)), hh)


#プレイヤースキル
for(i in 1:hh){
  X.login[ID[, 2]==i, 6] <- rnorm(1, 0, 1)   
}

##ガチャデータの初期値
#過去60日のガチャ率
buy.hist <- matrix(0, nrow=hh, ncol=60)

for(i in 1:hh){
  logit.bpast <- beta.s0 + runif(1, 0.5, 3.0)
  Pr.gpast <- exp(logit.bpast)/(1+exp(logit.bpast))
  buy.hist[i, ] <- rbinom(60, 1, Pr.gpast)
}
X.buy[index.f, 1]  <- rowSums(buy.hist)/60

#過去7日のガチャ有無
X.buy[index.f, 2] <- ifelse(rowSums(buy.hist[, 54:60]) > 0, 1, 0)   


#1回あたりのガチャ回数
cnt.zeros <- rpois(hh, 2.5)
g.cnt <- ifelse(cnt.zeros==0, 1, cnt.zeros)
g.cnt[index.pp==0] <- 0
X.buy[, 3] <- g.cnt

#1日当たりのプレイ数平均
buy.M <- matrix(0, nrow=hh, ncol=60) 

for(i in 1:hh){
  #平均ガチャ回数の記録
  
  play <- rpois(sum(login.hist[i, ]), pois)
  X.login[index.f[i], 4] <- log(sum(play)/sum(login.hist[i, ]))
  
  #プレイ履歴の記録
  index.play <- subset(1:length(login.hist[i, ]), login.hist[i, ]==1)
  play.M[i, index.play] <- play
}  
hist(X.login[index.f, 4], breaks=20, col="grey", xlab="1日あたりプレイ数平均",main="1日あたりプレイ数平均の分布")


#前回のガチャからの経過日数
buy.progress <- rep(0, hh)
index.buy <- subset(1:nrow(buy.hist), rowSums(buy.hist) > 0)
buy.last <- apply(buy.hist[index.buy, ], 1, function(x) tail(subset(1:length(x), x==1), 1))
buy.progress[index.buy] <- 61-buy.last
X.buy[index.f, 4] <- ifelse(buy.progress > 0, log(buy.progress), 0)   #前回からのログイン経過時間  


day.r <- c()
for(i in 1:hh){
  day.r <- c(day.r, (pt+1) - max(sample(1:pt, index.pp[i])))
}
day.rl <- log(day.r)
day.rl[is.infinite(day.rl)==TRUE] <- 0
X.buy[index.f, 5] <- day.rl


#1日あたりのプレイ数の対数
X.buy[index.f, 6] <- X.login[index.f, 4]

#URの確率アップキャラ
#確率アップ日とそうではない日を特定
cycle <- pt/15
r.no <- matrix(0, nrow=cycle, 5)
r.new1 <- matrix(0, nrow=cycle, 5)
r.new2 <- matrix(0, nrow=cycle, 5)

for(c in 1:cycle){
  r.no[c, ] <- (c-1) * 10 + ((c-1)*5+1):((c-1)*5+5)
  r.new1[c, ] <- (c-1) * 10 + ((c-1)*5+6):((c-1)*5+10)
  r.new2[c, ] <- (c-1) * 10 + ((c-1)*5+11):((c-1)*5+15)
}
r.cycle <- list(r.no, r.new1, r.new2)


#UR確率アップキャラを発生させる
#UR確率アップキャラのあらかじめ発生させておく
UR <- matrix(0, nrow=cycle*(length(r.cycle)-1), 9)

ur <- rep(0, length(beta.s6))
index.ur <- sample(1:length(beta.s6), 2)
ur[index.ur] <- 1
UR[1, ] <- ur

#前回のURキャラとかぶらないようにしておく
for(c in 1:(cycle*(length(r.cycle)-1))){
  for(i in 1:500){
    ur <- rep(0, length(beta.s6))
    index.ur <- sample(1:length(beta.s6), 2)
    ur[index.ur] <- 1
    UR[c, ] <- ur
    if(max(UR[c, ]+UR[c-1, ]) == 1) {break}  
  }
}


#データセットにURキャラを記録
for(t in 1:cycle){
  for(c in 1:length(r.cycle)){
    if(c==1) {next} else {
      tf <- r.cycle[[c]][t, ]
    }  
    #URキャラを記録
    X.buy[ID[, 3] %in% tf, 6:ncol(X.buy)] <- matrix(UR[2*t-2+(c-1), ], nrow=sum(ID[, 3] %in% tf), ncol=9, byrow=T)
  }
}
cbind(ID, X.buy[, 6:ncol(X.buy)])[1:120, ]   #データを確認 


##データの初期値を確認
round(cbind(ID[index.f, ], X.login[index.f, ]), 3)   #ログインデータ
round(cbind(ID[index.f, ], X.buy[index.f, ]), 3)   #ガチャデータ


####パネルごとに逐次的にデータを発生させる####
##データ更新用の保存配列
login.hist <- cbind(login.hist, matrix(0, nrow=hh, ncol=pt))
play.M <- cbind(play.M, matrix(0, nrow=hh, ncol=pt))
buy.hist <- cbind(buy.hist, matrix(0, nrow=hh, ncol=pt))


##ロジットの定義
#ガチャのロジット
logit.g <- beta.s0 + X.buy[ID[, 2]==1 & ID[, 3]==1, ] %*% beta.s   

#ログサム変数の定義
logsum <- log(1+exp(logit.g))
X.login[ID[, 2]==1 & ID[, 3]==1, 7] <- logsum

#ログインのロジット
logit.l <- alpha.s0 + X.login[ID[, 2]==1 & ID[, 3]==1, ] %*% alpha.s   #ログインのロジット

##確率を計算して応答変数を発生
Y <- matrix(0, nrow=hhpt, ncol=2)

#ログインの発生
Pr.l <- exp(logit.l)/(1+exp(logit.l))   #確率を計算
Y[ID[, 2]==1 & ID[, 3]==1, 1] <- rbinom(1, 1, Pr.l)

#ガチャの発生
#ログインがあった場合にのみ発生
Pr.g <- exp(logit.g)/(1+exp(logit.g))   #確率を計算
if(Y[ID[, 2]==1 & ID[, 3]==1, 1]==1) {
  Y[ID[, 2]==1 & ID[, 3]==1, 2] <- rbinom(1, 1, Pr.g)} else {
    Y[ID[, 2]==1 & ID[, 3]==1, 2] <- 0
}



##データの更新
#ログイン履歴の更新


