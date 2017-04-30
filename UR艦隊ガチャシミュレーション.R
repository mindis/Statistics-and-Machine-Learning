#####UR艦隊完成までの確率モデル#####
library(MASS)
library(reshape2)
library(plyr)

####シミュレーションモデルでUR艦隊に必要な期待金額を計算する####
##データの設定
UR <- 10   #URのカード種類数
P.ur <- 0.01   #URが出る確率
Pi.ur <- P.ur/UR   #それぞれのURカードの出る確率
P <- c(rep(Pi.ur, UR), 1-P.ur)

#シミュレーションの設定
T <-5000   #繰り返し数
t <- 3000   #最大ガチャ回数
Gacha <- list()
finish <- c()
UR.cnt <- c()

#UR艦隊が完成するまでガチャを回し続ける処理を繰り返す
for(i in 1:T){
  GachaM <- matrix(0, T, UR+1) 
  for(j in 1:t){
    GachaM[j, ] <- t(rmultinom(1, 11, P))
    if(j==1:17) next
    UR.c <- colSums(GachaM[1:j, 1:UR])
    UR.a <- ifelse(UR.c %% 2 == 1, UR.c-1, UR.c)
    if(sum(UR.a/2) >= 9) break
  }
  UR.cnt <- c(UR.cnt, sum(colSums(GachaM[, 1:UR])))
  finish <- c(finish, nrow(GachaM[rowSums(GachaM > 0), ]))
  print(c(i, UR.cnt[i]))
}

##結果の要約
#要約統計量
summary(UR.cnt)   #出たUR枚数の要約集計
table(UR.cnt)   #出たUR枚数の単純集計

summary(finish)   #完成までにかかるガチャ回数の要約集計
table(finish)   #完成までにかかるガチャ回数の単純集計
summary(finish)*3000   #かかった金額

#結果をプロット
hist(finish, col="grey")
hist(UR.cnt, col="grey")


####シミュレーション結果を用いてロジスティック回帰モデルで期待確率を近似する####
##ガチャ回数ごとにUR艦隊が完成している確率を計算する
max(finish) 
GC <- round(seq(1, max(finish), length=max(finish)), 0)   #検証するガチャ回数
flag <- c()


for(i in 1:length(GC)){
  flag <- c(flag, sum(ifelse(finish <= GC[i], 1, 0)))
}
Pr <- flag / T   #完成確率
round(Pr.c <- cbind(GC, Pr, flag, T), 3)   #ガチャ回数ごとの完成確率
Pr.min <- which.max(Pr.c[, "flag"] > 0)   #0％以上の最小値
plot(Pr.c[, 1:2], type="l", lwd=2)

#ロジスティック回帰モデルを当てはめる
fit <- glm(cbind(flag, T-flag) ~ GC, family="binomial", data=as.data.frame(Pr.c[Pr.min:nrow(Pr.c), ]))
summary(fit)

#結果を要約しプロット
logit <- fit$coef[1] + fit$coef[2]*1:max(finish)
Pr <- exp(logit) / (1 + exp(logit))
round(Pr, 3)

which.max(Pr>=0.5)   #50％の確率でUR艦隊が完成出来るガチャ回数
which.max(Pr>=0.5) * 3000   #かかる金額
which.max(Pr>=0.9)   #90％の確率でUR艦隊が完成出来るガチャ回数
which.max(Pr>=0.9) * 3000   #かかる金額

plot(1:length(Pr), Pr, type="l", lwd=2, xlab="ガチャ回数", main="UR艦隊完成率")   #結果をプロット
lines(Pr.c[, 1:2], lty=2)

####URカード種類を変化させてシミュレーションを実行####
##データの設定
UR <- seq(10, 30, length=11)   #検討するURのカード種類数
P.ur <- 0.01   #URが出る確率
Pi.ur <- P.ur/UR   #それぞれのURカードの出る確率

#シミュレーションの設定
T <-2000   #繰り返し数
t <- 3000   #最大ガチャ回数
Gacha <- list()
finish <- c()
UR.cnt <- c()

UR.cntL <- list()
finish.L <- list()
#UR艦隊が完成するまでガチャを回し続ける処理を繰り返す
for(u in 1:length(UR)){
  Pur <- c(rep(P.ur/UR[u], UR[u]), 1-P.ur)
  UR.cnt <- c()
  finish <- c()
  for(i in 1:T){
    GachaM <- matrix(0, T, UR[u]+1) 
    for(j in 1:t){
      GachaM[j, ] <- t(rmultinom(1, 11, Pur))
      if(j==1:17) next
      UR.c <- colSums(GachaM[1:j, 1:UR[u]])
      UR.a <- ifelse(UR.c %% 2 == 1, UR.c-1, UR.c)
      if(sum(UR.a/2) >= 9) break
    }
    UR.cnt <- c(UR.cnt, sum(colSums(GachaM[, 1:UR[u]])))
    finish <- c(finish, nrow(GachaM[rowSums(GachaM > 0), ]))
    print(c(i, UR.cnt[i]))
  }
  UR.cntL[[u]] <- UR.cnt 
  finish.L[[u]] <- finish 
  print(U)
}
hist(finish.L[[1]], xlim=c(50, 600), ylim=c(0, 1000))
par(new=T)
hist(finish.L[[11]], col=2, xlim=c(50, 600), ylim=c(0, 1000))
