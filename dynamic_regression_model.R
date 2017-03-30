#####動的回帰モデル#####
library(MASS)
library(dml)
library(KFAS)
library(reshape2)
library(dplyr)
####データの発生####
#set.seed(54390)
##デートとパラメーターの初期値の設定
n <- 700   #時点数
b0 <- 6.2   #ベース販売力の初期値
b1 <- -0.01   #価格弾力性の係数
b2 <- 0.3   #特別陳列の係数
b3 <- 0.2   #特別キャンペーンの係数
b4 <- 0.007   #価格ゲインの係数
b5 <- -0.017   #価格ロスの係数

#トレンドのシミュレーションデータの発生
tb <- b0
trend <- numeric()
s <- seq(0.8, 0.1, length=n)
for(i in 1:n){
  r <- rnorm(5, tb, 0.01)
  sort <- sort(r)
  bi <- rbinom(1, 1, s[i])
  bb <- ifelse(bi == 1, sort[4], sort[2])
  tb <- bb
  trend <- c(trend, bb)
}
plot(trend, type="l", lwd=1, ylim=c(4.0, 6.0), xlab="day")
summary(trend)

##説明変数のシミュレーションデータを発生
PRICE <- numeric()
DISP <- numeric()
CAMP <- numeric()
p <- seq(0.9, 0.2, length=700)   #価格の割引確率
for(i in 1:n){
  rn <- runif(2)   #一様乱数を発生
  
  #価格の設定
  if(rbinom(1, 1, p[i])==1) SP <- 108 else SP <- 108 * runif(1, 0.7, 0.95)
  PRICE <- c(PRICE, SP)
  
  #確率0.25で特別陳列あり
  DISP <- c(DISP, (rn[1] > 0.75))
  
  #確率0.15でキャンペーンあり
  CAMP <- c(CAMP, rn[2] > 0.85)
}
(X <- data.frame(PRICE, DISP, CAMP))
summary(X)

##参照価格のシミュレーションデータの発生
#トレンドのシミュレーションデータの発生
pi <- 0.95   #参照価格の繰越パラメータ
pp <- 115   #参照価格の初期値
rp <- numeric()
for(i in 1:n){
  if(i==1) rp <- c(rp, pp) else
  rp <- c(rp, (1-pi)*PRICE[i] + pi*rp[i-1])
}
summary(rp)
plot(rp, type="l", lwd=1, xlab="day", ylab="参照価格", ylim=c(75, 130))   #参照価格のプロット

##ゲイン変数とロス変数の設定
GL <- ifelse(rp-PRICE > 0, 1, 0)   #ゲイン・ロス指示変数
refPRICE <- rp   #参照価格
GLvalue <- abs(rp-PRICE)   #参照価格と価格との差

##販売数量を発生
#対数変換時のyの販売量
(y <- trend + b1*log(PRICE) + b2*DISP + b3*CAMP + Z*b4*log(GLvalue) + (1-Z)*b5*log(GLvalue) + rnorm(n, 0, 0.2))
yy <- round(exp(y), 0)
max(yy)
min(yy)
summary(yy)

#販売数量の時系列をプロット
plot(1:n, yy, type="l", xlab="day", ylab="販売数量")
lines(1:n, exp(trend), lwd=2)
 
##データをすべて結合
YX <- data.frame(yy, X, GL, refPRICE)
round(YX, 0)
