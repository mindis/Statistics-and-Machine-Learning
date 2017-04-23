#####ノンパラメトリック回帰#####
library(MASS)
library(fda)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(324)
##データの設定
x <- c(1:(52*3))   #半月単位の3年

##出力データを発生
tb <- 2.0   #初期値
trend <- numeric()
s <- seq(0.85, 0.2, length=length(x))
for(i in 1:length(x)){
  r <- rnorm(5, tb, 0.02)
  sort <- sort(r)
  bi <- rbinom(1, 1, s[i])
  bb <- ifelse(bi == 1, sort[4], sort[2])
  tb <- bb
  trend <- c(trend, bb)
}
plot(trend, lwd=1, xlab="half_month", ylab="p")

####3次スプライン基底####
##節点を決める
t <- seq(x[1], x[length(x)], length=16)   #節点

##節点から局所的な多項式を作る
b <- create.bspline.basis(rangeval=c(t[1], t[length(t)]), breaks=t)
plot(b)

#平滑化を実行
sm <- smooth.basis(x, trend, b); sm   #平滑化を実行

##結果をプロット
plot(trend, lwd=1, xlab="half_month", ylab="p")
lines(sm, lwd=2, col="red")


