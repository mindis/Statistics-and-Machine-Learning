#####離脱の潜在変数を含む階層ベイズ多変量離散時間ハザードモデル#####
library(MASS)
library(nlme)
library(glmm)
library(suvival)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

set.seed(9483)

####データの発生####
hh <- 1000   #サンプル数
pt <- 36   #観測期間
m <- 9

##IDの設定
u.id <- rep(1:hh, rep(pt, hh))
t.id <- rep(1:pt, hh)
ID <- data.frame(no=1:(hh*pt), id=u.id, time=t.id)

####説明変数の発生####
##階層モデルの説明変数
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


##個体内モデルの説明変数の発生
#勧誘できるURを月ごとに特定
ur <- matrix(0, nrow=pt, ncol=m)

for(i in 1:pt){
  for(j in 1:1000){
    ur[i, ] <- t(rmultinom(1, 2, rep(1/m, m)))
    if(i==1){ break
    } else {
      if(max(colSums(ur[(i-1):i, ]))==1) break
    }
  }
}

#全データ分のデータを発生
UR <- matrix(as.numeric(t(UR)), nrow=pt*hh, ncol=m, byrow=T)
colnames(UR) <- c("honoka", "kotori", "umi", "rin", "hanayo", "maki", "nico", "nozomi", "eri")
UR <- UR[, -1]   #基準メンバーを削除


#勧誘優待の有無を発生
camp <- rbinom(pt, 1, runif(1, 0.4, 0.55))
Camp <- rep(camp, hh)


####応答変数を発生####
##トレンド成分を発生
T <- 10000
for(t1 in 1:T){
  tb <- -1.0   #初期値
  trend <- numeric()
  s <- seq(0.85, 0.2, length=t)
  for(i in 1:pt){
    r <- rnorm(5, tb, 0.05)
    sort <- sort(r)
    bi <- rbinom(1, 1, s[i])
    bb <- ifelse(bi == 1, sort[4], sort[2])
    tb <- bb
    trend <- c(trend, bb)
  }
  if(max(pnorm(trend)) < 0.25 & min(pnorm(trend)) > 0.15) break
  print(t1)
}  
plot(pnorm(trend), type="l", lwd=1, xlab="月", ylab="p")
pnorm(trend)


##階層モデルのパラメータを設定
#回帰パラメータを設定
theta0 <- c(runif(1, -0.5, 0.5), runif(cont1, 0, 0.6), runif(bin1+multi1-1, -0.6, 0.6))
theta1 <- rbind(c(0.8, 0.2, -0.6, -0.4, 1.0, 0.6, -0.8, 0.5), matrix(runif(cont1*(m-1), 0, 0.6), nrow=cont1, ncol=m-1), 
                matrix(runif((bin1+multi1-1)*(m-1), -0.6, 0.7), nrow=bin1+multi1-1, ncol=m-1))
theta2 <- c(runif(1, 0.3, 0.8), runif(cont1, 0, 0.6), runif(bin1+multi1-1, -0.5, 0.8))
theta3 <- c(runif(1, -0.9, -0.4), runif(cont1, 0, 0.6), runif(bin1+multi1-1, -0.5, 0.8))
theta4 <- c(runif(1, 0.2, 0.5), runif(cont1, 0, 0.4), runif(bin1+multi1-1, -0.4, 0.5))
theta_T <- cbind(theta0, theta1, theta2, theta3, theta4)

#分散共分散行列を設定
Cov0 <- diag(runif(ncol(theta_T), 0.4, 0.8))

#多変量正規分布よりパラメータを発生
BETA0 <- ZX %*% theta_T + mvrnorm(hh, rep(0, ncol(theta_T)), Cov0)
BETA_T <- BETA0


##応答変数を発生させながら時間依存変数を発生させる
