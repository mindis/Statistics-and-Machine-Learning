#####動的離散選択モデル#####
library(MASS)
library(Matrix)
library(matrixStats)
library(RcppSMC)
library(SMC)
library(dml)
library(KFAS)
library(extraDistr)
library(reshape2)
library(dplyr)

####データの発生####
#set.seed(54390)

##データの設定
n <- 1000   #観測期間
time <- 1:n   #観測期間id
k <- 4   #説明変数数


####説明変数を発生####
#価格を生成
PRICE <- runif(n, 0, 1.0)   

#プロモーション有無を生成
DISP <- rbinom(n, 1, 0.4)   #特別陳列有無
AD <- rbinom(n, 1, 0.4)   #チラシ掲載有無

#データを結合
Data <- cbind(PRICE, DISP, AD)


####パラメータを生成し、応答変数を決定####
##トレンド成分の生成
b0 <- -0.8   #トレンドの初期値
v0 <- 0.015   #システムモデルの分散

#時間ごとに逐次的にトレンド成分を生成
trend <- rep(0, n)
trend[1] <- b0
s <- seq(0.7, 0.3, length=n-1)   
for(i in 2:n){
  diff <- rnorm(5, 0, v0)   #変化の候補を生成
  sortlist <- sort(diff)   #昇順に並び替える
  bi <- rbinom(1, 1, s[i-1])   #変化の仕方を決定
  trend[i] <- trend[i-1] + bi*sortlist[4] + (1-bi)*sortlist[2]
}
plot(1:n, trend, type="l", xlab="観測期間")
summary(trend); trend0 <- trend

##説明変数の動的パラメータを生成
#初期値を設定
beta1 <- beta2 <- beta3 <- rep(0, n)
beta1[1] <- -0.9   #価格の初期値
beta2[1] <- 0.8   #特別陳列の初期値
beta3[1] <- 0.7   #チラシ掲載有無の初期値 

#システムモデルの分散
v1 <- 0.015   
v2 <- 0.015
v3 <- 0.015

#時間ごとに逐次的に動的パラメータを生成
s1 <- seq(0.6, 0.4, length=n-1)
s2 <- seq(0.40, 0.7, length=n-1)
s3 <- seq(0.45, 0.6, length=n-1)
for(i in 2:n){
  diff1 <- rnorm(5, 0, v1); diff2 <- rnorm(5, 0, v2); diff3 <- rnorm(5, 0, v3)
  sortlist1 <- sort(diff1); sortlist2 <- sort(diff2); sortlist3 <- sort(diff3)
  bi1 <- rbinom(1, 1, s1[i-1]); bi2 <- rbinom(1, 1, s2[i-1]); bi3 <- rbinom(1, 1, s3[i-1])
  beta1[i] <- beta1[i-1] + bi1*sortlist1[2] + (1-bi1)*sortlist1[4]
  beta2[i] <- beta2[i-1] + bi2*sortlist2[4] + (1-bi2)*sortlist2[2]
  beta3[i] <- beta3[i-1] + bi3*sortlist3[4] + (1-bi3)*sortlist3[2]
}
plot(1:n, beta1, type="l", xlab="観測期間")
plot(1:n, beta2, type="l", xlab="観測期間")
plot(1:n, beta3, type="l", xlab="観測期間")

#パラメータを結合
beta0 <- cbind(beta1, beta2, beta3)


##動的ロジスティック回帰モデルから応答変数を生成
#ロジットと確率の定義
logit <- trend + rowSums(Data * beta0)   #ロジット
Pr0 <- exp(logit)/(1+exp(logit))   #確率

#ベルヌーイ分布から応答変数を生成
y <- rbinom(n, 1, Pr0)

#発生させた応答変数をプロット
plot(1:n, y, xlab="観測期間", ylab="購買有無", main="購買有無と購買確率の関連")
par(new=T)
plot(1:n, Pr0, xlim=c(0, n), ylim=c(0, 1), ylab="", xlab="", type="b", pch=4)


####粒子フィルタで動的ロジスティック回帰モデルを推定####
##ロジスティック回帰モデルの対数尤度関数(初期値推定数)
fr <- function(beta, X, y){
  
  #ロジットと確率を定義
  logit <- X %*% beta
  Pr <- exp(logit)/(1+exp(logit))
  
  #対数尤度関数の和
  Li <- y*log(Pr) + (1-y)*log(1-Pr)
  LL <- sum(Li)
  return(LL)
}

####粒子フィルタで動的パラメータを推定####
##粒子フィルタの設定
s <- 10000   #粒子数
tau <- diag(0.05^2, k)   #システム分散
LL <- rep(0, n)
BETA <- array(0, dim=c(s, k, n))

##システムモデルのパラメータの更新
betan <- mvrnorm(s, rep(0, k), diag(0.25, k))

##観測モデルの尤度を評価
#ロジットと確率を計算
logit <- betan[, 1] + rowSums(matrix(Data[1, ], nrow=s, ncol=k-1, byrow=T) * betan[, -1])
Pr <- exp(logit)/(1+exp(logit))

#尤度を評価
Li <- Pr^y[1] * (1-Pr)^(1-y[1])   #粒子ごとの尤度
LL[1] <- sum(Li)   #尤度の和

#尤度の負担率に応じてパラメータをリサンプリング
w <- Li/sum(Li) 
index <- as.numeric(rmnom(1, s, w))
resample <- rep(1:s, index)
BETA[, , 1] <- betan[resample, ]   #リサンプリングされたパラメータ


##2期目以降を粒子フィルタで逐次的に更新
for(i in 2:n){
  ##システムモデルのパラメータの更新
  betan <- BETA[, , i-1] + mvrnorm(s, rep(0, k), tau)
  
  ##観測モデルの尤度を評価
  #ロジットと確率を計算
  logit <- betan[, 1] + rowSums(matrix(Data[i, ], nrow=s, ncol=k-1, byrow=T) * betan[, -1])
  Pr <- exp(logit)/(1+exp(logit))
  
  #尤度を評価
  Li <- Pr^y[i] * (1-Pr)^(1-y[i])   #粒子ごとの尤度
  LL[i] <- sum(Li)   #尤度の和
  
  #尤度の負担率に応じてパラメータをリサンプリング
  w <- Li/sum(Li) 
  index <- as.numeric(rmnom(1, s, w))
  resample <- rep(1:s, index)
  BETA[, , i] <- betan[resample, ]   #パラメータをリサンプリング
}
#対数尤度の和
LLs <- sum(log(LL)) - n*log(s)


##固定ラグ平滑化により最終的なパラメータを確定
L <- 20   #ラグ数
LAG <- array(0, dim=c(s, k, L))
THETA <- BETA

for(i in L:n){
  print(i)
  
  lag <- i-L+1
  LAG <- THETA[, , lag:i]
  betan <- LAG[, , L-1] + mvrnorm(s, rep(0, k), tau)
  
  ##観測モデルの尤度を評価
  #ロジットと確率を計算
  logit <- betan[, 1] + rowSums(matrix(Data[i, ], nrow=s, ncol=k-1, byrow=T) * betan[, -1])
  Pr <- exp(logit)/(1+exp(logit))
  
  #尤度を評価
  Li <- Pr^y[i] * (1-Pr)^(1-y[i])   #粒子ごとの尤度
  
  #尤度の負担率に応じてパラメータをリサンプリング
  w <- Li/sum(Li) 
  index <- as.numeric(rmnom(1, s, w))
  resample <- rep(1:s, index)
  THETA[, , lag:i] <- LAG[resample, , ]   #パラメータをリサンプリング
}


####推定結果の確認と可視化####
##事後平均を確認
theta1 <- t(apply(BETA, c(2, 3), mean))
theta2 <- t(apply(THETA, c(2, 3), mean))
round(cbind(theta1, theta2, trend, beta0), 3)

#フィルタリング、平滑化、真値を可視化
matplot(theta1, type="l", xlab="観測期間", ylab="パラメータ")
matplot(theta2, type="l", xlab="観測期間", ylab="パラメータ")
matplot(cbind(trend, beta1, beta2, beta3), type="l", xlab="観測期間", ylab="パラメータ")

##確率と対数尤度を再計算
logit <- trend + rowSums(Data * beta0)
Pr1 <- exp(logit)/(1+exp(logit))
Li1 <- Pr1^y * (1-Pr1)^(1-y)   #粒子ごとの尤度
LL1 <- sum(log(Li1))

logit <- theta1[, 1] + rowSums(Data * theta1[, -1])
Pr2 <- exp(logit)/(1+exp(logit))
Li2 <- Pr2^y * (1-Pr2)^(1-y)   #粒子ごとの尤度
LL2 <- sum(log(Li2))

logit <- theta2[, 1] + rowSums(Data * theta2[, -1])
Pr3 <- exp(logit)/(1+exp(logit))
Li3 <- Pr3^y * (1-Pr3)^(1-y)   #粒子ごとの尤度
LL3 <- sum(log(Li3))
c(LL1, LL2, LL3)
round(cbind(y, Pr1, Pr2, Pr3), 3)
c(LL1, LL2, LL3)

##事後信用区間を推定
theta_lower <- matrix(0, nrow=n, ncol=k)
theta_upper <- matrix(0, nrow=n, ncol=k)
for(i in 1:n){
  theta_lower[i, ] <- apply(THETA[, , i], 2, function(x) quantile(x, 0.05))
  theta_upper[i, ] <- apply(THETA[, , i], 2, function(x) quantile(x, 0.95))
}


##トレンド成分を現系列にプロット
plot(1:n, y, xlab="観測期間", ylab="購買有無", main="購買有無と購買確率の関連")
par(new=T)
plot(1:n, Pr0, xlim=c(0, n), ylim=c(0, 1), ylab="", xlab="", type="b", pch=4)
par(new=T)
plot(1:n, exp(theta2[, 1])/(1+exp(theta2[, 1])), xlim=c(0, n), ylim=c(0, 1), 
     ylab="", xlab="", type="l", lwd=2, col=2)
par(new=T)
plot(1:n, exp(theta_lower[, 1])/(1+exp(theta_lower[, 1])), xlim=c(0, n), ylim=c(0, 1), 
     ylab="", xlab="", type="l", lwd=2, col=3)
par(new=T)
plot(1:n, exp(theta_upper[, 1])/(1+exp(theta_upper[, 1])), xlim=c(0, n), ylim=c(0, 1), 
     ylab="", xlab="", type="l", lwd=2, col=3)
par(new=T)
plot(1:n, exp(theta1[, 1])/(1+exp(theta1[, 1])), xlim=c(0, n), ylim=c(0, 1), 
     ylab="", xlab="", type="l", lwd=2, col=4)
par(new=T)
plot(1:n, exp(trend)/(1+exp(trend)), xlim=c(0, n), ylim=c(0, 1), ylab="", xlab="", 
     type="l", lwd=2, col=5)

