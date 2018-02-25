#####動的ポアソン回帰モデル#####
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

#set.seed(52698)

####データの発生####
##データの設定
d <- 2000   #観測期間
time <- 1:d   #観測期間id
k <- 4   #説明変数数

##時間ごとに逐次的に説明変数を生成
for(rp in 1:1000){
  #説明変数のトレンドを生成
  Data0 <- matrix(0, nrow=d, ncol=k)
  Data0[1, ] <- c(1.5, 0.6, -0.5, -0.5)
  v0 <- runif(k, 0.015, 0.03)   #システムモデルの分散
  s <- cbind(seq(0.8, 0.2, length=d-1), 
             seq(0.7, 0.4, length=d-1),
             seq(0.4, 0.7, length=d-1),
             seq(0.4, 0.7, length=d-1))
  
  for(i in 2:d){
    for(j in 1:k){
      diff <- rnorm(5, 0, v0[j])   #変化の候補を生成
      sortlist <- sort(diff)   #昇順に並び替える
      bi <- rbinom(1, 1, s[i-1, j])   #変化の仕方を決定
      Data0[i, j] <- Data0[i-1, j] + bi*sortlist[4] + (1-bi)*sortlist[2]
    }
  }
  Data1 <- Data0
  matplot(Data1, type="l", xlab="日数", ylab="説明変数の変動")
  
  #価格の割引率を設定
  x <- rnorm(d, 0.9, 0.15)
  Data1[, 2] <- Data0[, 2] * ifelse(x > 1, 1, x)
  
  #二値変数を生成
  u <- t(apply(Data0[, 3:k], 1, function(x) mvrnorm(1, x, diag(1, length(3:k)))))
  Data1[, 3:k] <- matrix(as.numeric(u > 0), nrow=d, ncol=length(3:k))
  Data <- cbind(1, Data1[, -1])
  
  
  ##説明変数の動的パラメータを生成
  #初期値を設定
  beta1 <- beta2 <- beta3 <- rep(0, d)
  beta1[1] <- -1.25   #価格の初期値
  beta2[1] <- 0.6   #特別陳列の初期値
  beta3[1] <- 0.6   #チラシ掲載有無の初期値 
  
  #システムモデルの分散
  v1 <- v2 <- v3 <- 0.015   
  
  #時間ごとに逐次的に動的パラメータを生成
  s1 <- seq(0.6, 0.4, length=d-1)
  s2 <- seq(0.4, 0.7, length=d-1)
  s3 <- seq(0.4, 0.6, length=d-1)
  for(i in 2:d){
    diff1 <- rnorm(5, 0, v1); diff2 <- rnorm(5, 0, v2); diff3 <- rnorm(5, 0, v3)
    sortlist1 <- sort(diff1); sortlist2 <- sort(diff2); sortlist3 <- sort(diff3)
    bi1 <- rbinom(1, 1, s1[i-1]); bi2 <- rbinom(1, 1, s2[i-1]); bi3 <- rbinom(1, 1, s3[i-1])
    beta1[i] <- beta1[i-1] + bi1*sortlist1[2] + (1-bi1)*sortlist1[4]
    beta2[i] <- beta2[i-1] + bi2*sortlist2[4] + (1-bi2)*sortlist2[2]
    beta3[i] <- beta3[i-1] + bi3*sortlist3[4] + (1-bi3)*sortlist3[2]
  }
  plot(1:d, beta1, type="l", xlab="観測期間")
  plot(1:d, beta2, type="l", xlab="観測期間")
  plot(1:d, beta3, type="l", xlab="観測期間")
  
  #パラメータを結合
  beta <- betat <- cbind(beta1, beta2, beta3)
  trend <- Data1[, 1]
  
  
  ##動的ポアソン回帰モデルから応答変数を生成
  pois_mu <- exp(rowSums(Data * cbind(trend, beta)))   #ポアソン分布の平均
  y <- rpois(d, pois_mu)   #ポアソン分布から応答変数を生成
  if(max(y[1:(d-500)]) > max(y[(d-500):d]) & max(y) > 30){
    break
  }
}
plot(1:d, y, type="l", xlab="日数", ylab="購買点数", xlim=c(0, d), ylim=c(0, max(y)))
par(new=T)
plot(1:d, pois_mu, xlim=c(0, d), ylim=c(0, max(y)), ylab="", xlab="", type="p", pch=4, col=4)


####粒子フィルタで動的ポアソン回帰モデルを推定####
##ポアソン回帰モデルの対数尤度関数
fr <- function(beta, Data, y, y_factorial){
  lambda <- exp(Data %*% beta)   #リンク関数
  LLi <- y*log(lambda)-lambda - y_factorial
  LL <- sum(LLi)
  return(LL)
}

##粒子フィルタの設定
s <- 10000   #粒子数
tau <- diag(0.05^2, k)   #システム分散
LL <- rep(0, d)
BETA <- array(0, dim=c(s, k, d))
y_factorial <- lfactorial(y)

##システムモデルのパラメータを更新
betan <- c(trend[1], beta[1, ]) + mvrnorm(s, rep(0, k), diag(0.25, k))


##観測モデルの尤度を評価
#尤度を評価
lambda <- exp(betan[, 1] + rowSums(matrix(Data[1, ], nrow=s, ncol=k-1, byrow=T) * betan[, -1]))
Li <- dpois(y[1], lambda)   #粒子ごとの尤度
LL[1] <- sum(Li)   #尤度の和


#尤度の負担率に応じてパラメータをリサンプリング
w <- Li/sum(Li) 
index <- as.numeric(rmnom(1, s, w))
resample <- rep(1:s, index)
BETA[, , 1] <- betan[resample, ]   #リサンプリングされたパラメータ


##2期目以降を粒子フィルタで逐次的に更新
for(i in 2:d){
  ##システムモデルのパラメータの更新
  betan <- BETA[, , i-1] + mvrnorm(s, rep(0, k), tau)

  ##観測モデルの尤度を評価
  #ロジットと確率を計算
  lambda <- exp(betan[, 1] + rowSums(matrix(Data[i, ], nrow=s, ncol=k-1, byrow=T) * betan[, -1]))
  Li <- dpois(y[i], lambda)   #粒子ごとの尤度
  LL[i] <- sum(Li)   #尤度の和
  
  #尤度の負担率に応じてパラメータをリサンプリング
  w <- Li/sum(Li) 
  index <- as.numeric(rmnom(1, s, w))
  resample <- rep(1:s, index)
  BETA[, , i] <- betan[resample, ]   #パラメータをリサンプリング
}
#対数尤度の和
LLs <- sum(log(LL)) - d*log(s)
LLs
sum(dpois(y, exp(rowSums(Data * cbind(trend, beta))), log=TRUE))


i <- 100
round(colMeans(BETA[, , i]), 3)
Data[i, ] %*% c(trend[i], beta[i, ])
Data[i,] %*% colMeans(BETA[, , i])

