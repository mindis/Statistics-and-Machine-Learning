#####多階層ネステッドロジットモデル#####
library(mlogit)
library(nnet)
library(MASS)
library(dplyr)
library(reshape2)
library(ggplot2)

#set.seed(31489)

####データの発生####
N <- 3000
choise <- 10   #選択肢数
k <- 3   #階層数
n1 <- 2   #ネスト1のネスト数
n2 <- 4   #ネスト2のネスト数

####説明変数の発生####
Gacha <- matrix(0, nrow=N, ncol=choise)
Prom <- matrix(0, nrow=N, ncol=choise)
NPS <- matrix(0, nrow=N, ncol=choise)
Roy <- matrix(0, nrow=N, ncol=choise)

for(i in 1:choise){
  p <- runif(2, 0.3, 0.5)
  Gacha[, i] <- rbinom(N, 1, p[1])
  Prom[, i] <- rbinom(N, 1, p[2])
  NPS[, i] <- rnorm(N, runif(1, 4, 8), runif(1, 1, 2))
  Roy[, i] <- rnorm(N, 0, 1)
}

#NPSを-5～5の範囲に収める
NPS <- round(ifelse(NPS > 10, 10, ifelse(NPS < 1, 1, NPS)), 0) - 5

##説明変数をベクトル変換する
#IDの設定
u.id <- rep(1:N, rep(choise, N))
c.id <- rep(1:choise, N)
ID <- data.frame(no=1:length(id), u.id=u.id, c.id=c.id)

#切片の設定
Value <- matrix(as.numeric(diag(choise)), nrow=N*choise, ncol=choise, byrow=T)[, -10]

#ベクトル変換
Gacha_vec <- as.numeric(t(Gacha))
Prom_vec <- as.numeric(t(Prom))
NPS_vec <- as.numeric(t(NPS))
Roy_vec <- as.numeric(t(Roy))

#データを結合
X <- data.frame(game=Value, Gacha=Gacha_vec, prom=Prom_vec, NPS=NPS_vec, Roy=Roy_vec)
XM <- as.matrix(X)


####応答変数の発生####
##ネスト構造の設定
nest01 <- c(rep(1, choise/2), rep(2, choise/2)) 
nest1 <- c(rep(1, 2), rep(2, 2))
nest2 <- c(rep(1, 2), rep(2, 3), rep(3, 3), rep(4, 2))

#ネスト行列を設定
nest1z <- matrix(0, nrow=n2, ncol=n1)
nest2z <- matrix(0, nrow=choise, ncol=n2)
for(i in 1:n1) {nest1z[nest1==i, i] <- 1}
for(i in 1:n2) {nest2z[nest2==i, i] <- 1}

##パラメータの設定
#相関パラメータの設定
rho01 <- 1/c(0.25, 0.6)
rho02 <- 1/c(0.3, 0.7, 0.2, 0.6)

#回帰パラメータの設定
beta00 <- runif(choise-1, -1.5, 2.4)
beta01 <- runif(1, 0.4, 1.0)
beta02 <- runif(1, 0.3, 0.9)
beta03 <- runif(1, 0, 0.15)
beta04 <- runif(1, 0, 0.4)
beta0 <- c(beta00, beta01, beta02, beta03, beta04)
#beta0[1:length(beta0)] <- 0

##効用関数とログサム変数の定義
#効用関数の定義
U <- matrix(XM %*% beta0, nrow=N, ncol=choise, byrow=T)

rho11 <- 1/matrix(c(0.25, 0.25, 0.6, 0.6), nrow=N, ncol=n2, byrow=T)

#ログサム変数の定義
#階層2のログサム変数
U_logm2 <- exp(U * matrix(as.numeric(rho02 %*% t(nest2z)), nrow=N, ncol=choise, byrow=T))
logsum2 <- log(U_logm2 %*% nest2z)

#クラスターごとの選択確率
rho22 <- rho21 / matrix(rho02, nrow=N, ncol=n2, byrow=T)
logsum2u <- exp(rho22 * logsum2)
CL2 <- logsum2u / matrix(rowSums(logsum2u), nrow=N, ncol=n2)


#階層1のログサム変数
U_logm1 <- exp(logsum2u / matrix(as.numeric(rho01 %*% t(nest01z)), nrow=N, ncol=n2, byrow=T))
logsum1 <- log(logsum2u %*% nest1z)

#クラスターごとの選択確率
rho11 <- 1/matrix(rho01, nrow=N, ncol=n1, byrow=T)
logsum1u <- exp(rho11 * logsum1)
CL1 <- logsum1u / matrix(rowSums(logsum1u), nrow=N, ncol=n1)


#最下層の確率を計算
U1 <- exp(U * matrix(as.numeric(rho02 %*% t(nest2z)), nrow=N, ncol=choise, byrow=T))
U2 <- U1 %*% nest2z
 
#ゲーム選択確率を計算
U_sums <- matrix(0, nrow=N, ncol=choise)
Pr <- matrix(0, nrow=N, ncol=choise)

for(i in 1:choise){
  U_sums[, i] <- U2[, nest2[i]]
  Pr[, i] <- CL1[, nest01[i]] * CL2[, nest2[i]] * U1[, i]/U_sums[, i]   #条件付き確率を計算
}
nest01
nest2

rowSums(Pr)

