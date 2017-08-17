#####潜在クラスを利用したデータ融合####
library(MASS)
library(bayesm)
library(MCMCpack)
library(mlogit)
library(flexmix)
library(psych)
library(gtools)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
##データの設定
hh <- 3000
h <- 1000   #両方プレイしているユーザー数
k <- 8   #セグメント数


####共通の説明変数の発生####
cont <- 2; bin <- 3; multi <- 4
X.cont <- matrix(rnorm(hh*cont), nrow=hh, ncol=cont)
X.bin <- matrix(0, nrow=hh, ncol=bin)
X.multi <- matrix(0, nrow=hh, ncol=multi)

#二値説明変数を設定
for(i in 1:bin){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(hh, 1, p)
}

#多値説明変数を設定
p <- runif(multi)
X.multi <- t(rmultinom(hh, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))] #冗長な変数は削除 


##説明変数をベクトル形式に変更
#切片をベクトル形式に変更
bv <- c(1, rep(0, k))
iv <- matrix(bv, nrow=hh*length(bv), ncol=k, byrow=T)
IV <- subset(iv, rowSums(iv) > 0)
IV <- IV[, -k]

#説明変数をベクトル形式に変更
index.z <- rep(1:hh, rep(k, hh))
Zi <- matrix(0, nrow=hh*k, ncol=(cont+bin+multi-1)*(k-1))
Xi <- cbind(X.cont, X.bin, X.multi)

for(i in 1:hh){
  x.bind <- c()
  for(j in 1:ncol(Xi)){
    x.diag <- diag(Xi[i, j], k)
    x.bind <- cbind(x.bind, x.diag[, -k])
  }
  Zi[index.z==i, ] <- x.bind
}

#データを結合
Zx1 <- cbind(inter=IV, Z=Zi)

####欠損パターンの発生####
#欠損するかどうかはランダムな欠損を仮定
#欠損パターンは3種類
pattern <- 3
alpha0 <- runif(pattern-1, -0.6, 0.6)
alpha1 <- runif((pattern-1)*cont, 0, 0.7)
alpha2 <- runif((pattern-1)*bin, -0.8, 0.8)
alpha3 <- runif((pattern-1)*(multi-1), -0.9, 0.9)
alphat <- c(alpha0, alpha1, alpha2, alpha3)

index.pattern <- as.numeric(t(matrix(1:ncol(Zx1), ncol=k-1, byrow=T)[, 1:(pattern-1)]))

##ロジットと確率の計算
logit1 <- matrix(Zx1[, index.pattern] %*% alphat, nrow=hh, ncol=pattern, byrow=T)
Pr1 <- exp(logit1)/matrix(rowSums(exp(logit1)), nrow=hh, ncol=pattern)

##応答変数の発生
Z1 <- t(apply(Pr1, 1, function(x) rmultinom(1, 1, x)))
z1 <- Z1 %*% 1:pattern
z.vec <- z1
n_seg <- colSums(Z1)

colMeans(Z1); colSums(Z1)   #結果を集計


##欠損パターンをベクトル形式に変更
Z1.vec <- matrix(0, nrow=nrow(Zx1), ncol=(k-1)*2)

for(i in 1:hh){
  r <- ((i-1)*k+1):((i-1)*k+k)
  Z1.vec[r, ] <- cbind(diag(Z1[i, 1], k)[, -k], diag(Z1[i, 2], k)[, -k])
}

Zx2 <- cbind(Zx1, type=Z1.vec)   #データを結合


####潜在クラスを発生####
##パラメータの設定
theta0 <- runif(k-1, -0.55, 0.55)
theta1 <- runif((k-1)*cont, 0, 0.6)
theta2 <- runif((k-1)*bin, -0.7, 0.7)
theta3 <- runif((k-1)*(multi-1), -0.7, 0.7)
theta4 <- runif((k-1)*(pattern-1), -1.0, 1.0)
thetat <- c(theta0, theta1, theta2, theta3, theta4)

##ロジットと確率の計算
logit2 <- matrix(Zx2 %*% thetat, nrow=hh, ncol=k, byrow=T)
Pr2 <- exp(logit2)/matrix(rowSums(exp(logit2)), nrow=hh, ncol=k)

##応答変数の発生
Z2 <- t(apply(Pr2, 1, function(x) rmultinom(1, 1, x)))
z2 <- Z2 %*% 1:k
colMeans(Z2); colSums(Z2)   #結果を集計
 
#欠損パターンと潜在クラスのクロス票
table(z1, z2)
round(table(z1, z2)/matrix(rowSums(table(z1, z2)), nrow=pattern, ncol=k), 3)


####観測変数の発生####
#欠損パターンの1はラブライブ!のみ2はアイマスのみ3は両方観測
k1 <- 9   #ラブライブデータの変数   
k2 <- 11   #アイマスデータの変数

#頻度を発生させる
ns11 <- rpois(hh, 40)
ns12 <- rpois(hh, 40)
ns21 <- rpois(hh, 30)
ns22 <- rpois(hh, 30)

##確率ベクトルを定義して、説明変数を発生
P11 <- matrix(0, nrow=k, ncol=k1) 
P12 <- matrix(0, nrow=k, ncol=k2)
P21 <- rep(0, k)
P22 <- rep(0, k)

X11 <- matrix(0, nrow=hh, ncol=k1)
X12 <- matrix(0, nrow=hh, ncol=k2)
X21 <- rep(0, hh)
X22 <- rep(0, hh)

##潜在クラスごとに応答変数を発生
for(j in 1:k){
  r <- subset(1:length(z2), z2==j)

  #確率を計算
  p11 <- rgamma(k1, 0.85, 0.2)
  p12 <- rgamma(k2, 1.0, 0.15)
  P11[j, ] <- p11 / sum(p11)
  P12[j, ] <- p12 / sum(p12)
  
  P21[j] <- runif(1, 0.15, 0.7)
  P22[j] <- runif(1, 0.2, 0.6)
  
  #応答変数を発生
  X11[r, ] <- t(apply(cbind(ns11[r], 1), 1, function(x) rmultinom(1, x[1], P11[j, ])))
  X12[r, ] <- t(apply(cbind(ns12[r], 1), 1, function(x) rmultinom(1, x[1], P12[j, ])))
  X21[r] <- apply(cbind(ns21[r], 1), 1, function(x) rbinom(1, x[1], P21[j]))
  X22[r] <- apply(cbind(ns22[r], 1), 1, function(x) rbinom(1, x[1], P22[j]))
}

##応答変数の集計
by(X11, as.character(z2), colMeans)
by(X12, as.character(z2), colMeans)

##z1=3以外のユーザーのデータを欠損させる
X11[z1==2, ] <- 10
X12[z1==1, ] <- 10
X21[z1==2] <- 10
X22[z1==1] <- 10
ns11[z1==2] <- 10*k1
ns12[z1==1] <- 10*k2
ns21[z1==2] <- 20
ns22[z1==1] <- 20

X <- data.frame(ll1=X11, im1=X12, ll2=X21, im2=X22)
XM <- as.matrix(X)


####EMアルゴリズムで欠損データのセグメント割当確率を予測####
####EMアルゴリズムに必要な対数尤度関数を定義####
##観測データの対数尤度と潜在変数zを計算するための関数
LLobz <- function(theta11, theta12, theta21, theta22, ns11, ns12, ns21, ns22, X11, X12, X21, X22, z, r, k, n_seg){
  
  #欠損パターンごとに尤度を計算
  LLind1 <- matrix(0, nrow=n_seg[1], ncol=k)
  LLind2 <- matrix(0, nrow=n_seg[2], ncol=k)
  LLind3 <- matrix(0, nrow=n_seg[3], ncol=k)
  
  #潜在クラスごとに尤度を計算して代入
  for(i in 1:k){
    Li11 <- apply(cbind(ns11, X11), 1, function(x) dmultinom(x[-1], x[1], theta11[i, ]))   #ラブライブの多項分布の尤度
    Li12 <- apply(cbind(ns12, X12), 1, function(x) dmultinom(x[-1], x[1], theta12[i, ]))   #アイマスの多項分布の尤度
    Li21 <- apply(cbind(ns21, X21), 1, function(x) dbinom(x[-1], x[1], theta21[i]))   #ラブライブの二項分布の尤度
    Li22 <- apply(cbind(ns22, X22), 1, function(x) dbinom(x[-1], x[1], theta22[i]))   #アイマスの二項分布の尤度
    
    #尤度を計算
    Li1 <- Li11[z==1] * Li21[z==1]
    Li2 <- Li12[z==2] * Li22[z==2]
    Li3 <- Li11[z==3] * Li12[z==3] * Li21[z==3] * Li22[z==3]
    
    #尤度を潜在クラスごとに代入
    LLind1[, i] <- Li1
    LLind2[, i] <- Li2
    LLind3[, i] <- Li3
  }
  
  #観測データの尤度を計算
  LLho1 <- matrix(r[z==1, ], nrow=n_seg[1], ncol=k, byrow=T) * LLind1   #ラブライブのみの観測データの尤度
  LLho2 <- matrix(r[z==2, ], nrow=n_seg[2], ncol=k, byrow=T) * LLind2   #アイマスのみの観測データの尤度
  LLho3 <- matrix(r[z==3, ], nrow=n_seg[3], ncol=k, byrow=T) * LLind3   #両方観測の観測データの尤度
  
  
  #潜在変数zの計算
  z1 <- LLho1 / matrix(apply(LLho1, 1, sum), nrow=n_seg[1], ncol=k)  
  z2 <- LLho2 / matrix(apply(LLho2, 1, sum), nrow=n_seg[2], ncol=k) 
  z3 <- LLho3 / matrix(apply(LLho3, 1, sum), nrow=n_seg[3], ncol=k) 
  
  #全データ行列に潜在変数zを代入
  Z <- matrix(0, nrow=sum(n_seg), ncol=k)
  Z[z==1, ] <- z1
  Z[z==2, ] <- z2
  Z[z==3, ] <- z3
  
  
  #観測データの対数尤度の和を計算
  LLosum1 <- sum(log(apply(r[z==1, ] * LLind1, 1, sum)))   
  LLosum2 <- sum(log(apply(r[z==2, ] * LLind2, 1, sum)))  
  LLosum3 <- sum(log(apply(r[z==3, ] * LLind3, 1, sum)))  
  LLosum <- LLosum1 + LLosum2 + LLosum3

  rval <- list(LLob=LLosum, Z=Z)
  return(rval)
}

##多項ロジットモデルの対数尤度関数
loglike <- function(x, Z, z, hh, seg){
  #パラメータの設定
  theta.z <- x
  
  #効用関数の設定
  U <- matrix(Z %*% theta.z, nrow=hh, ncol=seg, byrow=T)
  
  #対数尤度の計算
  d <- rowSums(exp(U))
  LLl <- rowSums(z1 * U) - log(d)
  LL <- sum(LLl)
  return(LL)
}

####EMアルゴリズムの設定と初期値の設定####
##EMアルゴリズムの設定
iter <- 0
dl <- 100
tol <- 1


##初期値の設定
#ベストな初期パラメータを選択
rp <- 100   #繰り返し数
Z.list <- list()
val <- c()
x <- list()
theta.x <- matrix(0, nrow=rp, ncol=ncol(Zx2))

#パラメータの格納用配列
beta11 <- matrix(0, nrow=k, ncol=k1)
beta12 <- matrix(0, nrow=k, ncol=k2)
beta21 <- rep(0, k)
beta22 <- rep(0, k)

#ベストな初期パラメータが選択されるまで反復される
for(m in 1:rp){
  for(j in 1:k){
    p11 <- abs(colSums(X11, na.rm=TRUE)/sum(X11, na.rm=TRUE) + runif(k1, -0.1, 0.15))
    p12 <- abs(colSums(X12, na.rm=TRUE)/sum(X12, na.rm=TRUE) + runif(k2, -0.1, 0.15))
    beta11[j, ] <- p11 / sum(p11)
    beta12[j, ] <- p12 / sum(p12)
    
    beta21[j] <- runif(1, 0.2, 0.6)
    beta22[j] <- runif(1, 0.2, 0.6)
  }
  beta <- list(beta11=beta11, beta12=beta12, beta21=beta21, beta22=beta22)
  x[[rp]] <- beta 
  
  #潜在クラス割当のパラメータの初期値
  theta0 <- runif(k-1, -0.55, 0.55)
  theta1 <- runif((k-1)*cont, 0, 0.6)
  theta2 <- runif((k-1)*bin, -0.7, 0.7)
  theta3 <- runif((k-1)*(multi-1), -0.7, 0.7)
  theta4 <- runif((k-1)*(pattern-1), -1.0, 1.0)
  theta <- c(theta0, theta1, theta2, theta3, theta4)
  theta.x[rp, ] <- theta
  
  #事前確率の初期値
  logit <- matrix(Zx2 %*% theta, nrow=hh, ncol=k, byrow=T)
  r <- exp(logit)/matrix(rowSums(exp(logit)), nrow=hh, ncol=k)
  
  #観測データの対数尤度を計算
  oll <- LLobz(beta11, beta12, beta21, beta22, ns11, ns12, ns21, ns22, X11, X12, X21, X22, z.vec, r, k, n_seg)
  val <- c(val, oll$LLob)
  Z.list[[rp]] <- oll$Z
  print(m)
}
