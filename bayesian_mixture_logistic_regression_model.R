#####ベイジアン混合ロジスティック回帰モデル####
library(MASS)
library(flexmix)
library(MCMCpack)
library(bayesm)
library(reshape2)
library(dplyr)
library(caret)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(452489)
k <- 3   #セグメント数
col <- 10   #変数数
n <- 500   #セグメントごとのサンプル数
N <- k*n   #全サンプル数
pt <- 5   #1人あたりの購買機会数
hh <- N/pt
ID <- rep(1:hh, rep(pt, hh))
seg.z <- rep(1:k, rep(n, k))


##説明変数の設定
#連続変数の発生
X1 <- matrix(runif(N*(col-4), 0, 1), N, (col-4))

#二値変数の発生
X2 <- matrix(0, N, (col-ncol(X1)))
for(i in 1:(col-ncol(X1))){
  bin <- rbinom(N, 1, runif(1, 0.2, 0.7))
  X2[, i] <- bin
}
X <- cbind(X1, X2)

##回帰係数の発生
for(t in 1:1000){
  lower <- c(-1.4, -1.0, -0.5)
  upper <- c(1.4, 1.0, 0.6)
  b1 <- matrix(0, nrow=k, ncol=col)
  b0 <- c()
  
  for(i in 1:k){
    b0[i] <- runif(1, -1.5, 1.5)
    b1[i, ] <- runif(col, lower[i], upper[i])
  }
  
  betat <- data.frame(b0, b=b1)   #真の回帰係数
  round(betat, 3)
  
  ##ロジスティック回帰のリンク関数と確率の計算
  Pr <- matrix(0, nrow=N, ncol=k)
  
  #確率を計算
  for(i in 1:k){
    logit <- b0[i] + X %*% b1[i, ]
    Pr[, i] <- exp(logit)/(1+exp(logit))
  }
  cbind(Pr[, 3], rbinom(N, 1, Pr[, 3]))
  
  ##ベルヌーイ乱数で二値データを発生
  Y <- c()
  for(i in 1:k){
    y <- rbinom(length(seg.z[seg.z==i]), 1, Pr[seg.z==i, i])
    Y <- c(Y, y)
  }
  
  YX <- data.frame(seg=seg.z, Y, X)   #すべてのデータを結合
  if(max(table(YX$seg, YX$Y)) < 400) break   #セグメント別でのクロス集計
}

##データの結合と集計
YX <- data.frame(seg=seg.z, Y, X)   #すべてのデータを結合
table(Y)   #全体でのyの単純集計
table(YX$seg, YX$Y)   #セグメント別でのクロス集計
round(table(YX$seg, YX$Y) / rowSums(table(YX$seg, YX$Y)), 3)   #セグメント別での比率クロス集計


#確率分布をプロット
#セグメントでの分布
hist(Pr[seg.z==1, 1], col="#0000ff40", xlim=c(0, 1.0), border = "#0000ff", xlab="rate", main="確率の分布")
hist(Pr[seg.z==2, 2], xlim=c(0, 1), col="#ff00ff40", border = "#ff00ff", xlab="rate", main="確率の分布")
hist(Pr[seg.z==3, 3], xlim=c(0, 1), col="#a5f0ff40", border = "#ff00ff", xlab="rate", main="確率の分布")

#全体での分布
PP <- c()
for(i in 1:k){
  PP <- c(PP, Pr[seg.z==i, i])
}
hist(PP, breaks=20, xlim=c(0, 1), col="#0000ff40", border = "#ff00ff",
     xlab="rate", main="確率の分布")


####ベイジアン混合ロジスティックモデルを推定####
####マルコフ連鎖モンテカルロ法推定のための準備####

##尤度と対数尤度を計算する関数
logl <- function(b, X, Y){
  #パラメータの設定
  alpha <- b[1]
  beta <- b[2:(col+1)]
  
  #尤度を定義して合計する
  logit <- alpha + as.matrix(X) %*% beta 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  val <- list(LL=LL, LLS=LLS, p=p)
  return(val)
}

##ロジスティック回帰モデルの対数尤度を定義
loglike <- function(b, X, Y){
  #パラメータの設定
  alpha <- b[1]
  beta <- b[2:(col+1)]
  
  #尤度を定義して合計する
  logit <- alpha + as.matrix(X) %*% beta 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}



##対数尤度を最大化する
b0 <- c(rep(0, col+1))   #初期パラメータの設定
res <- optim(b0, loglike, gr=NULL, X=X, Y=Y, method="BFGS", hessian=TRUE, control=list(fnscale=-1))

beta0 <- res$par[1]  
beta <- res$par[-1]
H <- res$hessian
invH <- solve(-H)


##MCMCアルゴリズムの設定
R <- 20000   #サンプリング回数
keep <- 2   #2回に1回の割合でサンプリング結果を利用

##事前分布の設定
#回帰係数の事前分布の設定
betas <- rep(0, col+1)
B0 <- 0.01*diag(col+1)

sbeta <- 0.2
rw <- t(chol(sbeta*invH))   #ランダムウォークの分散
rootBi <- t(chol(B0))   #事前分布の精度

#ディクレリ事前分布の設定
a <- rep(10, k)

##パラメータの保存用配列
#推定結果の格納用配列
BETA <- matrix(0, nrow=R/keep, ncol=(col+1)*k) 
THETA <- matrix(0, nrow=R/keep, ncol=k)
ZP <- array(0, dim=c(N, k, R/keep)) 
Z <- matrix(0, nrow=R/keep, ncol=N)

##初期値の設定
oldbeta <- matrix(res$par, nrow=k, ncol=col+1, byrow=T) + matrix(runif(k*(col+1)), nrow=k, ncol=col+1)
theta <- c(0.25, 0.35, 0.4)


####マルコフ連鎖モンテカルロ法で混合ロジスティックモデルを推定####
##潜在変数Zの発生
#尤度の計算
L <- matrix(0, nrow=N, ncol=k)
for(i in 1:k){
  ll <- logl(oldbeta[i, ], X, y)
  L[, i] <- exp(ll$LLS)
}

#個人ごとの尤度計算
LLh <- matrix(0, nrow=hh, ncol=k)
for(i in 1:hh){
  LLh[i, ] <- apply(L[ID==i, ], 2, prod)
}

#潜在確率の計算
R <- matrix(theta, nrow=hh, ncol=k, byrow=T)   #混合率
LLr <- R * LLh
z1 <- LLr / matrix(rowSums(LLr), nrow=hh, ncol=k)

#潜在変数zの発生
z <- matrix(0, nrow=N, ncol=k)
for(i in 1:hh){
  r <- ((i-1)*pt+1):((i-1)*pt+pt)
  zz <- t(rmultinom(pt, 1, z1[i, ]))
  z[r, ] <- zz
}
zn <- z %*% 1:3

##メトロポリスヘイスティングアルゴリズムでbetaを更新

