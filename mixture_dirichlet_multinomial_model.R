#####多項-ディクレリ混合分布モデル#####
library(MASS)
library(vcd)
library(gtools)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
##モデルの設定
k <- 4   #セグメント数
hh <- 500   #セグメントのサンプル数
n <- hh * k   #総サンプル数
pt <- rpois(n, 10)   #サンプルごとのセッション数
pt <- ifelse(pt < 2, 2, pt)
N <- sum(pt)   #総サンプル数

th <- 20   #セグメント別のパラメータ数
freq <- rpois(N, 25)   #個人ごとのの頻度

##IDの設定
id <- rep(1:n, pt)   #idの設定
t <- c()   #セッション数の設定
for(i in 1:n){
  t <- c(t, 1:pt[i])
}  
no <- 1:N   #整理番号を設定
ID <- data.frame(no, id, t)   #データを結合

#セグメントの設定
seg.z <- rep(1:4, rep(hh, k))

##セグメントごとに確率ベクトルを定義
#確率ベクトルの係数の設定
x <- matrix(rnorm(th*k, 0, 1), nrow=k, ncol=th, byrow=T)
a.vec <- matrix(0, nrow=k, ncol=th)
a <- c(0.2, 0.5, 0.9, 1.2)
for(i in 1:k){
  a.vec[i, ] <- rnorm(th, 0, a[i]) 
}

#確率ベクトルの計算
Pr <- matrix(0, nrow=k, ncol=th)
for(i in 1:k){
  Pr[i, ] <- exp(x[i, ]*a.vec[i, ]) / sum(exp(x[i, ]*a.vec[i, ]))
}
#発生させた確率を確認
round(Pr, 3); apply(Pr, 1, summary)


##発生させた確率にもとづき頻度データを生成
#頻度データを発生
Y <- matrix(0, nrow=N, ncol=th)
for(i in 1:N){
  Y[i, ] <- t(rmultinom(1, freq[i], Pr[seg.z[ID[i, 2]], ]))
}

#出現率を計算
YP <- matrix(0, nrow=k, ncol=th)  
for(i in 1:k){
  index.z <- subset(1:length(seg.z), seg.z==i)
  YP[i, ] <- colSums(Y[ID[, 2] %in% index.z, ]) / sum(Y[ID[, 2] %in% index.z, ])
}
round(YP, 3)
round(Pr, 3)   #セグメント別の確率
round(colSums(Y)/sum(Y), 3)   #全体での出現率


####EMアルゴリズムで混合多項-ディクレリモデルを推定####
##観測データの対数尤度と潜在変数zを計算するための関数
LLobz <- function(theta, r, Y, ID, n, k, freq, v){
  #尤度と対数尤度を計算
  LLind <- matrix(0, nrow=N, ncol=k)
  for(i in 1:k){
    Li <- apply(cbind(Y, freq), 1, function(x) dmultinom(x[1:v], x[v+1], theta[i, ]))
    LLind[, i] <- Li 
  }
  
  #観測データの対数尤度と潜在変数zの計算
  #混合率
  R <- matrix(r, nrow=n, ncol=k)
  
  #個人別の潜在確率の計算
  LLd <- matrix(0, nrow=n, ncol=k)
  LLl <- log(LLind)
  for(i in 1:n){
    LLd[i, ] <- apply(LLl[ID[, 2]==i, ], 2, sum)
  }
  
  #桁落ちを防ぐために対数尤度が-744以下の場合は対数尤度を嵩上げする
  LL.min <- apply(LLd, 1, min)
  index.loss <- subset(1:nrow(LLd), (LL.min + 743) < 0)
  lplus <- -matrix((LL.min[index.loss] + 743), nrow=length(index.loss), ncol=k)
  LLd[index.loss, ] <- LLd[index.loss, ] + lplus
  
  #潜在確率zの計算
  LLho <- R * exp(LLd)
  z <- LLho / matrix(rowSums(LLho), nrow=n, ncol=k)
  
  #観測データの対数尤度を計算
  LLosum <- sum(log(apply(matrix(r, nrow=n, ncol=k, byrow=T) * exp(LLd), 1, sum)))
  rval <- list(LLobz=LLosum, z=z, LL=LLd)
  return(rval)
}

##EMアルゴリズムの設定
#更新ステータス
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 1
iter <- 0

#初期値の設定
alpha <- rep(2, th)   #事前分布のパラメータ
r <- c(0.4, 0.3, 0.2, 0.2)   #混合率の初期値

#パラメータの初期値
theta.f <- matrix(0, nrow=k, ncol=th)
for(i in 1:k){
  minmax <- colSums(Y)
  pf <- runif(th, min(minmax), max(minmax))
  theta.f[i, ] <- pf/sum(pf)
}

#対数尤度の初期化
L <- LLobz(theta=theta.f, r=r, Y=Y, ID=ID, n=n, k=k, freq=freq, v=th)
LL1 <- L$LLob
z <- L$z

##EMアルゴリズム
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  #Eステップの計算
  z <- L$z   #潜在変数zの出力
  zpt <- matrix(0, nrow=N, ncol=k)
  for(i in 1:n){
    zpt[ID[, 2]==i, ] <- matrix(z[i, ], nrow=length(ID[ID[, 2]==i, 2]), ncol=k, byrow=T)
  }
  
  #Mステップの計算と最適化
  #thetaの推定
  theta <- matrix(0, nrow=k, ncol=th)
  for(j in 1:k){
    #完全データの対数尤度からthetaの推定量を計算
    theta.seg <- (colSums(zpt[, j]*Y) + r[j]*(alpha-1)) / (sum(zpt[, j]*Y) + r[j]*sum(alpha-1))
    theta[j, ] <- as.matrix(theta.seg)
  }
  #混合率を推定
  r <- apply(z, 2, sum) / n
  
  #観測データの対数尤度を計算
  L <- LLobz(theta=theta, r=r, Y=Y, ID=ID, n=n, k=k, freq=freq, v=th)
  LL <- L$LLob   #観測データの対数尤度
  iter <- iter+1   
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

####推定結果と要約
##推定されたパラメータ
round(theta, 3)   #推定されたパラメータ
round(Pr, 3)   #真のパラメータ
round(r, 3)   #混合率
round(data.frame(seg=apply(z, 1, which.max), z=z), 3)   #個人別のセグメントへの所属確率と所属セグメント

#推定されたパラメータのグラフ
max.theta <- max(theta)
par(mfrow=c(2,2)) 
barplot(theta[1, ], ylim=c(0, max.theta))
barplot(theta[2, ], ylim=c(0, max.theta))
barplot(theta[3, ], ylim=c(0, max.theta))
barplot(theta[4, ], ylim=c(0, max.theta))
par(mfrow=c(1, 1))

##適合度
round(LL, 3)   #最大化された対数尤度
round(AIC <- -2*LL + 2*(length(th*k)+k), 3)   #AIC
