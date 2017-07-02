#####共クラスタリング(確率的ブロックモデル)#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(gtools)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

#set.seed(318)

####データの発生####
#データの設定
N <- 1000   #ユーザー数
K <- 150   #アイテム数
sg_n <- 4   #ユーザーのセグメント数
sg_k <- 3   #アイテムのセグメント数

##パラメータとセグメントを決定し観測データを発生させる
#ユーザーセグメントを発生
alpha.s1 <- rep(5, sg_n)
pi.s1 <- rdirichlet(1, alpha.s1)
z.s1 <- t(rmultinom(N, 1, pi.s1))
zno1 <- z.s1 %*% 1:sg_n

#アイテムセグメントの発生
alpha.s2 <- rep(5, sg_k)
pi.s2 <- rdirichlet(1, alpha.s2)
z.s2 <- t(rmultinom(K, 1, pi.s2))
zno2 <- z.s2 %*% 1:sg_k

#観測変数のパラメータの設定
#ユーザーセグメント×アイテムセグメントのベータ事前分布のパラメータを発生
theta.s <- rbeta(sg_n*sg_k, 1.0, 2.0)
round(theta.m <- matrix(theta.s, nrow=sg_n, ncol=sg_k), 3)

#ベルヌーイ分布から観測行列を発生させる
X <- matrix(0, nrow=N, ncol=K)
for(i in 1:N){
  print(i)
  for(j in 1:K){
  X[i, j] <- rbinom(1, 1, theta.m[zno1[i], zno2[j]])
  }
}


####周辺化ギブスサンプリングで共クラスタリングを推定####
##MCMCアルゴリズムの設定
R <- 20000
keep <- 4
sbeta <- 1.5

##ハイパーパラメータの設定
#ディクレリ分布のパラメータ
alpha1 <- rep(5, sg_n)
alpha2 <- rep(5, sg_k)

#ベータ分布のパラメータ
alpha_b <- 1.0
theta_b <- 1.5

####混合多項分布モデルで潜在変数の割当の初期値を設定####
##潜在変数の割当の初期値を設定
##観測データの対数尤度と潜在変数zを計算するための関数
LLobz <- function(theta, r, Y, ID, n, k, freq, v){
  #尤度と対数尤度を計算
  LLind <- matrix(0, nrow=n, ncol=k)
  for(i in 1:k){
    Li <- apply(cbind(Y, freq), 1, function(x) dmultinom(x[1:v], x[v+1], theta[i, ]))
    LLind[, i] <- Li 
  }
  
  #観測データの対数尤度と潜在変数zの計算
  #混合率
  R <- matrix(r, nrow=n, ncol=k, byrow=T)
  
  #個人別の潜在確率の計算
  LLd <- matrix(0, nrow=n, ncol=k)
  LLl <- log(LLind)
  for(i in 1:n){
    if(NROW(ID[ID[, 2]==i, ]==1)) { 
      LLd[i, ] <- LLl[i, ]  
      } else {
      LLd[i, ] <- apply(LLl[ID[, 2]==i, ], 2, sum)
    }
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


####EMアルゴリズムでユーザーの潜在変数の初期値を決定する####
##EMアルゴリズムの設定
#更新ステータス
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 1
iter <- 0
ID <- cbind(1:N, 1:N)

##初期値の設定
alpha <- rep(5, K)   #事前分布のパラメータ
r <- c(0.25, 0.25, 0.25, 0.25)   #混合率の初期値

#パラメータの初期値
theta.f <- matrix(0, nrow=sg_n, ncol=K)
for(i in 1:sg_n){
  minmax <- colSums(X)
  pf <- runif(K, min(minmax), max(minmax))
  theta.f[i, ] <- pf/sum(pf)
}

#対数尤度の初期化
L <- LLobz(theta=theta.f, r=r, Y=X, ID=ID, n=N, k=sg_n, freq=rowSums(X), v=K)
LL1 <- L$LLob
z <- L$z

##EMアルゴリズム
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  #Eステップの計算
  z <- L$z   #潜在変数zの出力
  zpt <- matrix(0, nrow=N, ncol=sg_n)
  for(i in 1:N){
    zpt[ID[, 2]==i, ] <- matrix(z[i, ], nrow=length(ID[ID[, 2]==i, 2]), ncol=sg_n, byrow=T)
  }
  
  #Mステップの計算と最適化
  #thetaの推定
  theta <- matrix(0, nrow=sg_n, ncol=K)
  for(j in 1:sg_n){
    #完全データの対数尤度からthetaの推定量を計算
    theta.seg <- (colSums(zpt[, j]*X)) / (sum(zpt[, j]*X))
    theta[j, ] <- as.matrix(theta.seg)
  }
  #混合率を推定
  r <- apply(z, 2, sum) / N
  
  #観測データの対数尤度を計算
  L <- LLobz(theta=theta, r=r, Y=X, ID=ID, n=N, k=sg_n, freq=rowSums(X), v=K)
  LL <- L$LLob   #観測データの対数尤度
  iter <- iter+1   
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

##所属セグメントを決定
Z1 <- as.numeric(t(apply(z, 1, function(x) rmultinom(1, 1, x))) %*% 1:sg_n)


####EMアルゴリズムでアイテムの潜在変数の初期値を決定する####
##EMアルゴリズムが正常に動くまでサンプリングを繰り返す
for(rp in 1:1000){
  #更新ステータス
  dl <- 100   #EMステップでの対数尤度の差の初期値
  tol <- 1
  iter <- 0
  ID <- cbind(1:K, 1:K)
  
  ##初期値の設定
  alpha <- rep(5, K)   #事前分布のパラメータ
  r <- c(0.25, 0.25, 0.25)   #混合率の初期値
  
  ##アイテムを基軸にするとユーザー数(変数)が多いのでユーザーをランダムサンプリングする
  N_samp <- 200
  index.u <- sample(1:nrow(X), N_samp)
  X_item <- t(X[index.u, ])
  
  #パラメータの初期値
  theta.f <- matrix(0, nrow=sg_k, ncol=N_samp)
  for(i in 1:sg_k){
    minmax <- colSums(X_item)
    pf <- runif(N_samp, min(minmax), max(minmax))
    theta.f[i, ] <- pf/sum(pf)
  }
  
  #対数尤度の初期化
  L <- LLobz(theta=theta.f, r=r, Y=X_item, ID=ID, n=K, k=sg_k, freq=rowSums(X_item), v=N_samp)
  LL1 <- L$LLob
  z <- L$z
  
  
  ##EMアルゴリズム
  LL.vec <- c()
  while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
    #Eステップの計算
    z <- L$z   #潜在変数zの出力
    zpt <- matrix(0, nrow=N_samp, ncol=sg_k)
    for(i in 1:K){
      zpt[ID[, 2]==i, ] <- matrix(z[i, ], nrow=length(ID[ID[, 2]==i, 2]), ncol=sg_k, byrow=T)
    }
    
    #Mステップの計算と最適化
    #thetaの推定
    theta <- matrix(0, nrow=sg_k, ncol=N_samp)
    for(j in 1:sg_k){
      #完全データの対数尤度からthetaの推定量を計算
      theta.seg <- (colSums(zpt[, j]*X_item)) / (sum(zpt[, j]*X_item))
      theta[j, ] <- as.matrix(theta.seg)
    }
    #混合率を推定
    r <- apply(z, 2, sum) / K
    
    #観測データの対数尤度を計算
    L <- LLobz(theta=theta, r=r, Y=X_item, ID=ID, n=K, k=sg_k, freq=rowSums(X_item), v=N_samp)
    LL <- L$LLob   #観測データの対数尤度
    iter <- iter+1   
    dl <- LL-LL1
    LL1 <- LL
    LL.vec <- c(LL.vec, LL)
    print(LL)
    if(is.nan(LL)==TRUE) break
  }
  print(rp)
  LL.vec <- LL.vec[is.nan(LL.vec)==FALSE]
  if(max(LL.vec)==LL.vec[length(LL.vec)]) break
}

##所属セグメントを決定
Z2 <- as.numeric(t(apply(z, 1, function(x) rmultinom(1, 1, x))) %*% 1:sg_k)

##統計量の初期値を計算
#ユーザーセグメントの統計量
M1k <- as.numeric(table(Z1))

Mkl <- matrix(0, nrow=sg_n, sg_k)
Xkl <- array(0, dim=c(N, K, sg_n*sg_k))
Zkl <- array(0, dim=c(N, K, sg_n*sg_k))
for(i in 1:sg_n){
  for(j in 1:sg_k){
    r <- (i-1) + j
    Zkl[, , r] <- as.matrix(ifelse(Z1==i, 1, 0), nrow=N, ncol=1) %*% ifelse(Z2==j, 1, 0)
    Xkl[, , r] <- X * Zkl[, , r]
    Mkl[i, j] <- sum(Xkl[, , r])
  }
}

####周辺化ギブスサンプリングで潜在変数zをサンプリング####
##ユーザーの潜在変数をサンプリング

