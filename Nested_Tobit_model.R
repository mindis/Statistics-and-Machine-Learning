#####入れ子型トービットモデル#####
library(MASS)
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(lattice)

####多変量正規乱数を発生させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #固有値分解で強制的に正定値行列に修正する
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}


####データの発生####
##データの設定
N <- 2000   #サンプル数
k <- 12    #説明変数数
k.cont <- 5    #連続変数数
k.bin <- 3   #二値変数数
k.multi <- 4    #多値変数数

##IDの設定
ID <- 1:N

##説明変数の発生
#連続変数
X.cont <- matrix(rnorm(N*k.cont), nrow=N, ncol=k.cont)

#二値変数
X.bin <- matrix(0, nrow=N, ncol=k.bin)
for(i in 1:k.bin){
  p.bin <- runif(1, 0.3, 0.6)
  X.bin[, i] <- rbinom(N, 1, p.bin)
}

#多値変数
p.multi <- runif(k.multi)
X.multi <- t(rmultinom(N, 1, p.multi))
X.multi <- X.multi[, -which.min(colSums(X.multi))]

##パラメータの設定
#プロビットモデルのパラメータの設定
