#####二変量プロビットモデルによるサンプルセレクションモデル#####
library(MASS)
library(sampleSelection)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  X.Sigma <- eigen(Sigma)
  Lambda <- diag(X.Sigma$values)
  
  #新しい相関行列の定義と対角成分を1にする
  Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda)
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
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


####説明変数の発生####
N <- 5000   #サンプル数

##利用歴があるかどうかを決定する説明変数
cont1 <- 3; bin1 <- 4; multi1 <- 4
Z.cont <- matrix(rnorm(N*cont1), nrow=N, ncol=cont1)
Z.bin <- matrix(0, nrow=N, ncol=bin1)
Z.multi <- matrix(0, nrow=N, ncol=multi1)

#二値説明変数を設定
for(i in 1:bin1){
  p <- runif(1, 0.3, 0.7)
  Z.bin[, i] <- rbinom(N, 1, p)
}

#多値説明変数を設定
p <- runif(multi1)
Z.multi <- t(rmultinom(N, 1, p))
Z.multi <- Z.multi[, -which.min(colSums(X.multi))] #冗長な変数は削除

#データを結合
Z <- cbind(1, Z.cont, Z.bin, Z.multi)


##利用したならどれくらい利用したのかを決定する説明変数
#連続変数の発生
cont2 <- 3
X.cont <- matrix(rnorm(N*cont2, 0, 1), nrow=N, ncol=cont2)

#二値変数の発生
bin2 <- 5
X.bin <- matrix(0, nrow=N, ncol=bin2)
for(i in 1:bin2){
  X.bin[, i] <- rbinom(N, 1, runif(1, 0.35, 0.65))
}

#データの結合
X <- cbind(1, X.cont, X.bin)


####二変量正規分布から応答変数を発生####
#分散共分散行列を設定
sigma0 <- 1.5
rho0 <- 0.6
Cov0 <- matrix(c(1, rho0*sqrt(sigma0), rho0*sqrt(sigma0), sigma0), nrow=2, ncol=2)


for(i in 1:1000){
  print(i)
  #パラメータを設定
  alpha0 <- c(runif(1, -0.8, -0.5), runif(cont1, 0, 0.8), runif(bin1+multi1-1, -0.5, 1.0))
  beta0 <- c(runif(1, 0.3, 0.5), runif(cont2, 0, 0.5), runif(bin2, -0.3, 0.5))
  
  #応答変数の平均構造
  y1 <- Z %*% alpha0
  y2 <- X %*% beta0
  
  #多変量正規分布から応答変数を発生
  Y <- t(apply(cbind(y1, y2), 1, function(x) mvrnorm(1, x, Cov0)))
  y2 <- exp(Y[, 2])
  y1 <- ifelse(Y[, 1] > 0, 1, 0)
  
  if(max(y2) < 150) break
}
summary(cbind(y1, y2))
cor(Y)
