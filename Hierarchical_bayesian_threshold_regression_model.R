#####Modeling the Stock and Threshold Effect of Personal Selling#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(actuar)
library(extraDistr)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
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
#データの設定
hh <- 2000  
pt <- rpois(hh, rgamma(hh, 20, 0.5))
hhpt <- sum(pt)

#IDの設定
user_id <- rep(1:hh, pt)
pt_id <- c()
for(i in 1:hh){
  pt_id <- c(pt_id, 1:pt[i])
}

##説明変数の生成
#個体内モデルの説明変数の生成
k1 <- 3; k2 <- 3; k3 <- 4
x1 <- matrix(runif(hhpt*k1, 0, 1), nrow=hhpt, ncol=k1)
x2 <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
P <- rpois(hhpt, 3.25)
C <- rpois(hhpt, 3.0)
x <- cbind(1, x1, x2, x3, P, C)   #データを結合

#階層モデルの説明変数の生成
#ユーザーの説明変数
k1 <- 3; k2 <- 4; k3 <- 4
u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #データを結合


##パラメータを生成
#パラメータ数
k1 <- ncol(x)
k2 <- ncol(u)

#階層モデルの回帰パラメータを生成
theta1 <- array(0, dim=c(k2, k1, 4))
theta2 <- array(0, dim=c(k2, k1, 2))
theta3 <- array(0, dim=c(k2, k1, 2))
for(j in 1:4){
  theta1[, 1, j] <- c(runif(1, 4.0, 5.25), runif(k2-1, -0.5, 0.5))
  theta1[, 2:(k1-2), j] <- rbind(runif(k1-3, -0.4, 0.4), matrix(runif((k1-3)*(k2-1), -0.4, 0.4), nrow=k2-1, ncol=k1-3))
  theta1[, (k1-1):k1, j] <- rbind(runif(2, -0.35, 0.35), matrix(runif(2*(k2-1), -0.35, 0.35), nrow=k2-1, ncol=2))
}

beta <- u %*% theta1[, , 1]
hist(rowSums(x * beta[user_id, ]))
