#####無限混合ガウス分布モデル#####
library(MASS)
library(mclust)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####多変量正規分布の乱数を発生させる関数を定義####
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
hh <- 5000
seg <- 5
k <- 4

#セグメント割当の設定
seg_id <- rep(1:seg, rep(hh/seg, seg))

##パラメータの設定
#平均構造の発生
mu0 <- matrix(runif(k*seg, -5, 5), nrow=seg, ncol=k)

#分散共分散行列の発生
Cov0 <- list()
for(j in 1:seg){
  Cor <- corrM(k, -0.55, 0.7, 0.1, 0.3)
  Cov0[[j]] <- covmatrix(k, Cor, 0.5, 2)$covariance
}

##多変量正規分布からデータを発生
Data <- matrix(0, nrow=hh, ncol=k)
for(j in 1:seg){
  Data[seg_id==j, ] <- mvrnorm(length(seg_id[seg_id==j]), mu0[j, ], Cov0[[j]])
}

#二変量ごとの結果を可視化
plot(Data[, 1:2], col=seg_id[seg_id %in% 1:5], pch=20, xlab="データ1の値", ylab="データ2の値", main="混合二変量正規分布のプロット")
plot(Data[, c(1, 3)], col=seg_id[seg_id %in% 1:5], pch=20, xlab="データ1の値", ylab="データ3の値", main="混合二変量正規分布のプロット")
plot(Data[, c(1, 4)], col=seg_id[seg_id %in% 1:5], pch=20, xlab="データ1の値", ylab="データ4の値", main="混合二変量正規分布のプロット")
plot(Data[, 2:3], col=seg_id[seg_id %in% 1:5], pch=20, xlab="データ2の値", ylab="データ3の値", main="混合二変量正規分布のプロット")
plot(Data[, c(2, 4)], col=seg_id[seg_id %in% 1:5], pch=20, xlab="データ2の値", ylab="データ4の値", main="混合二変量正規分布のプロット")
plot(Data[, 3:4], col=seg_id[seg_id %in% 1:5], pch=20, xlab="データ3の値", ylab="データ4の値", main="混合二変量正規分布のプロット")


####マルコフ連鎖モンテカルロ法で無限混合ガウス分布モデルを推定####
##アルゴリズムの設定
R <- 10000
keep <- 2
sbeta <- 1.5
iter <- 0

##事前分布の設定


