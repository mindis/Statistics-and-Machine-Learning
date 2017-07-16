#####Paired combinatorial logit model#####
library(MASS)
library(mlogit)
library(nnet)
library(flexmix)
library(caret)
library(reshape2)
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


####データの発生####
####データの設定####
hh <- 2000   #サンプル数
k <- 5   #選択数

####説明変数の発生####
#条件付きの説明変数の発生
X1.cont <- matrix(rnorm(2*(hh*member), 0, 1), nrow=hh, ncol=k*2)

X1.bin <- matrix(0, nrow=hh, ncol=k*2)
for(i in 1:(k*2)){
  X1.bin[, i]  <- rbinom(hh, 1, runif(1, 0.35, 0.6))
}

#多項型の説明変数の発生
X2.cont <- matrix(rnorm(hh, 0, 1), nrow=hh, ncol=1)
X2.bin <- matrix(rbinom(hh, 1, runif(1, 0.35, 0.7)), nrow=hh, ncol=1)


##説明変数をベクトル形式のデータフォーマットに変更
#IDを設定
id <- rep(1:hh, rep(k, hh))
choise <- rep(1:k, hh)
ID <- data.frame(no=1:length(id), id=id, choise=choise)

#切片の設定
p <- c(1, rep(0, k))
Pop <- matrix(p, nrow=hh*length(p), ncol=k, byrow=T)
Pop <- subset(Pop, rowSums(Pop) > 0)[, -k]


#多項型説明変数をベクトル形式に設定
X2v.cont <- matrix(0, hh*k, ncol=k)
X2v.bin <- matrix(0, hh*k, ncol=k)

for(i in 1:hh){
  index.v <- ((i-1)*k+1):((i-1)*k+k)
  v.cont <- diag(X2.cont[i, ], k)
  v.bin <- diag(X2.bin[i, ], k)
  X2v.cont[index.v, ] <- v.cont 
  X2v.bin[index.v, ] <- v.bin
}
X2v.cont <- X2v.cont[, -k]
X2v.bin <- X2v.bin[, -k]

#条件付き説明変数をベクトル形式に設定
X1v.cont <- matrix(0, nrow=hh*k, ncol=2)
X1v.bin <- matrix(0, nrow=hh*k, ncol=2)

for(i in 1:2){
  index.r <- ((i-1)*k+1):((i-1)*k+k)
  X1v.cont[, i] <- as.numeric(t(X1.cont[, index.r]))
  X1v.bin[, i] <- as.numeric(t(X1.bin[, index.r]))
}

##データを結合
X <- data.frame(pop=Pop, c1=X1v.cont, b1=X1v.bin, c2=X2v.cont, b2=X2v.bin)
round(XM <- as.matrix(X), 3)


####PCLモデルに基づき応答変数を発生####
##パラメータの設定
#回帰パラメータの設定
b0 <- runif(k-1, -1.0, 2.7)
b1 <- runif(2, 0, 1.3)
b2 <- runif(2, -1.0, 1.4)
b3 <- runif(k-1, 0, 1.5)
b4 <- runif(k-1, -1.3, 1.3)
b <- c(b0, b1, b2, b3, b4)
beta.t <- b


#類似度パラメータの設定
Cov <- corrM(col=k, lower=0, upper=0.7)   #類似度パラメータを発生

#効用関数の設定
logit <- matrix(XM %*% b, nrow=hh, ncol=k, byrow=T)

#PCLモデルに基づく確率の計算
#PCLモデルの分母部分
(1-Cov[1, 2]) * (exp(logit[1, 1]/(1-Cov[1, 2])) + exp(logit[1, 2]/(1-Cov[1, 2])))^(1-Cov[1, 2])
