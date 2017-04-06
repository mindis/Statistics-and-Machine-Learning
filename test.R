######ランダムフォレスト######
library(MASS)
library(randomForest)
library(kernlab)
library(reshape2)
library(mlogit)
library(nnet)
library(plyr)

####データの発生####
#set.seed(321)
#多変量正規分布からの乱数発生
k <- 4   #群の数
val <- 10   #説明変数の数
n <- 700   #群ごとの学習に使うためのデータ数
nt <- 700   #群ごとのテストに使うためのデータ数

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

#群ごとの相関行列を作成(群ですべて同じ)
corM <- corrM(col=val, lower=-0.3, upper=0.4)
eigen(corM)

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
    diag(c) <- m   #対角行列を元の分散に戻す
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

#分散共分散行列を作成(群ですべて同じ)
Sigma1 <- covmatrix(col=val, corM=corM, lower=10, upper=20)
Sigma2 <- covmatrix(col=val, corM=corM, lower=10, upper=20)
Sigma3 <- covmatrix(col=val, corM=corM, lower=10, upper=20)
Sigma4 <- covmatrix(col=val, corM=corM, lower=10, upper=20)
S <- list(Sigma1$covariance, Sigma2$covariance, Sigma3$covariance, Sigma4$covariance)

#群ごとの変数の平均を作成
mu1 <- c(rnorm(val, 12, 10))
mu2 <- c(rnorm(val, 15, 10))
mu3 <- c(rnorm(val, 7, 10))
mu4 <- c(rnorm(val, 19, 10))
mu <- list(mu1, mu2, mu3, mu4)

##多変量正規分布からの乱数を発生させる
k; n; nt
X <- matrix(0, 0, val)
for(kk in 1:k){
  xx <- mvrnorm(n=n+nt, mu[[kk]], S[[kk]])
  X <- rbind(X, xx)
}

#教師データのベクトルを作成
y <- rep(1:4, rep(n+nt, 4))

#データを結合
YX <- as.data.frame(cbind(y, X))
by(YX, YX[, 1] , function(x) summary(x))   #群ごとの要約関数
by(YX[, 2:7], YX[, 1] , function(x) cor(x))   #群ごとの相関
plot(YX[, 2:7], col=YX[, 1])   #散布図

#テストデータと学習データを分離
#学習データ
YXd <- data.frame(id=1:nrow(YX), YX)   #idをつける
index <- sample(1:nrow(YXd), nrow(YXd)*0.5,  replace=FALSE)   #ランダムに標本を二分割

YXl <- YXd[index, ]   #データの半分を学習データに使う
table(YXl[, 2])

YXt <- YXd[-index, ]   #残りは検証用に使う
table(YXt[, 2])

##500サンプルをのブートストラップ標本をに抽出
index <- sample(1:nn, 500, replace=TRUE)
YXb <- YXt[index, ]

##木を成長させる
C <- sample(3:12, round(sqrt(val), 0), replace=FALSE)   #ランダムに変数を選ぶ
YXbc <- YXb[, c(1:2, C)]

##最良の変数を選ぶ
#多項ロジットを当てはめる
yy <- YXbc[, 2] 
res <- apply(YXbc[, 3:5], 2, function(x) multinom(yy ~ -1 + x))   #多項ロジットを当てはめる
op <- which.min(c(res[[1]]$AIC, res[[2]]$AIC, res[[3]]$AIC))   #最小のAICを選ぶ
Pr <- predict(res[[op]], YXbc[, op+1], type="probs")   #選択確率

yosoku <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
bootpre <- data.frame(YXbc[, 1:2], yosoku)

