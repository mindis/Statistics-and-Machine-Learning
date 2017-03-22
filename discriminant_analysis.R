#####判別モデル#####
library(MASS)
library(plyr)
library(reshape2)
####多群の正準判別モデル####
####データの発生####
#set.seed(4235)
#多変量正規分布からの乱数発生
k <- 4   #群の数
val <- 6   #説明変数の数
n <- 400   #群ごとの学習に使うためのデータ数
nt <- 300   #群ごとのテストに使うためのデータ数

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
corM <- corrM(col=6, lower=-0.2, upper=0.2)
eigen(corM)

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- m[i] * m[j]
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

#分散共分散行列を作成(群ですべて同じ)
Sigma1 <- covmatrix(col=6, corM=corM, lower=10, upper=20)
Sigma2 <- covmatrix(col=6, corM=corM, lower=10, upper=20)
Sigma3 <- covmatrix(col=6, corM=corM, lower=10, upper=20)
Sigma4 <- covmatrix(col=6, corM=corM, lower=10, upper=20)
S <- list(Sigma1$covariance, Sigma2$covariance, Sigma3$covariance, Sigma4$covariance)

#群ごとの変数の平均を作成
mu1 <- c(rnorm(6, 12, 10))
mu2 <- c(rnorm(6, 18, 10))
mu3 <- c(rnorm(6, 6, 10))
mu4 <- c(rnorm(6, 24, 10))
mu <- list(mu1, mu2, mu3, mu4)

##多変量正規分布からの乱数を発生させる
k; n; nt
X <- matrix(0, 0, 6)
for(kk in 1:k){
    xx <- mvrnorm(n=n+nt, mu[[kk]], S[[kk]])
    X <- rbind(X, xx)
  }

#教師データのベクトルを作成
y <- rep(1:4, rep(700, 4))

#データを結合
YX <- as.data.frame(cbind(y, X))
by(YX, YX[, 1] , function(x) summary(x))   #群ごとの要約関数
by(YX[, 2:7], YX[, 1] , function(x) cor(x))   #群ごとの相関
plot(YX[, 2:7], col=YX[, 1])   #散布図

#テストデータと学習データを分離
#学習データ
YXl <- rbind(YX[1:400, ], YX[701:1100, ], YX[1401:1800, ], YX[2101:2500, ])
table(YXl[, 1])

#テストデータ
YXt <- rbind(YX[401:700, ], YX[1101:1400, ], YX[1801:2100, ], YX[2501:2800, ])
table(YXt[, 1])

####正準判別モデルで学習データから分類器を学習する####
#群内平均
mucat <- matrix(0, k, val)
for(i in 1:k){
muc <- apply(YXl[YXl[, 1]==i, 2:7], 2, mean)
mucat[i, ] <- (muc)
}

#群全体の平均
muall <- colMeans(YXl[, 2:7])
names(muall)[1:6] <- c("v1", "v2", "v3", "v4", "v5", "v6")

#群全体の分散共分散行列
covall <- var(YXl[, 2:7])
rownames(covall)[1:6] <- c("v1", "v2", "v3", "v4", "v5", "v6")
colnames(covall)[1:6] <- c("v1", "v2", "v3", "v4", "v5", "v6")

#群間の分散共分散行列
covB <- 1/k * t(mucat - muallm) %*% (mucat - muallm)
covB
mucat - muallm
##固有値問題を解いて群間の分離度を最大化する解を得る
covA <- solve(covall)
covA
covB
M <- eigen(covA %*% as.matrix(covB))   #固有値問題を解く
R <- M$values   #行列のランクは群数-1
a <- M$vectors   #固有ベクトルが判別関数の係数
a
(cont <- R/sum(R))   #寄与率
(cumcont <-cumsum(cont))   #累積寄与率

#第2固有値までを用いて、2次元空間上に射影した行列をプロット
#学習データに対する合成変量とプロット
Z <- as.matrix(YXl[, 2:7]) %*% t(a[1:3, ])
plot(as.data.frame(Z), col=YXl[, 1])

#テストデータに対する合成変量とプロット
Z <- as.matrix(YXt[, 2:7]) %*% t(a[1:3, ])
plot(as.data.frame(Z))

##判別性能を確認

