#####ベイジアン混合多項プロビットモデル#####
library(MASS)
library(bayesm)
library(condMVNorm)
library(MCMCpack)
library(gtools)
library(MNP)
library(reshape2)
library(plyr)
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

#多変量回帰モデルの相関行列を作成
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
#set.seed(8437)
##データの設定
seg <- 3   #セグメント数
hh <- 500   #セグメントごとの購入者数
H <- seg*hh
choise <- 5   #選択可能数
st <- 5   #基準ブランド
k <- 5   #回帰係数の数
g <- rep(1:seg, rep(hh, seg))

##説明変数の発生
#通常価格の発生
PRICE <- matrix(runif(H*choise, 0.7, 1), nrow=H, ncol=choise)   

#ディスカウント率の発生
DISC <- matrix(runif(H*choise, 0, 0.3), nrow=H, ncol=choise)

#特別陳列の発生
DISP <- matrix(0, nrow=H, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(H, 1, r)
}

#特別キャンペーンの発生
CAMP <- matrix(0, nrow=H, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(H, 1, r)
}

##セグメントごとの分散共分散行列の設定
Cov <- list()
for(i in 1:seg){
  corM <- corrM(col=choise-1, lower=-0.55, upper=0.60)   #相関行列を作成
  Sigma <- covmatrix(col=choise-1, corM=corM, lower=1, upper=4)   #分散共分散行列
  Cov[[i]] <- Sigma$covariance
}


##セグメントの設定
z.seg <- rep(1:seg, rep(hh, seg))

##パラメータの設定
beta1 <- c(-9.0, -7.0, -5.0)   #価格のパラメータ
beta2 <- c(8.2, 5.5, 9.2)   #割引率のパラメータ
beta3 <- c(0.7, 2.3, 3.2)   #特別陳列のパラメータ
beta4 <- c(0.5, 1.2, 3.3)   #キャンペーンのパラメータ

#ブランド1～4の相対ベース販売力
beta0 <- matrix(0, nrow=seg, ncol=choise-1)
for(i in 1:seg){
  beta0[i, ] <- c(runif(1, -1.2, 2.8), runif(1, -1.5, 3.6), runif(1, -1.0, 4.3), runif(1, -1.5, 4.0))   
}

#回帰係数を結合
betat <- matrix(0, nrow=seg, ncol=ncol(beta0)+k-1)
for(i in 1:seg){
  betat[i, ] <- c(beta0[i, ], beta1[i], beta2[i], beta3[i], beta4[i])
}
round(betat, 2)


##相対効用を発生させ、選択されたブランドを決定
#基準ブランドとの相対説明変数
PRICE.r <- PRICE[, -5] - PRICE[, 5]
DISC.r <- DISC[, -5] - DISC[, 5]
DISP.r <- DISP[, -5] - DISP[, 5]
CAMP.r <- CAMP[, -5] - CAMP[, 5]

#相対効用の平均構造を発生
U.mean <- matrix(0, nrow=H, ncol=choise-1)
for(s in 1:seg){
  index <- subset(1:length(g), g==s)
  for(b in 1:ncol(beta0)){
    U.mean[index, b] <- beta0[s, b] + PRICE.r[g==s, b]*beta1[s] + DISC.r[g==s, b]*beta2[s] + 
                        DISP.r[g==s, b]*beta3[s] + CAMP.r[g==s, b]*beta4[s]
  }
}

#誤差構造を加えた相対効用
U <- matrix(0, nrow=H, ncol=choise-1)
for(s in 1:seg){
  index <- subset(1:length(g), g==s)
  U[index, ] <- t(apply(U.mean[g==s, ], 1, function(x) mvrnorm(1, x, Cov[[s]])))
}

#効用最大化原理に基づき選択ブランドを決定
Y <- apply(U, 1, function(x) ifelse(max(x) < 0, 5, which.max(x)))

#購買を0、1行列に変更
BUY <- matrix(0, nrow=H, ncol=choise)
for(i in 1:H){
  BUY[i, Y[i]] <- 1
}

table(Y)   #選択ブランドの集計
round(data.frame(Y, U=U), 1)   #効用と選択ブランドを比較


####マルコフ連鎖モンテカルロ法で多項プロビットモデルを推定####
####MCMC推定のための推定準備####

##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##MCMCアルゴリズムの設定
R <- 20000
keep <- 2

##回帰モデルを推定するために説明変数をベクトル形式に変更設定
#切片の設定
p <- c(1, rep(0, choise-1))
bp <- matrix(p, nrow=H*choise, ncol=choise-1, byrow=T)
BP <- subset(bp, rowSums(bp) > 0)

#説明変数の設定
PRICE.v <- as.numeric(t(PRICE.r))
DISC.v <- as.numeric(t(DISC.r))
DISP.v <- as.numeric(t(DISP.r))
CAMP.v <- as.numeric(t(CAMP.r))

X <- data.frame(BP=BP, PRICE.v, DISC.v, DISP.v, CAMP.v)   #データの結合
XM <- as.matrix(X)   #データ形式を行列に変換

#IDのインデックスを設定
index <- rep(1:H, rep(choise-1, H))   

#Xをベクトル化しておく
XV <- matrix(0, nrow=H, ncol=(choise-1)*ncol(XM))
for(i in 1:H){
  XV[i, ] <- as.numeric(XM[index==i, ])
}

##事前分布の設定
a <- rep(10, seg)   #ディクレリ事前分布

##サンプリング結果の保存用配列
Util <- array(0, dim=c(hh, choise-1, R/keep))
BETA <- matrix(0, nrow=R/keep, length(beta0)+k-1)
SIGMA <- array(0, dim=c(choise-1, choise-1, R/keep))
Z <- matrix(0, nrow=R/keep, ncol=hh)
THETA <- matrix(0, nrow=R/keep, ncol=seg)

##初期パラメータの設定
#多項プロビットモデルで初期値を設定
data1 <- list(p=choise, y=Y[g==1], X=XM[1:(length(g[g==1])*(choise-1)), ])
mcmc1 <- list(R=10000, keep=1)

out <- rmnpGibbs(Data=data1, Mcmc=mcmc1)   #多項プロビットモデルを推定
betaf <- colMeans(out$betadraw[5000:nrow(out$betadraw), ])   
sigmaf <- matrix(colMeans(out$sigmadraw[5000:nrow(out$sigmadraw), ]), nrow=choise-1, ncol=choise-1)

#セグメント別に初期値を設定
#betaの初期値
betaf.M <- matrix(betaf, nrow=seg, ncol=length(betaf), byrow=T)
betaold <- betaf.M + matrix(runif(seg*length(betaf), -0.5, 0.5), nrow=seg, ncol=length(betaf))

#sigmaの初期値
sigmaold <- list()
for(i in 1:seg){
  sigmaold[[i]] <- sigmaf
}

#thetaの初期値
theta <- c(0.3, 0.3, 0.4)


####マルコフ連鎖モンテカルロ法で混合多項プロビットモデルを推定####
##セグメント別に多項プロビットモデルを推定
##推定のための準備
#説明変数のセグメントの指示変数を作成
gx <- matrix(0, nrow=H, ncol=choise-1)
for(s in 1:(choise-1)){
  gx[, s] <- g
}
gx <- as.numeric(t(gx))

burnin <- 3000   #バーンイン期間
R <- 10000

#推定結果の格納用配列
beta.s <- matrix(0, nrow=seg, ncol=length(betaf))
sigma.s <- list() 
BETA.S <- list()
SIGMA.S <- list()
  
##多項プロビットモデルをギブスサンプリングでセグメント別に推定
for(s in 1:seg){
#データの設定
  data1 <- list(p=choise, y=Y[g==s], X=XM[gx==s, ])
  mcmc1 <- list(R=10000, keep=1)
  
  #多項プロビットモデルを推定
  out <- rmnpGibbs(Data=data1, Mcmc=mcmc1)   #多項プロビットモデルを推定
  beta.s[s, ] <- colMeans(out$betadraw[burnin:nrow(out$betadraw), ])   
  sigma.s[[s]] <- matrix(colMeans(out$sigmadraw[burnin:nrow(out$sigmadraw), ]), nrow=choise-1, ncol=choise-1)
  BETA.S[[s]] <- out$betadraw
  SIGMA.S[[s]] <- out$sigmadraw
}

####混合多項プロビットモデルの推定結果と要約####
##パラメータの識別性を確保
BETA.SI <- list()
SIGMA.SI <- list()

#共分散行列の(1, 1)要素を1に固定する制約
for(s in 1:seg){
  BETA.SI[[s]] <- BETA.S[[s]] / matrix(SIGMA.S[[s]][, 1], nrow=R, ncol=ncol(beta.s))
  SIGMA.SI[[s]] <- SIGMA.S[[s]] / matrix(SIGMA.S[[s]][, 1], nrow=R, ncol=(choise-1)^2)
}

##パラメータの要約
round(colMeans(BETA.SI[[1]][burnin:R, ]), 3)

  
SIGMA.S[[1]]

matrix(SIGMA.S[[1]])[, 1] / 

round(beta.s, 3)
round(betat, 3)
sigma.s
Cov

##潜在変数Zのサンプリング
#多項プロビットモデルの確率の計算
Pr <- list()
for(s in 1:seg){
  Pr[[s]] <- t(apply(XV, 1, function(x) abs(mnpProb(betaold[s, ] / sigmaold[[s]][1, 1], 
                                                    sigmaold[[s]] / sigmaold[[s]][1, 1], 
                                                    matrix(x, nrow=choise-1, ncol=ncol(XM)), r=25))))
}

#応答変数Yと対応するセグメントごとの確率を計算
Pr.y <- matrix(0, nrow=H, ncol=seg)
for(s in 1:seg){0
  Pr.y[, s] <- rowSums(Pr[[s]] * BUY)
}

#潜在変数Zを発生
z1 <- Pr.y * matrix(theta, nrow=H, ncol=length(theta), byrow=T)
zp <- z1 / rowSums(z1)
z <- t(apply(zp, 1, function(x) rmultinom(1, 1, x)))
Z <- apply(z, 1, which.max)

##ギブスサンプリングで多項プロビットモデルを推定
#固有値分解で強制的に正定値行列に修正する
for(i in 1:seg){
  UDU <- eigen(sigmaold[[s]])
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  sigmaold[[s]] <- vec %*% diag(D) %*% t(vec)
}

for(s in 1:seg){
  index.z <- subset(1:length(Z), Z==s)
  XZ <- XM[index %in% index.z, ]
  
  #データの設定
  Data1 <- list(p=choise, y=Y[index.z], X=XZ)
  Mcmc1 <- list(beta0=betaold[s, ], sigma0=sigmaold[[s]], R=1000, keep=1)
  
  #多項プロビットモデルをギブスサンプリング
  out=rmnpGibbs(Data=Data1,Mcmc=Mcmc1)
  betaold[s, ] <- colMeans(out$betadraw[500:1000, ])
  sigmaold[[s]] <- matrix(colMeans(out$sigmadraw[500:1000, ]), nrow=choise-1, ncol=choise-1)
}

##thetaをギブスサンプリング
theta <- table(Z) / sum(table(Z))
round(Pr.y[1:20, ], 3)
theta
betaold
round(betat, 2)


Pr <- list()
for(s in 1:seg){
  Pr[[s]] <- t(apply(XV, 1, function(x) abs(mnpProb(betat[s, ] / Cov[[s]][1, 1], 
                                                    Cov[[s]] / Cov[[s]][1, 1], 
                                                    matrix(x, nrow=choise-1, ncol=ncol(XM)), r=25))))
}

#応答変数Yと対応するセグメントごとの確率を計算
Pr.y <- matrix(0, nrow=H, ncol=seg)
for(s in 1:seg){
  Pr.y[, s] <- rowSums(Pr[[s]] * BUY)
}

cbind(Y, round(Pr.y, 2))

cbind(Y, round(Pr[[1]], 2), round(Pr[[2]], 2))

betat
betaold
