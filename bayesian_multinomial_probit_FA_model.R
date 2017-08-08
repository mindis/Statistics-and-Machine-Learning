#####多項プロビットモデル#####
library(MASS)
library(bayesm)
library(pyhch)
library(condMVNorm)
library(MCMCpack)
library(glmm)
library(lme4)
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
hh <- 2500   #プレイヤー数
choise <- 10   #選択可能数
st <- choise   #基準ブランド
k <- 5   #回帰係数の数

##IDの設定
id <- rep(1:hh, rep(choise-1, hh))
c <- rep(1:(choise-1), hh)
ID <- data.frame(no=1:(hh*(choise-1)), id=id, c=c)
id_r <- matrix(1:(hh*(choise-1)), nrow=hh, ncol=choise-1, byrow=T)


##説明変数の発生
#通常価格の発生
PRICE <- matrix(runif(hh*choise, 0.7, 1), nrow=hh, ncol=choise)   

#ディスカウント率の発生
DISC <- matrix(runif(hh*choise, 0, 0.3), nrow=hh, ncol=choise)

#特別陳列の発生
DISP <- matrix(0, nrow=hh, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hh, 1, r)
}

#特別キャンペーンの発生
CAMP <- matrix(0, nrow=hh, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hh, 1, r)
}

##分散共分散行列の設定
corM <- corrM(col=choise-1, lower=-0.55, upper=0.75)   #相関行列を作成
Sigma <- covmatrix(col=choise-1, corM=corM, lower=1, upper=1)   #分散共分散行列
Cov <- Sigma$covariance

##パラメータの設定
beta1 <- -6.3   #価格のパラメータ
beta2 <- 7.2   #割引率のパラメータ
beta3 <- 2.0   #特別陳列のパラメータ
beta4 <- 1.8   #キャンペーンのパラメータ
beta0 <- runif(choise-1, -0.4, 2.1)   #ブランド1～4の相対ベース販売力
betat <- c(beta0, beta1, beta2, beta3, beta4)

#基準ブランドとの相対説明変数
PRICE.r <- PRICE[, -5] - PRICE[, 5]
DISC.r <- DISC[, -5] - DISC[, 5]
DISP.r <- DISP[, -5] - DISP[, 5]
CAMP.r <- CAMP[, -5] - CAMP[, 5]

##回帰モデルを推定するために説明変数をベクトル形式に変更設定
#切片の設定
p <- c(1, rep(0, choise-1))
bp <- matrix(p, nrow=hh*choise, ncol=choise-1, byrow=T)
BP <- subset(bp, rowSums(bp) > 0)


#説明変数の設定
PRICE.v <- as.numeric(t(PRICE.r))
DISC.v <- as.numeric(t(DISC.r))
DISP.v <- as.numeric(t(DISP.r))
CAMP.v <- as.numeric(t(CAMP.r))

X <- data.frame(BP=BP, PRICE.v, DISC.v, DISP.v, CAMP.v)   #データの結合
XM <- as.matrix(X)

##相対効用を発生させ、選択されたブランドを決定
U.mean <- matrix(XM %*% betat, nrow=hh, ncol=choise-1, byrow=T)   #相対効用の平均構造
U <- U.mean + mvrnorm(hh, rep(0, choise-1), Cov)   #誤差構造を加えた効用

#効用最大化原理に基づき選択ブランドを決定
Y <- apply(U, 1, function(x) ifelse(max(x) < 0, choise, which.max(x)))


#購買を0、1行列に変更
BUY <- matrix(0, hh, choise)
for(i in 1:hh){
  BUY[i, Y[i]] <- 1
}

table(Y)   #選択ブランドの集計
round(cbind(Y, U, U.mean), 2)   #効用と選択ブランドを比較


####マルコフ連鎖モンテカルロ法で多項プロビットモデルを推定####
####MCMC推定のための推定準備####

##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##多変量正規分布の条件付き期待値と分散を計算する関数
cdMVN <- function(mu, Cov, dependent, U){
  
  #分散共分散行列のブロック行列を定義
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent]
  Cov21 <- Cov[-dependent, dependent]
  Cov22 <- Cov[-dependent, -dependent]
  
  #条件付き分散と条件付き平均を計算
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mu[, dependent] + t(CDinv %*% t(U[, -dependent] - mu[, -dependent]))   #条件付き平均を計算
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #条件付き分散を計算
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
factors <- 3
llike <- array(0, dim=c(R/keep))   #対数尤度の保存用


##説明変数を多次元配列化
X.array <- array(0, dim=c(choise-1, ncol(X), hh))
for(i in 1:hh){
  X.array[, , i] <- XM[ID[, 2]==i, ]
}
YX.array <- array(0, dim=c(choise-1, ncol(X)+1, hh))


#推定プロセスの格納配列
UM <- matrix(0, nrow=hh, ncol=choise-1)
util.M <- matrix(0, nrow=hh, ncol=choise-1)   

##事前分布の設定
nu <- choise   #逆ウィシャート分布の自由度
V <- solve((1/10)*diag(choise-1))    #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(X))  #回帰係数の平均の事前分布
Adelta <- solve(100 * diag(rep(1, ncol(X))))   #回帰係数の事前分布の分散

inv.facov <- solve(diag(100, factors))
alpha_d <- 1
beta_d <- 100


##サンプリング結果の保存用配列
Util <- array(0, dim=c(hh, choise-1, R/keep))
BETA <- matrix(0, nrow=R/keep, length(beta0)+k-1)
SIGMA <- matrix(0, nrow=R/keep, ncol=(choise-1)^2)
FA.A <- matrix(0, nrow=R/keep, ncol=(choise-1)*factors)
FA.D <- matrix(0, nrow=R/keep, ncol=choise-1)
FA.F <- array(0, dim=c(hh, factors, R/keep))


##初期値の設定
#回帰係数の初期値
oldbeta <- c(runif(choise-1, 0, 3), -3.0, 3.0, runif(2, 0, 2))   

#分散共分散行列の初期値
corM.f <- corrM(col=choise-1, lower=0, upper=0)   #相関行列を作成
Sigma.f <- covmatrix(col=choise-1, corM=corM.f, lower=1, upper=1)   #分散共分散行列
oldcov <- Sigma.f$covariance

#効用の平均構造の初期値
old.utilm <- matrix(XM %*% oldbeta, nrow=hh, ncol=choise-1, byrow=T)

#効用の初期値
old.util <- old.utilm + mvrnorm(nrow(old.utilm), rep(0, choise-1), oldcov)

#因子負荷量と独自因子の初期値
A <- matrix(runif((choise-1)*factors, -1, 1), nrow=choise-1, ncol=factors)
D <- diag(runif(choise-1, 0, 0.5))


####マルコフ連鎖モンテカルロ法で多項プロビットモデルを推定####
for(rp in 1:R){
  
  ##選択結果と整合的な潜在効用を発生させる
  #条件付き期待値と条件付き分散を計算
  S <- rep(0, choise-1)
  
  for(j in 1:(choise-1)){
    MVR <- cdMVN(mu=old.utilm, Cov=oldcov, dependent=j, U=old.util)   #条件付き分布を計算
    UM[, j] <- MVR$CDmu   #条件付き期待値を取り出す
    S[j] <- sqrt(MVR$CDvar)    #条件付き分散を取り出す
    
    #潜在変数を発生させる
    #切断領域の設定
    max.u <- apply(cbind(old.util[, -j], 0), 1, max)
    max.u <- ifelse(Y==choise, 0, max.u)
    
    #切断正規分布より潜在変数を発生
    old.util[, j] <- ifelse(Y==j, rtnorm(mu=UM[, j], sigma=S[j], a=max.u, b=100), 
                            rtnorm(mu=UM[, j], sigma=S[j], a=-100, b=max.u))
    old.util[, j] <- ifelse(is.infinite(old.util[, j]), ifelse(Y==j, max.u + runif(1), max.u - runif(1)), old.util[, j])
  }
  util.v <- as.numeric(t(old.util))
  
  ##betaの分布のパラメータの計算とmcmcサンプリング
  #z.vecとX.vecを結合して多次元配列に変更
  YX.bind <- cbind(util.v, XM)
  for(i in 1:hh){
    YX.array[, , i] <- YX.bind[id_r[i, ], ]
  }
  
  ##回帰モデルのギブスサンプリングでbetaとsigmaを推定
  #betaのギブスサンプリング
  invcov <- solve(oldcov)
  xvx.vec <- rowSums(apply(X.array, 3, function(x) t(x) %*% invcov %*% x))
  XVX <- matrix(xvx.vec, nrow=ncol(X), ncol=ncol(X), byrow=T)
  XVY <- rowSums(apply(YX.array, 3, function(x) t(x[, -1]) %*% invcov %*% x[, 1]))
  
  #betaの分布の分散共分散行列のパラメータ
  inv_XVX <- solve(XVX + Adelta)
  
  #betaの分布の平均パラメータ
  B <- inv_XVX %*% (XVY + Adelta %*% Deltabar)   #betaの平均
  b1 <- as.numeric(B)
  
  #多変量正規分布から回帰係数をサンプリング
  oldbeta <- mvrnorm(1, b1, inv_XVX)
  
  ##Covの分布のパラメータの計算とmcmcサンプリング
  #逆ウィシャート分布のパラメータを計算
  R.error <- matrix(util.v - XM %*% oldbeta, nrow=hh, ncol=choise-1, byrow=T)
  IW.R <- V + matrix(rowSums(apply(R.error, 1, function(x) x %*% t(x))), nrow=choise-1, ncol=choise-1)
  
  #逆ウィシャート分布の自由度を計算
  Sn <- nu + hh
  
  #逆ウィシャート分布からCovをサンプリング
  Cov_hat <- rwishart(Sn, solve(IW.R))$IW
  oldcov <- cov2cor(Cov_hat)
  
  ##潜在効用とパラメータを更新
  #潜在効用と潜在効用の平均を更新
  old.utilm <- matrix(XM %*% oldbeta, nrow=hh, ncol=choise-1, byrow=T)
  Z <- old.util - old.utilm 
  Z <- scale(Z)
  
  ##潜在効用の誤差項から因子分析モデルを推定
  #多変量正規分布から潜在変数f(共通因子)をサンプリング
  ADA <- t(A) %*% solve(A %*% t(A) + D)
  F_mean <- Z %*% t(ADA)   #共通因子の平均
  F_var <- diag(factors) - ADA %*% A    #共通因子の分散共分散行列
  Fi <- t(apply(F_mean, 1, function(x) mvrnorm(1, x, F_var)))   #多変量正規分布から共通因子をサンプリング 
  
  
  #ガンマ分布から独自因子dをサンプリング
  Z.error <- Z - Fi %*% t(A)
  Zv.R <- matrix(rowSums(apply(Z.error, 1, function(x) x %*% t(x))), nrow=choise-1, ncol=choise-1)
  
  gamma_alpha <- (hh + alpha_d)/2   #alphaを計算
  gamma_beta <- (diag(Zv.R) + beta_d)/2   #betaを計算
  D <- diag(rgamma(length(gamma_beta), gamma_alpha, gamma_beta))   #ガンマ分布から独自因子をサンプリング
  
  #多変量正規分布から因子負荷量Aをサンプリング
  FF <- t(Fi) %*% Fi
  FZ <- t(Fi) %*% Z
  d_sigma <- list()
  
  for(i in 1:(choise-1)){
    d_sigma[[i]]  <- 1/diag(D)[i] * inv.facov
    A_mu <- solve(d_sigma[[i]] + FF) %*% FZ[, i]
    A_cov <- solve(inv.facov + diag(D)[i]*FF) 
    A[i, ] <- mvrnorm(1, A_mu, A_cov)
  }
  A[1, 1] <- 0
  A[2, 2] <- 0
  A[3, 3] <- 0
  
  ##サンプリング結果を保存
  if(rp%%keep==0){
    print(rp)

    mkeep <- rp/keep
    Util[, , mkeep] <- old.util
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    FA.A[mkeep, ] <- as.numeric(A)
    FA.D[mkeep, ] <- diag(D)
    FA.F[, , mkeep] <- Fi
    print(round(cbind(oldcov, Cov), 2))
    print(round(rbind(oldbeta, betat), 2))
    print(round(t(A), 2))
  }
}


####関数で推定####
Data1 <- list(p=choise, y=Y, X=XM)
Mcmc1 <- list(R=10000, keep=2)

#多項プロビットモデルを推定
out <- rmnpGibbs(Data=Data1,Mcmc=Mcmc1)
BETA.out <- out$betadraw
SIGMA.out <- out$sigmadraw

####推定結果の要約と適合度の確認####
burnin <- 10000/keep   #バーンイン期間

##サンプリング結果を可視化
#回帰係数のプロット
matplot(BETA[, 1:4], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 5:8], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")

#分散供分散行列の可視化
matplot(SIGMA[, 1:4], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 5:9], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 10:13], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 14:18], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 19:22], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 23:27], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 28:31], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 32:36], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 37:40], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 41:45], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 46:49], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 50:54], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 55:58], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 59:63], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 64:67], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 68:72], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 73:76], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 77:81], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")

#因子負荷量の可視化
matplot(FA.A[, 1:4], type="l", main="因子負荷量のサンプリング結果", ylab="パラメータ推定値")
matplot(FA.A[, 5:9], type="l", main="因子負荷量のサンプリング結果", ylab="パラメータ推定値")
matplot(FA.A[, 10:13], type="l", main="因子負荷量のサンプリング結果", ylab="パラメータ推定値")
matplot(FA.A[, 14:18], type="l", main="因子負荷量のサンプリング結果", ylab="パラメータ推定値")
matplot(FA.A[, 19:22], type="l", main="因子負荷量のサンプリング結果", ylab="パラメータ推定値")
matplot(FA.A[, 23:27], type="l", main="因子負荷量のサンプリング結果", ylab="パラメータ推定値")


##推定値の事後平均の比較
#betaの要約統計量
round(colMeans(BETA.out[burnin:nrow(BETA.out), ] / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(関数推定)の事後平均
round(colMeans(BETA[burnin:nrow(BETA), ]), 3)   #betaの事後平均
round(betat, 3)   #真の値
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(BETA[burnin:nrow(BETA), ], 2, sd), 2)   #事後標準偏差

#sigmaの要約統計量
round(colMeans(SIGMA.out[burnin:nrow(SIGMA.out), ]  / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(関数推定)の事後平均
round(colMeans(SIGMA[burnin:nrow(SIGMA), ]), 3)   #betaの事後平均
round(as.numeric(Cov), 3)   #真の値
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, sd), 2) #事後標準偏差
