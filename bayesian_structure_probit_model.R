#####構造プロビットモデル#####
library(MASS)
library(bayesm)
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
hh <- 1000   #プレイヤー数
member <- 10   #選択可能数
k <- 5   #回帰係数の数

##IDの設定
id <- rep(1:hh, rep(member-1, hh))
c <- rep(1:(member-1), hh)
ID <- data.frame(no=1:(hh*(member-1)), id=id, c=c)
id_r <- matrix(1:(hh*(member-1)), nrow=hh, ncol=member-1, byrow=T)


##説明変数の発生
cont <- 2
X1.cont <- matrix(rnorm(hh*member*2, 0, 1), nrow=hh, ncol=2*member)

bin <- 2
X1.bin <- matrix(0, nrow=hh, ncol=2*member)
for(i in 1:(2*member)){
  X1.bin[, i]  <- rbinom(hh, 1, runif(1, 0.35, 0.6))
}



##説明変数をベクトル形式のフォーマットに変更
#IDを設定
id <- rep(1:hh, rep(member-1, hh))
m <- rep(1:(member-1), hh)
ID <- data.frame(no=1:length(id), id=id, m=m)

#切片の設定
p <- c(1, rep(0, member-1))
Pop <- matrix(p, nrow=hh*member, ncol=member-1, byrow=T)
POP <- subset(Pop, rowSums(Pop) > 0)

#構造変数の設定
c <- 3 
g1 <- rep(c(rep(1, (member-1)/c), rep(0, 2*(member-1)/c)), hh)
g2 <- rep(c(rep(0, (member-1)/c), rep(1, (member-1)/c), rep(0, (member-1)/c)), hh)
g3 <- rep(c(rep(0, 2*(member-1)/c), rep(1, (member-1)/c)), hh)
G <- cbind(g1, g2)


#説明変数の設定
#相対効用に変更
X1r.cont1 <- X1.cont[, 1:(member-1)] - X1.cont[, member]
X1r.cont2 <- X1.cont[, (member+1):(2*member-1)] - X1.cont[, (2*member)]
X1r.bin1 <- X1.bin[, 1:(member-1)] - X1.bin[, member]
X1r.bin2 <- X1.bin[, (member+1):(2*member-1)] - X1.bin[, (2*member)]

#ベクトル形式に変更
X1v.cont1 <- as.numeric(t(X1r.cont1))
X1v.cont2 <- as.numeric(t(X1r.cont2))
X1v.bin1 <- as.numeric(t(X1r.bin1))
X1v.bin2 <- as.numeric(t(X1r.bin2))

##データを結合
X <- data.frame(pop=POP, G, c1=X1v.cont1, c2=X1v.cont2, b1=X1v.bin1, b2=X1v.bin2)
XM <- as.matrix(X)

##分散共分散行列の設定
corM <- corrM(col=member-1, lower=-0.6, upper=0.70)   #相関行列を作成
Sigma <- covmatrix(col=member-1, corM=corM, lower=1, upper=1)   #分散共分散行列
Cov <- Sigma$covariance


##多項データを発生
##パラメータの設定
##妥当な応答変数が発生するまで繰り返す
for(i in 1:10000){
  print(i)
  
  #回帰係数のパラメ-タの設定
  b0 <- runif(member-1, 0.1, 1.7)
  b1 <- runif(c-1, 0.2, 1.3)
  b2 <- runif(cont, 0, 1.1)
  b3 <- runif(bin, -0.9, 1.2)
  b <- c(b0, b1, b2, b3)
  beta.t <- b
  
  ##相対効用を発生させる
  err <- mvrnorm(hh, rep(0, member-1), Cov)   #誤差構造
  U.mean <- matrix(XM %*% b, nrow=hh, ncol=member-1, byrow=T)   #相対効用の平均構造
  U <- U.mean + err   #誤差構造を加えた構造
  
  ##効用最大化原理に基づき選択メンバーを決定
  y <- apply(cbind(U, 0), 1, which.max)
  
  #基準メンバーが適当な人数に選ばれるまでループさせる
  if(sum(y==member) > 15 & sum(y==member) < 100) {break} else {next}
}

#選択メンバーを0、1行列に変更
Y <- matrix(0, hh, member)
for(i in 1:hh){
  Y[i, y[i]] <- 1
}

table(y)   #選択ブランドの集計
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
llike <- array(0, dim=c(R/keep))   #対数尤度の保存用


##説明変数を多次元配列化
X.array <- array(0, dim=c(member-1, ncol(X), hh))
for(i in 1:hh){
  X.array[, , i] <- XM[ID[, 2]==i, ]
}
YX.array <- array(0, dim=c(member-1, ncol(X)+1, hh))


#推定プロセスの格納配列
UM <- matrix(0, nrow=hh, ncol=member-1)
util.M <- matrix(0, nrow=hh, ncol=member-1)   

##事前分布の設定
nu <- member   #逆ウィシャート分布の自由度
V <- solve((1/10)*diag(member-1))    #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(X))  #回帰係数の平均の事前分布
Adelta <- solve(100 * diag(rep(1, ncol(X))))   #回帰係数の事前分布の分散

##サンプリング結果の保存用配列
Util <- array(0, dim=c(hh, member-1, R/keep))
BETA <- matrix(0, nrow=R/keep, length(b))
SIGMA <- matrix(0, nrow=R/keep, ncol=(member-1)^2)

##初期値の設定
#回帰係数の初期値
oldbeta <- c(runif(member-1, 0, 2), runif(c-1, 0, 1.5), runif(2, 0, 2.0), runif(2, -2, 2))   


#分散共分散行列の初期値
corM.f <- corrM(col=member-1, lower=0, upper=0)   #相関行列を作成
Sigma.f <- covmatrix(col=member-1, corM=corM.f, lower=1, upper=1)   #分散共分散行列
oldcov <- Sigma.f$covariance

#効用の平均構造の初期値
old.utilm <- matrix(XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)

#効用の初期値
old.util <- old.utilm + mvrnorm(nrow(old.utilm), rep(0, member-1), oldcov)


####マルコフ連鎖モンテカルロ法で多項プロビットモデルを推定####
for(rp in 1:R){
  
  ##選択結果と整合的な潜在効用を発生させる
  #条件付き期待値と条件付き分散を計算
  S <- rep(0, member-1)
  
  for(j in 1:(member-1)){
    MVR <- cdMVN(mu=old.utilm, Cov=oldcov, dependent=j, U=old.util)   #条件付き分布を計算
    UM[, j] <- MVR$CDmu   #条件付き期待値を取り出す
    S[j] <- sqrt(MVR$CDvar)    #条件付き分散を取り出す
    
    #潜在変数を発生させる
    #切断領域の設定
    max.u <- apply(cbind(old.util[, -j], 0), 1, max)
    max.u <- ifelse(y==member, 0, max.u)
    
    #切断正規分布より潜在変数を発生
    old.util[, j] <- ifelse(y==j, rtnorm(mu=UM[, j], sigma=S[j], a=max.u, b=100), 
                            rtnorm(mu=UM[, j], sigma=S[j], a=-100, b=max.u))
    old.util[, j] <- ifelse(is.infinite(old.util[, j]), ifelse(y==j, max.u + runif(1), max.u - runif(1)), old.util[, j])
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
  R.error <- matrix(util.v - XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)
  IW.R <- V + matrix(rowSums(apply(R.error, 1, function(x) x %*% t(x))), nrow=member-1, ncol=member-1, byrow=T)
  
  #逆ウィシャート分布の自由度を計算
  Sn <- nu + hh
  
  #逆ウィシャート分布からCovをサンプリング
  Cov_hat <- rwishart(Sn, solve(IW.R))$IW
  oldcov <- cov2cor(Cov_hat)
  
  ##潜在効用とパラメータを更新
  #潜在効用と潜在効用の平均を更新
  old.utilm <- matrix(XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)
  
  ##サンプリング結果を保存
  if(rp%%keep==0){
    print(rp)
    mkeep <- rp/keep
    Util[, , mkeep] <- old.util
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    print(round(cbind(oldcov, Cov), 2))
    print(round(rbind(oldbeta, beta.t), 2))
  }
}


####関数で推定####
Data1 <- list(p=member, y=y, X=XM)
Mcmc1 <- list(R=10000, keep=2)

#多項プロビットモデルを推定
out <- rmnpGibbs(Data=Data1,Mcmc=Mcmc1)
BETA.out <- out$betadraw
SIGMA.out <- out$sigmadraw

####推定結果の要約と適合度の確認####
burnin <- 5000/keep   #バーンイン期間

##サンプリング結果を可視化
#回帰係数のプロット
matplot(BETA[, 1:3], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 4:6], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 7:9], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 10:12], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 13:15], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")

#分散供分散行列の可視化
matplot(SIGMA[, 1:4], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 5:8], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 9:12], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 13:16], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")


##推定値の事後平均の比較
#betaの要約統計量
round(colMeans(BETA.out[burnin:nrow(BETA.out), ] / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(関数推定)の事後平均
round(colMeans(BETA[burnin:nrow(BETA), ]), 3)   #betaの事後平均
round(beta.t, 3)   #真の値
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(BETA[burnin:nrow(BETA), ], 2, sd), 2)   #事後標準偏差

#sigmaの要約統計量
round(colMeans(SIGMA.out[burnin:nrow(SIGMA.out), ]  / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(関数推定)の事後平均
round(colMeans(SIGMA[burnin:nrow(SIGMA), ]), 3)   #betaの事後平均
round(as.numeric(Cov), 3)   #真の値
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, sd), 2)   #事後標準偏差
