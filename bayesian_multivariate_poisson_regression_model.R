#####多変量ポアソン回帰モデル#####
library(MASS)
library(nlme)
library(glmm)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(98437)

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
  P <- X.Sigma$vector
  
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


####データの発生####
N <- 4000   #サンプル数
k <- 4   #応答変数数

####説明変数の発生####
cont1 <- 3; bin1 <- 3; multi1 <- 4
X.cont <- matrix(rnorm(N*cont1), nrow=N, ncol=cont1)
X.bin <- matrix(0, nrow=N, ncol=bin1)
X.multi <- matrix(0, nrow=N, ncol=multi1)

#二値説明変数を設定
for(i in 1:bin1){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(N, 1, p)
}

#多値説明変数を設定
p <- runif(multi1)
X.multi <- t(rmultinom(N, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))]   #冗長な変数は削除

#データを結合
X <- cbind(1, X.cont, X.bin, X.multi)


####応答変数の発生####
##パラメータの設定
#回帰パラメータの設定
beta01 <- runif(k, 0.2, 1.2)
beta02 <- matrix(runif(k*cont1, 0, 0.6), nrow=cont1, ncol=k)
beta03 <- matrix(runif(k*bin1, -0.5, 0.6), nrow=bin1, ncol=k)
beta04 <- matrix(runif(k*(multi1-1), -0.5, 0.8), nrow=multi1-1, ncol=k)
beta0 <- rbind(beta01, beta02, beta03, beta04)
rownames(beta0) <- 1:nrow(beta0)

#分散共分散パラメータの設定
corM <- corrM(col=k, lower=-0.5, upper=0.8, eigen_lower=0.01, eigen_upper=0.2)   #相関行列を作成
Sigma <- covmatrix(col=k, corM=corM, lower=0.5, upper=0.75)   #分散共分散行列
Cov0 <- Sigma$covariance
cov2cor(Cov0)

##応答変数の発生
#多変量正規分布よりlambdaを発生
mu <- X %*% beta0
lambda <- t(apply(mu, 1, function(x) mvrnorm(1, x, Cov0)))

#ポアソン分布より応答変数を発生
lambda_exp <- exp(lambda)
Y <- apply(lambda_exp, 2, function(x) rpois(N, x))

#発生させた乱数の分布
hist(Y[, 1], col="grey", main="発生させた応答変数の分布1", xlab="応答変数")
hist(Y[, 2], col="grey", main="発生させた応答変数の分布1", xlab="応答変数")
hist(Y[, 3], col="grey", main="発生させた応答変数の分布1", xlab="応答変数")
hist(Y[, 4], col="grey", main="発生させた応答変数の分布1", xlab="応答変数")


####マルコフ連鎖モンテカルロ法で多変量ポアソン回帰モデルを推定####
##ポアソン分布の対数尤度関数を定義
fr <- function(theta, Y, Y_factorial){
  #パラメータの設定
  lambda <- exp(theta)
  
  #対数尤度を計算
  LLi <- rowSums(Y*log(lambda)-lambda - Y_factorial)
  LL <- sum(LLi)
  LL_val <- list(LLi=LLi, LL=LL)
  return(LL_val)
}

##アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- c()   #対数尤度の保存用
Y_factorial <- lfactorial(Y)   #Yの階乗の対数を計算しておく

##事前分布の設定
#階層モデルの事前分布
Deltabar <- matrix(0, nrow=ncol(X), ncol=k)
Adelta <- 0.01 * diag(ncol(X))
nu <- k + ncol(X)
V <- nu * diag(k)

##サンプリング結果の保存用配列
BETA <- matrix(0, nrow=R/keep, ncol=k*ncol(X))
SIGMA <- matrix(0, nrow=R/keep, ncol=k^2)
THETA <- matrix(0, nrow=R/keep, ncol=N)

##初期値の設定
#回帰パラメータの初期値を設定
beta1 <- runif(k, -0.6, 1.2)
beta2 <- matrix(runif(k*cont1, 0, 0.6), nrow=cont1, ncol=k)
beta3 <- matrix(runif(k*bin1, -0.6, 0.8), nrow=bin1, ncol=k)
beta4 <- matrix(runif(k*(multi1-1), -0.6, 0.9), nrow=multi1-1, ncol=k)
oldbeta <- rbind(beta01, beta02, beta03, beta04)
rownames(oldbeta) <- 1:nrow(beta0)

#分散共分散パラメータの初期値を設定
corM <- corrM(col=k, lower=-0.4, upper=0.4, eigen_lower=0.01, eigen_upper=0.2)   #相関行列を作成
Sigma <- covmatrix(col=k, corM=corM, lower=0.3, upper=0.45)   #分散共分散行列
oldcov <- Sigma$covariance
cov_inv <- solve(oldcov)
cov2cor(oldcov)

#lambdaの発生
#多変量正規分布よりlambdaを発生
theta_mu <- X %*% beta0
oldtheta <- t(apply(theta_mu, 1, function(x) mvrnorm(1, x, Cov0)))


####MCMCでパラメータをサンプリング####
for(rp in 1:R){

  ##MH法でサンプルごとにthetaをサンプリング
  thetad <- oldtheta
  thetan <- thetad + 0.3 * mvrnorm(1, rep(0, k), diag(k))
  
  #事前分布の誤差を計算
  er_new <- thetan - theta_mu
  er_old <- thetad - theta_mu
  
  #対数尤度と対数事前分布を計算
  lognew <- fr(thetan, Y, Y_factorial)$LLi
  logold <- fr(thetad, Y, Y_factorial)$LLi
  logpnew <- apply(er_new, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  logpold <- apply(er_old, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  
  ##MHサンプリング
  #サンプリングを採択するかどうかを決定
  rand <- runif(N)   #一様乱数から乱数を発生
  LLind.diff <- exp(lognew + logpnew - logold - logpold)   #採択率を計算
  alpha <- ifelse(LLind.diff > 1, 1, LLind.diff)
  alpha <- matrix(alpha, nrow=N, ncol=k)
  
  #alphaに基づきbetaを採択
  oldtheta.r <- ifelse(alpha > rand, thetan, thetad)   #alphaがrandを上回っていたら採択
  adopt <- sum(oldtheta[, 1]!=oldtheta.r[, 1])/N   #採択率
  oldtheta <- oldtheta.r   #パラメータを更新
  
  
  ##多変量回帰モデルによる階層モデルのサンプリング
  out <- rmultireg(Y=oldtheta, X=X, Bbar=Deltabar, A=Adelta, nu=nu, V=V)
  oldbeta <- out$B
  oldcov <- out$Sigma
  
  #階層モデルのパラメータを更新
  cov_inv <- solve(oldcov)
  theta_mu <- X %*% oldbeta
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep  
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    
    print(rp)
    print(round(sum(lognew), 2))
    print(round(cbind(oldbeta, beta0), 3))
    print(round(cbind(oldcov, Cov0), 3))
    print(round(cbind(cov2cor(oldcov), cov2cor(Cov0)), 3))
    print(round(adopt, 3))
  }
}

matplot(BETA[, 1:5], type="l")
matplot(BETA[, 6:10], type="l")
matplot(BETA[, 11:15], type="l")
matplot(BETA[, 16:20], type="l")
matplot(BETA[, 2:10], type="l")
matplot(SIGMA[, 1:4], type="l")
matplot(SIGMA[, 5:8], type="l")
matplot(SIGMA[, 9:12], type="l")
matplot(SIGMA[, 13:16], type="l")
