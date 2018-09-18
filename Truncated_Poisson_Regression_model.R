#####切断ポアソン回帰モデル#####
options(warn=0)
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(2506787)

####データの発生####
##データの設定
hh <- 30000   #サンプル数
k <- 11   #説明変数数


##素性ベクトルを生成
k1 <- 3; k2 <- 4; k3 <- 5
x1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
x2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hh, 1, pr)
}
x3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #データを結合
k <- ncol(x)   #説明変数数

##応答変数の生成
repeat {
  #パラメータの生成
  beta <- betat <- c(0.5, rnorm(k-1, 0, 0.5))
  
  #切断ポアソン分布から応答変数を生成
  lambda <- as.numeric(exp(x %*% beta))   #期待値
  y <- rtpois(hh, lambda, a=0, b=Inf)
  
  if(max(y) > 15 & max(y) < 30){
    break
  }
}
hist(y, breaks=25, col="grey", main="アクセス頻度の分布", xlab="アクセス頻度")


####最尤法で切断ポアソン回帰モデルを推定####
##切断ポアソン回帰モデルの対数尤度
loglike <- function(beta, y, x, y_lfactorial){
  lambda <- as.numeric(exp(x %*% beta))   #期待値
  LL <- sum(y*log(lambda) - lambda - log(1-exp(-lambda)) - y_lfactorial)   #対数尤度関数
  return(LL)
}

##切断ポアソン回帰モデルを準ニュートン法で最尤推定
#パラメータを推定
y_lfactorial <- lfactorial(y)   #yの対数階乗
beta <- rep(0, ncol(x))   #初期値
res <- optim(beta, loglike, gr=NULL, y, x, y_lfactorial, method="BFGS", hessian=TRUE,   #準ニュートン法
             control=list(fnscale=-1, trace=TRUE))

#推定結果
beta <- res$par
rbind(beta, betat)   #真のパラメータ
(tval <- beta/sqrt(-diag(solve(res$hessian))))   #t値
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(hh)*length(beta)) #BIC

#観測結果と期待値の比較
lambda <- exp(x %*% beta)
mu <- lambda*exp(lambda) / (exp(lambda) - 1)   #期待値
round(data.frame(y, mu, lambda), 3)   #観測結果との比較

##ポアソン回帰との比較
out <- glm(y ~ x[, -1], family="poisson")
rbind(tpois=beta, pois=as.numeric(out$coefficients))   #回帰係数
res$value; as.numeric(logLik(out))   #対数尤度
sum((y - mu)^2); sum((y - as.numeric(out$fitted.values))^2)   #二乗誤差


####ハミルトニアンモンテカルロ法で切断ポアソン回帰モデルを推定####
##対数事後分布を計算する関数
loglike <- function(beta, y, x, inv_tau, y_lfactorial){
  
  #切断ポアソン回帰モデルの対数尤度
  lambda <- as.numeric(exp(x %*% beta))   #期待値
  Lho <- sum(y*log(lambda) - lambda - log(1-exp(-lambda)) - y_lfactorial)   #対数尤度関数
  
  #多変量正規分布の対数事前分布
  log_mvn <- -1/2 * as.numeric(beta %*% inv_tau %*% beta)
  
  #対数事後分布
  LL <- Lho + log_mvn
  return(LL)
}

##HMCでパラメータをサンプリングするためのサンプリング
#切断ポアソン回帰モデルの対数事後分布の微分関数


lambda <- as.numeric(exp(x %*% beta))   #期待値
lambda_exp <- exp(-lambda)

y*1/lambda - lambda_exp / (1-lambda_e







