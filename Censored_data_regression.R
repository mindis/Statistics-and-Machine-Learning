#####打ち切りデータのモデリング#####
library(MASS)
library(reshape2)
library(plyr)

####片側打ち切りモデル####
####データの発生####
n <- 3000   #サンプル数
p <- 15   #説明変数数
b <- runif(p, -1.5, 4.5)   #回帰係数
b0 <- 8.4   #切片
sigma <- 8   #標準偏差
X <- matrix(runif(n*p, -1.0, 5.0), n, p)   #説明変数

betaT <- c(b, b0, sigma)

#真のデータを発生
D <- trunc(X %*% b + b0 + rnorm(n, 0, sigma))   #真の需要関数
S <- trunc(X %*% b + b0 + runif(n, 0, 2.5))   #真の供給関数

#購買データを発生(需要が供給を上回っている場合供給を購買データとする)
B <- ifelse(D > S, S, D)
(BDS <- data.frame(B, D, S))

#打ち切りデータの指示変数を作成
z1 <- subset(1:n, BDS$D < BDS$S)   #需要が満たされているデータ
z2 <- subset(1:n, BDS$D > BDS$S)   #需要が供給を上回っているデータ
length(z1); length(z2)

####打ち切りデータモデルを推定#####
##対数尤度関数の定義
fr <- function(theta, D, B, X, z1, z2, p){
  beta <- theta[1:p]
  beta0 <- theta[p+1]
  sigma <- exp(theta[p+2])   #非負制約
  Xb <- beta0 + as.matrix(X) %*% as.vector(beta)
  
  #非打ち切りデータの尤度
  L1 <- sum(-log(sigma^2) - ((D - Xb)[z1])^2 / sigma^2)
  
  #打ち切りデータの尤度
  Lt <- 1-pnorm((B - Xb)[z2] / sigma)
  if(sum(Lt==0)!=0){i <- subset(1:length(Lt), Lt==0); Lt[i] <- 10^-100}
  L2 <- sum(log(Lt))   
  
  #対数尤度を合計
  LL <- sum(L1 + L2)
  return(LL)
}

##初期値の設定
fitf <- lm(B ~ X)
betaf <- fitf$coef[2:16]
betaf0 <- fitf$coef[1]
betaff <- as.numeric(c(betaf, betaf0, 1))

##対数尤度の最大化
fit <- optim(betaff, fn=fr, gr=NULL, D, B, X, z1, z2, p, 
             method="BFGS", hessian=T, control=list(fnscale=-1))

##結果と統計量
round(b <- c(fit$par[1:16], exp(fit$par[17])), 3)   #推定されたパラメータ
round(betaT, 3)   #真の係数


c(b[1:16], log(b[17]))/sqrt(-diag(solve(fit$hessian)))   #t値
(AIC <- -2*fit$value + 2*length(fit$par))   #AIC
(BIC <- -2*fit$value + log(nrow(X))*length(fit$par))   #BIC


####両側打ち切りモデル####
####データの発生####
n <- 3000   #サンプル数
p <- 20   #説明変数数
b <- c(rnorm(12, 0.21, 0.13), rnorm(8, -0.20, 0.11))   #回帰係数
b0 <- 0.6   #切片
sigma <- 0.5   #標準偏差
X1 <- matrix(trunc(rnorm(n*p, 3, 1)), n, p)   #説明変数
X2 <- ifelse(X1 > 5, 5, X1)   #5以上の数値を5にする
X <- ifelse(X2 < 1, 1, X2)   #1以下の数値を1にする

#真のデータの発生
D <- b0 + X %*% b + rnorm(n, 0, sigma)   #真の需要関数
min(D)
max(D)


X %*% b
