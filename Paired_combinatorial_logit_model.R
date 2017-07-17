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
hh <- 10000   #サンプル数
k <- 5   #選択数

####説明変数の発生####
#条件付きの説明変数の発生
X1.cont <- matrix(rnorm(2*(hh*k), 0, 1), nrow=hh, ncol=k*2)

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
b0 <- runif(k-1, -1.2, 2.0)
b1 <- runif(2, 0, 1.3)
b2 <- runif(2, -1.2, 1.4)
b3 <- runif(k-1, 0, 1.5)
b4 <- runif(k-1, -1.3, 1.5)
b <- c(b0, b1, b2, b3, b4)
beta.t <- b

#類似度パラメータの設定
tau <- matrix(0, nrow=k, ncol=k-1) 
Cov <- corrM(col=k, lower=0, upper=0.8)   #類似度パラメータを発生

for(i in 1:k){
  r <- 1:k
  tau[i, ] <- Cov[i, r[-i]]
}

#PCLモデルに基づく確率の計算
logit <- matrix(XM %*% b, nrow=hh, ncol=k, byrow=T)   #効用関数の設定
P1 <- matrix(0, nrow=hh, ncol=k)   #確率の分子部分のパラメータの格納用配列
P2 <- matrix(0, nrow=hh, ncol=k-1)   #確率の分母部分のパラメータの格納用配列

#確率計算の必要部分を計算
for(i in 1:k){
  #分子部分の計算
  tau.m1 <- matrix(1-tau[i, ], nrow=hh, ncol=k-1, byrow=T)
  tau.m2 <- matrix(-tau[i, ], nrow=hh, ncol=k-1, byrow=T)
  logit1 <- matrix(logit[, i], nrow=hh, ncol=k-1)
  logit2 <- logit[, (1:k)[-i]]
  P1[, i] <- rowSums(tau.m1 * (exp(logit1/tau.m1) + exp(logit2/tau.m1))^tau.m2 * exp(logit1/tau.m1))
  
  if(i==k) {break}
  r <- (1:k)[-(1:i)]
  tau.m3 <- matrix(cbind(0, 1-tau)[i, r], nrow=hh, ncol=length(r), byrow=T)
  logit3 <- matrix(logit[, i], nrow=hh, ncol=length(r))
  logit4 <- logit[, ((i+1):k)]
  P2[, i] <- rowSums(tau.m3 * (exp(logit3/tau.m3) + exp(logit4/tau.m3))^tau.m3)
}

#確率の計算
Pr <- P1 / matrix(rowSums(P2), nrow=hh, ncol=k)

Pr1 <- matrix(exp(XM %*% b), hh, k, byrow=T) / rowSums(matrix(exp(XM %*% b), hh, k, byrow=T))
round(cbind(Pr, Pr1), 3)
round(Pr, 3)

##応答変数を発生
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
colSums(Y)
round(data.frame(Y=Y, P=Pr), 3)

####Paired combinatorial logit modelを最尤推定####
loglike <- function(x, Y, X, k, P1, P2){
  
  #パラメータの設定
  b <- x[1:ncol(X)]
  tau.v <- abs(x[(ncol(X)+1):length(x)])
  
  #類似度パラメータの設定
  tau1 <- matrix(0, nrow=k, ncol=k-1)
  K <- diag(k)
  K[lower.tri(K)] <- tau.v
  K_t <- t(K) + K - diag(k)
  
  for(i in 1:k){
    r <- 1:k
    tau[i, ] <- K_t[i, r[-i]]
  }
  
  #PCLモデルに基づく確率の計算
  logit <- matrix(X %*% b, nrow=hh, ncol=k, byrow=T)   #効用関数の設定
  
  #確率計算の必要部分を計算
  for(i in 1:k){
    #分子部分の計算
    tau.m1 <- matrix(1-tau[i, ], nrow=hh, ncol=k-1, byrow=T)
    tau.m2 <- matrix(-tau[i, ], nrow=hh, ncol=k-1, byrow=T)
    logit1 <- matrix(logit[, i], nrow=hh, ncol=k-1)
    logit2 <- logit[, (1:k)[-i]]
    P1[, i] <- rowSums(tau.m1 * (exp(logit1/tau.m1) + exp(logit2/tau.m1))^tau.m2 * exp(logit1/tau.m1))
    
    if(i==k) {break}
    r <- (1:k)[-(1:i)]
    tau.m3 <- matrix(cbind(0, 1-tau)[i, r], nrow=hh, ncol=length(r), byrow=T)
    logit3 <- matrix(logit[, i], nrow=hh, ncol=length(r))
    logit4 <- logit[, ((i+1):k)]
    P2[, i] <- rowSums(tau.m3 * (exp(logit3/tau.m3) + exp(logit4/tau.m3))^tau.m3)
  }
  
  #対数尤度を計算
  LLl <- rowSums(Y*log(P1)) - log(rowSums(P2))
  LL <- sum(LLl)
  return(LL)
}

##多項ロジットモデルの対数尤度関数
LL_logit <- function(x, X, Y, hh, k){
  #パラメータの設定
  theta <- x
  
  #効用関数の設定
  U <- matrix(X %*% theta, nrow=hh, ncol=k, byrow=T)
  
  #対数尤度の計算
  d <- rowSums(exp(U))
  LLl <- rowSums(Y * U) - log(d)
  LL <- sum(LLl)
  return(LL)
}


##制約付き準ニュートン法でPCLモデルを最尤推定
#パラメータの初期値初期値を多項ロジットモデルで決定
#準ニュートン法で最尤推定
theta <- runif(ncol(XM), -1, 1)
res.z <- optim(theta, LL_logit, gr=NULL,  X=XM, Y=Y, hh=hh, k=k, method="BFGS", hessian=FALSE,
               control=list(fnscale=-1))   
theta <- res.z$par   #パラメータの初期値

P1 <- matrix(0, nrow=hh, ncol=k)   #確率の分子部分のパラメータの格納用配列
P2 <- matrix(0, nrow=hh, ncol=k-1)   #確率の分母部分のパラメータの格納用配列

for(i in 1:1000){
  print(i)
  
  #類似度パラメータの初期値
  tau.first1  <- corrM(col=k, lower=0, upper=0.7)   #類似度パラメータを発生
  tau.first <- tau.first1[lower.tri(tau.first1)]

  #初期値をベクトルに結合
  x <- c(theta, tau.first)
  
  #パラメータの上限値と下限値を設定
  l <- c(theta+2, rep(0, sum(1:(k-1))))
  u <- c(theta-2, rep(1, sum(1:(k-1))))

  #準ニュートン法で最尤推定
  res <- try(optim(x, loglike, gr=NULL, X=XM, Y=Y, k=k, P1=P1, P2=P2, method="SANN", hessian=TRUE, control=list(fnscale=-1)), 
             silent=FALSE)
  if(class(res)=="try-error") {next} else {break}   #エラー処理
}

####パラメータ推定結果と要約####
##推定されたパラメータと真のパラメータの比較
round(beta <- res$par, 3)
round(c(beta.t, Cov[lower.tri(Cov)]), 3)

#類似度パラメータを相関パラメータに変換
rho.m <- diag(k)
rho <- beta[(ncol(XM)+1):length(beta)]
rho.m[lower.tri(rho.m)] <- rho
round(rho.M <- abs(t(rho.m) + rho.m - diag(k)), 3)
round(Cov, 3)

#最大化された対数尤度とAIC
round(res$value, 3)

round(tval <- beta/sqrt(-diag(solve(res$hessian))), 3)   #t値
round(AIC <- -2*res$value + 2*length(beta), 3)   #AIC
round(BIC <- -2*res$value + log(hh)*length(beta), 3)   #BIC

