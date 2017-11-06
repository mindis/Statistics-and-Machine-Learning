#####ベイジアン正則化多項ロジットモデル#####
library(MASS)
library(matrixStats)
library(flexmix)
library(glmnet)
library(FAdist)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
hh <- 4000
select <- 6
k <- 100   #説明変数数

##説明変数の発生
freq <- rpois(hh, 150)   #ポアソン分布から頻度を発生
p <- rdirichlet(hh, runif(k, 0.2, 1.0))   #ディレクリ分布から出現確率を発生
X <- scale(t(apply(cbind(freq, p), 1, function(x) rmultinom(1, x[1], x[-1]))))   #多項分布から説明変数を発生


#説明変数をベクトル化
#IDのベクトル化
u.id <- rep(1:hh, rep(select, hh))
i.id <- rep(1:select, hh)
ID <- data.frame(no=1:(hh*select), u.id=u.id, i.id=i.id)

#切片のベクトル化
BP <- matrix(diag(select), nrow=hh*select, ncol=select, byrow=T)[, -select]
x <- diag(1, select)[, -select]

#説明変数のベクトル化
X_vec <- matrix(0, nrow=hh*select, ncol=ncol(X)*(select-1))

for(i in 1:hh){
  print(i)
  x_diag0 <- c()
  for(j in 1:ncol(X)){
    x_diag0 <- cbind(x_diag0, diag(X[i, j], select-1))
  }
  X_vec[ID$u.id==i, ] <- rbind(x_diag0, 0)
}
XM_vec <- cbind(BP, X_vec)


##応答変数の発生
for(i in 1:1000){
  print(i)
  
  #パラメータの設定
  b00 <- runif(select-1, -1.0, 1.0)
  b01 <- matrix(runif((select-1)*k, -1, 1), nrow=k, ncol=select-1) * matrix(rbinom(k*(select-1), 1, 0.25), nrow=k, ncol=select-1)
  b0 <- rbind(b00, b01)

  #ロジットと確率の計算
  logit <- cbind(cbind(1, X) %*% b0, 0)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #多項分布から応答変数を発生
  y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
  if(min(colMeans(y)) > 0.05 & max(colMeans(y)) < 0.4) break
}
rownames(b0) <- NULL


####マルコフ連鎖モンテカルロ法でL1正則化多項ロジットモデルを推定####
##多項ロジットモデルの対数尤度を定義
fr <- function(beta, y, X, hh, select){
  
  #ロジットと確率の計算
  logit <- t(x %*% t(beta))
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #対数尤度を設定
  LLi <- rowSums(y*log(Pr)) 
  return(LLi)
}

##アルゴリズムの設定
R <- 20000
keep <- 4
sbeta <- 1.5
iter <- 0

##初期値の設定
beta0 <- scale(colSums(y))
oldbeta <- mvrnorm(hh, beta0[-select], diag(0.2, select-1))
oldtheta <- matrix(0, nrow=k+1, ncol=select-1)
oldcov <- diag(0.1, select-1)
inv_cov <- solve(oldcov)
mu <- cbind(1, X) %*% oldtheta

##サンプリング結果の保存用配列
THETA <- array(0, dim=c(k+1, select-1, R/keep))
COV <- array(0, dim=c(select-1, select-1, R/keep))

##アルゴリズム推定用配列
lognew <- rep(0, hh)
logold <- rep(0, hh)
logpnew <- rep(0, hh)
logpold <- rep(0, hh)
lambda <- rep(1, select-1)
er_new <- matrix(0, nrow=hh, select-1)
er_old <- matrix(0, nrow=hh, select-1)

#IDの設定
id_vec <- c()
for(i in 1:hh){
  id_vec <- c(id_vec, rep(i, sum(ID$id==i)*select))
}

####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##多項ロジットモデルのパラメータをサンプリング
  #新しいパラメータをサンプリング
  betad <- oldbeta
  betan <- betad + mvrnorm(hh, rep(0, select-1), diag(0.05, select-1))
  
  #誤差を設定
  er_new <- betan - mu
  er_old <- betad - mu
  
  #対数尤度と対数事前分布を計算
  lognew <- fr(betan, y, x, hh, select)
  logold <- fr(betad, y, x, hh, select)
  logpnew <- -0.5 * rowSums(er_new %*% inv_cov * er_new)
  logpold <- -0.5 * rowSums(er_old %*% inv_cov * er_old)
  
  #メトロポリスヘイスティング法でパラメータの採択を決定
  rand <- runif(hh)   #一様分布から乱数を発生
  LLind_diff <- exp(lognew + logpnew - logold - logpold)   #採択率を計算
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=hh, ncol=select-1)
  oldbeta <- flag*betan + (1-flag)*betad   #alphaがrandを上回っていたら採択
  
  
  ##ベイジアンlassoで階層モデルの回帰パラメータをサンプリング
  for(j in 1:(select-1)){
    res <- blasso(X=X, y=oldbeta[, j], beta=oldtheta[-1, j], lambda2=lambda[j], s2=diag(oldcov)[j], T=2)
    oldtheta[, j] <- c(res$mu[2], res$beta[2, ])
    lambda[j] <- res$lambda2[2]
  }
  
  ##階層モデルの分散共分散行列を推定
  mu <- cbind(1, X) %*% oldtheta
  oldcov <- var(oldbeta - mu)
  inv_cov <- solve(oldcov)
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- oldtheta
    COV[, , mkeep] <- oldcov
    print(rp)
    print(sum(lognew))
    print(round(lambda, 3))
    print(round(t(cbind(oldtheta, b0)[1:15, ]), 3))
  }
}

matplot(t(THETA[1, , ]), type="l")
matplot(t(COV[1, , ]), type="l")

THETA[1, , 1]
