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
hh <- 3000
select <- 6
k <- 300   #説明変数数

##説明変数の発生
freq <- rpois(hh, 200)   #ポアソン分布から頻度を発生
p <- rdirichlet(hh, runif(k, 0.2, 1.0))   #ディレクリ分布から出現確率を発生
X <- t(apply(cbind(freq, p), 1, function(x) rmultinom(1, x[1], x[-1])))   #多項分布から説明変数を発生


#説明変数をベクトル化
#IDのベクトル化
u.id <- rep(1:hh, rep(select, hh))
i.id <- rep(1:select, hh)
ID <- data.frame(no=1:(hh*select), u.id=u.id, i.id=i.id)

#切片のベクトル化
BP <- matrix(diag(select), nrow=hh*select, ncol=select, byrow=T)[, -select]

#説明変数のベクトル化
X_vec <- matrix(0, nrow=hh*select, ncol=ncol(X)*(select-1))

for(i in 1:hh){
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
  bo1 <- matrix(runif((select-1)*k, -1, 1), nrow=k, ncol=select-1) * matrix(rbinom(k, 1, 0.15), nrow=k, ncol=select-1)
  b01 <- as.numeric(t(ifelse(abs(bo1) > 0.25, bo1, 0)))
  b0 <- c(b00, b01)
  
  #ロジットと確率の計算
  logit <- matrix(XM_vec %*% b0, nrow=hh, ncol=select, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #多項分布から応答変数を発生
  y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
  if(min(colMeans(y)) > 0.1 & max(colMeans(y)) < 0.4) break
}

####マルコフ連鎖モンテカルロ法でL1正則化多項ロジットモデルを推定####
fr <- function(theta, lambda, y, X, hh, select){
  
  #ロジットと確率の計算
  logit <- matrix(X %*% theta, nrow=hh, ncol=select, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #罰則付き対数尤度を定義
  LLi <- rowSums(y*log(Pr)) - lambda*sum(abs(theta))
  LL <- sum(LLi)
  return(LL)
}

##対数尤度を最大化
theta <- runif((select-1)*(k+1), -0.2, 0.2)
lambda <- 0.001
res <- try(optim(theta[1:15], fr, gr=NULL, lambda, y, XM_vec[, 1:15], hh, select, method="BFGS", hessian=FALSE, 
                 control=list(fnscale=-1, trace=TRUE)), silent=TRUE)
res$par
b0
fr(b0, lambda, y, XM_vec, hh, select)


