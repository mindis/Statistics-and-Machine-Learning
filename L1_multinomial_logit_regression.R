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
hh <- 6000
select <- 10
k <- 30   #説明変数数

##説明変数の発生
freq <- rpois(hh, 200)   #ポアソン分布から頻度を発生
p <- rdirichlet(hh, runif(k, 0.2, 1.0))   #ディレクリ分布から出現確率を発生
X <- scale(t(apply(cbind(freq, p), 1, function(x) rmultinom(1, x[1], x[-1]))))   #多項分布から説明変数を発生


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
for(i in 1:5000){
  print(i)
  
  #パラメータの設定
  b00 <- runif(select-1, -0.7, 0.7)
  bo1 <- matrix(runif((select-1)*k, 0, 0.9), nrow=k, ncol=select-1) * matrix(rbinom(k, 1, 0.3), nrow=k, ncol=select-1)
  b01 <- as.numeric(t(ifelse(abs(bo1) > 0.2, bo1, 0)))
  b0 <- c(b00, b01)
  
  #ロジットと確率の計算
  logit <- matrix(XM_vec %*% b0, nrow=hh, ncol=select, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #多項分布から応答変数を発生
  y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
  if(min(colMeans(y)) > 0.05 & max(colMeans(y)) < 0.5) break
}
round(b0, 3)
colMeans(y)

####総当り座標最適化によるL1正則化多項ロジットモデルを推定####
loglike <- function(b1, b2, lambda, y, X1, X2, hh, select){

  #ロジットと確率の計算
  logit <- matrix(X1 %*% b1 + X2 %*% b2, nrow=hh, ncol=select, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #罰則付き対数尤度を定義
  LLi <- rowSums(y*log(Pr)) - lambda*sum(abs(b1)) - lambda*sum(abs(b2))
  LL <- sum(LLi)
  return(LL)
}

#真の対数尤度
LLt <- loglike(b1=b0[1], b2=b0[-1], lambda=lambda, y=y, X1=as.matrix(XM_vec[, 1]), X2=XM_vec[, -1], hh=hh, select=select)


##対数尤度を最大化
##初期値の設定
lambda <- 0.0025
theta <- rep(0, (select-1)*(k+1))

##アルゴリズムの設定
max.iter <- 30
iter <- 1
tol <- 10
diff <- 100
L1 <- 0

##L1正則化多項ロジットモデルを推定
while(diff >= tol & iter <= max.iter){
  for(i in 1:length(b0)){
    #パラメータの設定
    b1 <- theta[i]
    b2 <- theta[-i]

    #罰則付き対数尤度を最大化
    res <- optim(b1, loglike, gr=NULL, b2, lambda, y, as.matrix(XM_vec[, i]), XM_vec[, -i], hh, select, method="Nelder-Mead",
                 hessian=FALSE, control=list(fnscale=-1))
  
    #パラメータの更新
    theta[i] <- res$par
    print(c(i, LLt, res$value))
  }
  
  #アルゴリズムのパラメータの更新
  iter <- iter+1
  LL <- res$value
  diff <- abs(LL-L1)
  L1 <- LL
}

round(cbind(theta, b0), 3)
colMeans(y)
loglike()

sum(abs(theta) >= 0.1)
