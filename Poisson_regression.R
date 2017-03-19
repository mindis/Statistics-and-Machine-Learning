#####ポアソン回帰モデル#####
####シミュレーションデータの発生####
#set.seed(438297)
options(digits = 5)
col <- 10   #変数数
n <- 5000   #サンプル数
rn <- rnorm(n, 0, 1.5)
x1 <- matrix(rn, nrow=n/col, ncol=col, byrow=T)
x2 <- x1[, 1] > -0.3 & x1[, 2] < 0 
Xm <- as.data.frame(cbind(x1, x2))
round(head(Xm, 10), 3)
n <- nrow(Xm)   #サンプル数

#乱数発生用行列
X <- cbind(1, Xm)
names(X)[1] <- "intercept"
names(X)[12] <- "v11"
names(X)

#ポアソン乱数を発生させる
beta <- c(0.9, rnorm(11, 0, 0.4))
yr <- as.matrix(X) %*% as.vector(beta)
lambda <- exp(yr)   #ポアソン平均
round(lambda, 3)   #データの確認
y <- rpois(n, lambda)   #ポアソン乱数の発生
round(y, 3)
mean(y)
#データを結合
Xy <- cbind(y, Xm)[, -2]
round(head(Xy, 20), 3)
X1 <- X[, -1]

####ポアソン回帰モデル推定する
##対数尤度の設定
#パラメータの設定
fr <- function(b, x, y){
  alpha <- b[1]
  beta <- b[2:12]
  
  #尤度を定義して合計する
  lambda <- exp(alpha + as.matrix(x) %*% as.vector(beta))
  lambdal <- log(lambda)   #リンク関数
  LLi <- y*log(lambda)-lambda - lfactorial(y)
  LL <- sum(LLi)
  return(LL)
}

##対数尤度を最大化する
b0 <- c(rep(0, 12))   #初期パラメータの設定
res <- optim(b0, fr, gr=NULL, x=X1, y=y, method="BFGS", hessian=TRUE, control=list(fnscale=-1))

#結果
(b <- res$par)   #推定されたパラメータ
beta   #真のパラメータ
(tval <- b/sqrt(-diag(solve(res$hessian))))   #t値
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(n)*length(b))   #BIC

#適合度
lambdae <- round(exp(b[1] + as.matrix(X1) %*% as.vector(b[2:12])), 3)   #推定された計数の期待値
round(cbind(y, lambda, lambdae, y-lambdae), 3)   #真の計数との誤差
round(rbind(beta, b, beta-b), 5)   #パラメータの誤差

#関数を使うなら
res2 <- glm(y ~ ., data=cbind(y, X1), family = poisson(link=log))
summary(res2)
round(coef(res2), 3)
round(b, 3)
