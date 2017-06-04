#####負の多項分布モデル#####
library(MASS)
library(vcd)
library(gtools)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(534)
N <- 5000
k <- 20

#ガンマ分布からパラメータλを発生させる
theta.t <- matrix(0, nrow=k, ncol=2)
Y <- matrix(0, nrow=N, ncol=k)
lambda <- matrix(0, nrow=N, ncol=k)

for(i in 1:k){
  #ガンマ分布からlambdaを発生
  theta.t[i, ] <- c(runif(1, 0.001, 3), 1.5)
  lambda[, i] <- rgamma(N, shape=theta.t[i, 1], scale=theta.t[i, 2])
  
  #lambdaからポアソン乱数を発生させる
  poison <- rpois(N, lambda[, i])  
  Y[, i] <- poison
  
}
#発生させたデータを確認
round(data.frame(Y=Y, lambda=lambda), 1)
summary(Y)
round(apply(Y, 2, sd), 2)

#データの分布
hist(Y[, 1], breaks=25, col="grey", xlab="頻度", main="負の二項分布")   #Yの分布
hist(lambda[, 1], breaks=25, col="grey", xlab="頻度", main="lambdaの分布")   #lambdaの分布


##浸透率と頻度を計算
#浸透率の計算
perm_rate <- c()
for(i in 1:k){
 perm_rate <- c(perm_rate, mean(ifelse(Y[, i]==0, 0, 1)))
}

#頻度の計算
summary(rowSums(Y))
freq.all <- colSums(Y)/N   #全体での頻度
freq.cond <- colSums(Y)/colSums(apply(Y, 2, function(x) ifelse(x > 0, 1, 0)))

#浸透率と頻度の結合して比較
round(data.frame(perm_rate, freq.all, freq.cond), 3)   
plot(perm_rate, freq.cond, xlab="浸透率", ylab="購買頻度")


####負の多項分布モデルを最尤推定####
fr <- function(b0, Y, N, k){
  #パラメータの設定
  alpha <- exp(b0[1:k])
  scale <- exp(b0[k+1])
  
  #対数尤度を計算
  LL <- sum(
            Y * log(scale/(1+scale)) +
            matrix(alpha*log(1/(1+scale)), nrow=N, ncol=k, byrow=T) +
            log(gamma(Y+matrix(alpha, nrow=N, ncol=k, byrow=T))) -
            lfactorial(Y) - log(matrix(gamma(alpha), nrow=N, ncol=k, byrow=T)) 
            )
  return(LL)
}

##対数尤度を最大化
b0 <- c(runif(k, 0.1, 2), 1.0)   #初期パラメータの設定
res <- optim(b0, fr, gr=NULL, Y, N, k, method="BFGS", hessian=TRUE, control=list(fnscale=-1))   #準ニュートン法で最適化

####推定結果と可視化####
##推定されたパラメータ
b <- res$par
round(exp(b), 2)   #推定結果
round(c(theta.t[, 1], theta.t[1, 2]), 2)   #真のパラメータ

alpha <- exp(b[1:k])
scale <- exp(b[k+1])

##適合度とAIC
res$value   #最大化された対数尤度
(tval <- b/sqrt(-diag(solve(res$hessian))))   #t値
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(length(poison))*length(b))   #BIC


##推定結果から浸透率と購買頻度を計算
#浸透率の計算
rate.nmd <- (1-(scale+1)^-alpha) / (1-(scale+1)^-sum(alpha))  
round(data.frame(rate.nmd, perm_rate), 3)

#購買頻度の計算
freq.nmd <- (alpha*scale) / ((1-(scale+1)^-sum(alpha))*rate.nmd)   #購買頻度
round(data.frame(freq.nmd, freq.cond), 2)


##浸透率と購買頻度の関係を可視化
#浸透率と購買頻度の関係を関数スムージングで推定
lo <- loess(freq.nmd ~ rate.nmd)
x <-  seq(min(rate.nmd), max(rate.nmd), length=500)
pred.lo <- predict(lo, x)

#関数スムージング(基準線)と浸透率と購買頻度の実測値をプロット
plot(perm_rate, freq.cond, xlab="観測された浸透率", ylab="観測された購買頻度", main="浸透率と購買頻度の関係", 
     pch=3, cex=1.25)
lines(x, pred.lo, type="l", col=2)

