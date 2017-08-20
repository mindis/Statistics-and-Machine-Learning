#####対数ロジスティック加速故障モデル#####
library(MASS)
library(survival)
library(actuar)
library(FAdist)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
hh <- 5000   #サンプル数
censor_time <- 150   #打ち切り時間

####説明変数の発生####
#レベルの対数
lv.weib <- round(rweibull(hh*2, 1.8, 280), 0)
index.lv <- sample(subset(1:length(lv.weib), lv.weib > 80), hh)
lv <- scale(lv.weib[index.lv])

#スコアの対数
score.norm <- exp(rnorm(hh*2, 12.5, 0.5))
index.score <- sample(subset(1:length(score.norm), score.norm > 150000), hh)
score <- scale(score.norm[index.score])
SCORE <- score

#どのメンバーの勧誘回だったか
member <- 9
prob <- 1/(member)
scout <- t(rmultinom(hh, 2, rep(prob, member)))


#メンバーで勧誘が重複しなくなるまで乱数を発生させ続ける
for(i in 1:10000){
  if(max(scout)==1) break
  index.scout <- subset(1:nrow(scout), apply(scout, 1, max) > 1)
  scout[index.scout, ] <- t(rmultinom(length(index.scout), 2, rep(prob, member)))
  print(i)
}
SCOUT <- scout[, -which.min(colSums(scout))]


#平均ガチャガチャ経過時間
time.weib <- round(rweibull(hh, 2.5, 35), 3)
time.avg <- scale(time.weib)

#累積ガチャ回数
gamma <- runif(hh, 40, 120)
ammout.pois <- rpois(hh, gamma)
ammout <- ammout.pois/100


##データの結合
X <- data.frame(lv, score=SCORE, scout=SCOUT, time=time.avg, ammout=ammout)
XM <- as.matrix(X)


####応答変数の発生####
##パラメータの設定
shape <- runif(1, 9.0, 10.5)
beta0 <- 4.1
beta1 <- c(runif(2, -0.5, 0.3), runif(member-1, -1.0, 0.7), runif(1, -1.0, -0.5), runif(1, -0.9 -0.6))

##応答変数の発生
scale <- exp(beta0 + XM %*% beta1)   #スケールパラメータ
y <- rllogis(hh, shape=10, scale=scale)   #対数ロジスティック分布からガチャ間隔を発生
sum(y <= censor_time)   #150日以内に収まっているユーザー数

hist(y[y <= censor_time], breaks=30, main="ガチャ間隔の分布", xlab="時間", col="grey")   #分布を可視化

##打ち切り指示変数を設定
Z <- as.numeric(y <= 150)   #打ち切り指示変数
y[Z==0] <- NA


####対数ロジスティックモデルを推定####
##対数ロジスティックモデルの対数尤度
loglike <- function(x, y, Z, X, censor_time){
  #パラメータの設定
  shape <- x[1]
  theta0 <- x[2]
  theta1 <- x[3:(ncol(X)+2)]
  
  #対数尤度を計算
  scale <- exp(theta0 + X %*% theta1)   #スケールパラメータの計算
  LL.f <- log(dllogis(y[Z==1], shape=shape, scale=scale[Z==1]))   #非打ち切りデータの対数尤度
  LL.S <- log(1 - pllogis(censor_time, shape=shape, scale=scale[Z==0]))   #右側打ち切りデータの対数尤度
  LL <- sum(LL.f) + sum(LL.S)   #対数尤度の和
  return(LL)
}

##準ニュートン法で対数尤度を最大化
for(i in 1:1000){
  #初期値の設定
  x <- c(runif(1, 5.0, 10.5), runif(1, 2.0, 4.5), runif(2, -0.8, 0.8), runif(member-1, -1.2, 1.2), runif(1, -1.2, -0.3), 
         runif(1, -1.2 -0.3))
  
  #準ニュートン法で対数尤度を最大化
  res <- try(optim(x, loglike, y=y, Z=Z, X=XM, censor_time=censor_time, method="BFGS", hessian=TRUE, 
                   control=list(fnscale=-1, trace=TRUE)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break} #エラー処理
}


####推定結果の確認と要約####
##推定されたパラメータ
theta <- res$par
round(rbind(theta, thetat=c(shape, beta0, beta1)), 3)   #推定されたパラメータと真のパラメータの比較
round(exp(theta[2:length(theta)]), 3)   #パラメータを指数変換

#パラメータを格納
gamma <- theta[1]   #スケールパラメータ
theta0 <- theta[2]   #切片
theta1 <- theta[3:(ncol(XM)+2)]   #回帰係数

##適合度を計算
round(res$value, 3)   #最大化された対数尤度
round(tval <- theta/sqrt(-diag(solve(res$hessian))), 3)   #t値
round(AIC <- -2*res$value + 2*length(theta), 3)   #AIC
round(BIC <- -2*res$value + log(hh)*length(theta), 3) #BIC

##推定結果を可視化
hist(rllogis(hh, shape=gamma, scale=exp(beta0)), main="推定されたshapeパラメータおよび切片での生存時間の分布", 
     xlab="生存時間", col="grey")



