#####ワイブル比例ハザードモデル#####
####Webサイトの離脱分析####
library(MASS)
library(survival)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
N <- 2000   #サンプル数
page_cnt <- 10 

##ページ閲覧回数と閲覧履歴の発生
#ページ閲覧回数の発生
lam_lower <- 5
lam_upper <- 9
p_cnt.zero <- rpois(N, runif(N, lam_lower, lam_upper))
p_cnt <- ifelse(p_cnt.zero==0, 1, p_cnt.zero)
hist(p_cnt, breaks=15, col="grey", xlab="ページ閲覧数", main="ページ閲覧数の分布")

#ページ閲覧履歴の発生
p_rate <- runif(page_cnt)
p_hist <- matrix(0, N, page_cnt)

for(i in 1:N){
  p_hist[i, ] <- t(rmultinom(1, p_cnt[i], p_rate))
}
p_hist.r <- p_hist / rowSums(p_hist)

#離脱時のページと閲覧時間を発生
p_last <- t(rmultinom(N, 1, p_rate))

##累計アクセス数の発生
#ファーストランディングかどうか
pf <- 0.5
fl <- rbinom(N, 1, pf)

#2回目以降のアクセスなら累計アクセス数を発生
index.f <- subset(1:length(fl), fl==0)
pois <- rpois(length(index.f), 3)
access_cnt <- ifelse(pois==0, 1, pois)

#2回目以降のアクセスに累計アクセス数を代入
ac <- rep(0, N)
ac[index.f] <- access_cnt 

##前回からのアクセス経過時間(単位=日)の発生
mu <- 0.5
sigma <- 0.8
t <- exp(rnorm(length(index.f), mu, sigma))

#2回目以降のアクセスの場合前回からの経過時間を代入
tp <- rep(0, N)
tp[index.f] <- t 
 
##冗長な変数を削除してデータを結合
index.h <- which.min(colSums(p_hist))
ph_hist <- p_hist[, -index.h]
ph_hist.r <- p_hist.r[, -index.h]
ph_last <- p_last[, -index.h]

X <- data.frame(f=fl, t=tp, a=ac, l=ph_last, h=ph_hist.r, c=p_cnt)   #データの結合
round(X, 3)


##パラメータの設定
#回帰モデルのパラメーター
for(i in 1:50000){
  beta.f <- runif(1, -0.2, 0.5)
  beta.t <- runif(1, 0.02, 0.065)
  beta.a <- runif(1, -0.1, 1.1)
  beta.l <- runif(page_cnt-1, -0.4, 0.75)
  beta.h <- runif(page_cnt-1, -0.7, 1.1)
  beta.c <- runif(1, 0.02, 0.06)
  betat <- c(beta.f, beta.t, beta.a, beta.l, beta.h, beta.c)
  
  #ワイブル分布のパラメータ
  alpha <- runif(1, 0.5, 0.9)   #尺度パラメータ
  scale <- runif(1, 0, 3)
  lambda <- exp(scale + as.matrix(X) %*% betat)
  
  ##ワイブル乱数の発生
  y <- rweibull(N, shape=alpha, scale=lambda)
  if(max(y) > 30 & max(y) < 180 & min(y) > 0.001) break
  print(c(round(min(y), 3), i))
}

#0.25分以下のアクセス時間は取り除く
index.t <- subset(1:length(y), y < 0.25)
Y <- y[-index.t]
max(Y); min(Y)
alpha; scale

#lambdaと説明変数の0.25分以下の部分を取り除く
lambda.w <- lambda[-index.t]
XW <- X[-index.t, ]

#経過時間の分布を確認
hist(Y, breaks=50, col="grey", xlab="経過時間", main="ワイブル比例ハザードモデルの経過時間分布")
hist(rweibull(N, shape=alpha, scale=scale), breaks=50, col="grey", xlab="経過時間", main="ワイブル分布")


##コンバージョンした場合打ち切りに設定
for(t in 1:1000){
  beta.c <- c(-1.75, betat + rnorm(length(betat), 0, 0.3))
  logit <- as.matrix(cbind(1, XW)) %*% beta.c
  p <- exp(logit)/(1+exp(logit))

#コンバージョンを発生
z <- c()
  for(i in 1:length(p)){
    zbin <- rbinom(1, 1, p[i])
    z <- c(z, zbin) 
  }
  if(sum(z) > length(Y)/10 & sum(z) < length(Y)/5) break
  print(t)
}

sum(z)   #コンバージョン数
round(cbind(Y, z), 3)
mean(Y[z==1]); median(Y[z==1])   #コンバージョンした場合の滞在時間平均および滞在時間中央値
mean(Y[z==0]); median(Y[z==0])   #コンバージョンしていない場合の滞在時間平均および滞在時間中央値

Z <- 1-z   #打ち切り指示変数に変換

####ワイブル比例ハザードモデルを最尤推定#####
##ワイブル比例ハザードモデルの対数尤度
llike <- function(theta, Y, X, Z, k){
  a <- exp(theta[1])
  g <- theta[2]
  beta <- theta[3:(k+2)]
  
  lambda <- exp(g + as.matrix(X) %*% beta)   #線形結合
  LL <- sum(Z*(log(lambda)+log(a)+(a-1)*log(Y)) - lambda*Y^a)   #対数尤度を計算
  return(LL)
}


##対数尤度を最大化
for(i in 1:10000){
  print(i)
  #初期パラメータの設定
  beta0 <- c(0, 1, runif(ncol(X), -0.5, 0.5))
  
  #準ニュートン法で対数尤度を最大化
  res <- try(optim(beta0, llike, Y=Y, X=XW, Z=Z, k=ncol(XW), method="BFGS", 
                   hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #エラー処理
  
}


####結果の確認と要約####
round(alpha.res <- exp(res$par[1]), 3)   #形状パラメータの推定値
round(scale.res <- res$par[2], 3)   #スケールパラメータ(切片)の推定値
round(theta <- res$par[3:length(res$par)], 3)   #回帰係数の推定値
round(exp(theta), 3)   #ハザード比

##統計量とAIC
round(res$value, 3)   #最大対数尤度
round(tval <- res$par/sqrt(-diag(solve(res$hessian))), 3)   #t値
round(AIC <- -2*res$value + 2*length(res$par), 3)   #AIC
round(BIC <- -2*res$value + log(N)*length(res$par), 3)   #BIC


####関数を用いてワイブル加速モデルを当てはめる####
YX <- data.frame(time=Y, Z=Z, XW)   #データの設定

#ワイブル加速モデルを推定
model2<-survreg(Surv(time, Z) ~ f+t+a+l.1+l.2+l.3+l.4+l.5+l.6+l.7+l.8+l.9+h.1+h.2+h.3+h.4+h.5+h.6+
                  h.7+h.8+h.9+c, data=YX, dist="weibull")
summary(model2)

##推定結果と関数での推定の比較
#形状パラメータの比較
round(alpha.func <- 1/model2$scale, 3)   #関数で推定
round(alpha.res, 3)   #推定結果

#スケールパラメータおよび回帰係数の比較
round(as.numeric(beta.func　<- model2$coef[1:length(model2$coef)]), 3)
round(-model2$scale*res$par[2:length(res$par)], 3)

