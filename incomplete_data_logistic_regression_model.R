#####目的変数に欠測を含んだロジットモデル#####
library(MASS)
library(bayesm)
library(R2WinBUGS)
library(rstan)
library(reshape2)
library(plyr)
library(lattice)
library(ggplot2)

####データの発生####
##データの設定
col <- 15   #パラメータ数
N <- 4000   #サンプル数

##説明変数の発生
#連続変数の発生
cont <- 7   #連続変数のパラメータ数
X.cont <- matrix(rnorm(N*cont, -1, 1), N, cont)

#二値変数の発生
bin <- 3   #二値変数のパラメータ数
X.bin <- matrix(0, N, bin)
for(i in 1:bin){
  r <- runif(1, 0.2, 0.8)
  X.bin[, i] <- rbinom(N, 1, r)
}

#多値変数の発生
multi <- 5   #多値変数のパラメータ数
m <- runif(5)
X.ma <- t(rmultinom(N, 1, m))
zm <- which.min(colSums(X.ma))
X.multi <- X.ma[, -zm]

#データの結合
round(X <- data.frame(cont=X.cont, bin=X.bin, multi=X.multi), 2)

##回帰係数の設定
alpha0 <- 0.6
beta.cont <- runif(cont, 0, 0.65)
beta.bin <- runif(bin, -0.7, 0.9)
beta.multi <- runif(multi-1, -0.8, 1.2)
betaT <- c(alpha0, beta.cont, beta.bin, beta.multi)

##応答変数の発生
#確率の計算
logit <- alpha0 + as.matrix(X) %*% betaT[-1]   #ロジット
P <- exp(logit)/(1+exp(logit))   #確率の計算
hist(P, col="grey", main="確率の分布")

#ベルヌーイ乱数で応答変数を発生
Y <- rbinom(N, 1, P)
round(cbind(Y, P), 2)   #応答変数と確率の比較


##教師データをランダムな欠損に従って削除する
#ランダムな欠損を仮定
alpha.na <- 0.8
beta1.na <- runif(cont, 0, 0.7)
beta2.na <- runif(bin, -0.8, 1.4)

#ロジットと確率の計算
logit.na <- alpha.na + as.matrix(X[, 1:cont]) %*% beta1.na + as.matrix(X.bin) %*% beta2.na
P.na <- exp(logit.na)/(1+exp(logit.na))

#ベルヌーイ乱数で応答変数を欠損させる
z.na <- rbinom(N, 1, P.na)

#欠損ベクトルの作成
Y.na <- ifelse(z.na==1, NA, Y)
YZ <- data.frame(Y, z.na, Y.na)

####EMアルゴリズムで欠損のあるデータを推定####
##ロジスティック回帰モデルの対数尤度を定義
loglike <- function(b, X, Y){
  #パラメータの設定
  alpha <- b[1]
  beta <- b[2:(col)]
  
  #尤度を定義して合計する
  logit <- alpha + as.matrix(X) %*% as.vector(beta) 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}


b <- beta
Z <- z
##完全データでのロジスティック回帰モデルの対数尤度
fr <- function(b, X.comp, X.na, Y.comp, Z, k){
  beta0 <- b[1]
  beta <- b[2:(k+1)]
  
  #完全データの対数尤度
  logit.comp <- beta0 + as.matrix(X.comp) %*% beta 
  P.comp <- exp(logit.comp) / (1 + exp(logit.comp))
  LLS.comp <- Y.comp*log(P.comp) + (1-Y.comp)*log(1-P.comp)  
  LL.comp <- sum(LLS.comp)
  
  #不完全データの対数尤度
  logit.na <- beta0 + as.matrix(X.na) %*% beta
  P.na <- exp(logit.na) / (1 + exp(logit.na))
  LLs.na <- Z[, 1]*log(P.na) + Z[, 2]*log(1-P.na)
  LL.na <- sum(LLs.na)
  LL <- LL.comp + LL.na
  return(LL)
}

##観測データでの尤度と潜在変数zの計算
obsll <- function(b, X.comp, Y.comp, X.na, r, N.na, k){
  beta0 <- b[1]
  beta <- b[2:(k+1)]
  
  #完全データの対数尤度
  logit.comp <- beta0 + as.matrix(X.comp) %*% beta 
  P.comp <- exp(logit.comp) / (1 + exp(logit.comp))
  LLS.comp <- Y.comp*log(P.comp) + (1-Y.comp)*log(1-P.comp)  
  LL.comp <- sum(LLS.comp)
  
  #不完全データの尤度を計算
  logit.na <- beta0 + as.matrix(X.na) %*% beta
  P.na <- exp(logit.na) / (1 + exp(logit.na))
  
  #観測データの対数尤度と潜在変数zの計算
  #混合率
  R <- matrix(r, N.na, 2, byrow=T)
  
  #潜在変数zの計算
  LLr <- R * cbind(P.na, 1-P.na)
  z0 <- matrix(apply(LLr, 1, sum), N.na, 2)   #zの分母
  z1 <- LLr/z0   #zの計算
  
  #観測データの対数尤度
  LL.na <- sum(log(apply(matrix(r, N.na, 2, byrow=T) * cbind(P.na, 1-P.na), 1, sum)))
  LLobz <- LL.na
  rval <- list(LLobz=LLobz, z1=z1)
  return(rval)
}

##データの設定
#データを完全データと欠損データにソート
index.na <- subset(1:nrow(YZ), is.na(YZ$Y.na)==TRUE)

#完全データ
YZ.comp <- YZ[-index.na, ]
Y.comp <- Y[-index.na]
X.comp <- X[-index.na, ]

#欠損のあるデータ
YZ.na <- YZ[index.na, ]
Y.na <- Y[index.na]
X.na <- X[index.na, ]

#結合応答変数
YZ <- rbind(YZ.comp, YZ.na)


##EMアルゴリズムの設定と初期値の設定
iter <- 0
dl <- 100   #EMステップでの対数尤度の初期値を設定
tol <- 0.1

#パラメータの初期値の設定
#対数尤度を最大化する
b0 <- c(rep(0, ncol(X.comp)+1))   #初期パラメータの設定
res <- optim(b0, loglike, gr=NULL, X=X.comp, Y=Y.comp, method="BFGS", hessian=TRUE, control=list(fnscale=-1))
beta <- res$par + runif(length(res$par), 0.25, 0.25)
r <- c(0.5, 0.5)   #混合率の初期値

#観測データの尤度と潜在変数zの初期値
obsllz <- obsll(b=beta, X.comp=X.comp, Y.comp=Y.comp, X.na=X.na, r=r, N.na=nrow(X.na), k=ncol(X.comp))
LL1 <- obsllz$LLobz
z <- obsllz$z1
z

####EMアルゴリズムによる不完全データロジスティック回帰モデル####
##完全データでの回帰係数を推定(Mステップ)
while(dl >= tol){   #dlがtol以上なら繰り返す
  res <- optim(beta, fr, X.comp=X.comp, Y.comp=Y.comp, X.na=X.na, Z=z, k=ncol(X.comp), method="BFGS", 
               hessian=FALSE, control=list(fnscale=-1))
  
  beta <- res$par
  r <- apply(z, 2, sum) / nrow(X.na)   #混合率の計算
  
  ##Eステップ
  obsllz <- obsll(b=beta, X.comp=X.comp, Y.comp=Y.comp, X.na=X.na, r=r, N.na=nrow(X.na), k=ncol(X.comp))
  LL <- obsllz$LLobz
  z <- obsllz$z1
  
  iter <- iter+1
  dl <- abs(LL- LL1)
  LL1 <- LL
  print(LL)
}

round(beta, 3)
round(betaT, 3)

round(cbind(YZ.na, z), 3)
