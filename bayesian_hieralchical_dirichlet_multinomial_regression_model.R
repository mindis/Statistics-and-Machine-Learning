#####階層ベイズディクレリ多項回帰モデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(gtools)
library(matrixStats)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(6812)

####データの発生####
hh <- 1000   #ユーザー数
pt <- 10   #観測期間
hhpt <- hh*pt   #総サンプル数
select <- 10   #選択肢数
st <- 10   #基準選択肢

##IDの設定
id <- rep(1:hh, rep(pt, hh))
time <- rep(1:pt, hh)
ID <- data.frame(no=1:hhpt, id, time)

#ベクトル型IDの設定
u.id <- rep(1:hh, rep(pt*select, hh))
i.id <- rep(1:select, pt*hh)
time <- rep(rep(1:pt, rep(select, pt)), hh)
ID_vec <- data.frame(no=1:length(u.id), id=u.id, time=time, item=i.id) 

####説明変数の発生####
##離散選択モデルの説明変数の発生
##条件付きの説明変数の発生
X1.cont <- matrix(rnorm(hhpt*select*2, 0, 1), nrow=hhpt*select, ncol=2)

X1.bin <- matrix(0, nrow=hhpt*select, ncol=3)
for(i in 1:3){
  X1.bin[, i]  <- rbinom(hhpt*select, 1, runif(1, 0.35, 0.6))
}

X1.multi <- t(rmultinom(hhpt*select, 1, c(0.2, 0.2, 0.3, 0.3)))[, -4]

#切片の設定
Pop0 <- matrix(diag(1, select), nrow=hhpt*select, ncol=select, byrow=T)
Pop <- Pop0[, -st]

#データの結合
X <- data.frame(bp=Pop, c=X1.cont, b=X1.bin, m=X1.multi)
XM <- as.matrix(X)
round(XM, 2)


##多項型の説明変数の発生
X2.cont <- rnorm(hh, 0, 1)
X2.pois <- rpois(hh, 2)
X2.bin <- matrix(rbinom(hh*2, 1, 0.4), nrow=hh, ncol=2)

#データの結合
Z <- data.frame(i=1, c=X2.cont, p=X2.pois, b=X2.bin)
ZX <- as.matrix(Z)


####応答変数の発生####
##購買数量を発生
freq <- rep(0, hhpt)
for(i in 1:hh){
  par <- rgamma(hh, 25, 1.25)
  freq[ID$id==i] <- rpois(pt, par)   #購買数
}

##階層回帰モデルからパラメータを発生
#回帰パラメータを設定
gamma00 <- c(runif(select-1, -0.5, 0.6), runif(ncol(X1.cont), 0, 0.4), runif(ncol(X1.bin)+ncol(X1.multi), -0.4, 0.4))
gamma01 <- runif(ncol(XM), 0, 0.4)
gamma02 <- runif(ncol(XM), -0.2, 0.2)
gamma03 <- matrix(runif(ncol(XM*2), -0.4, 0.3), nrow=2, ncol=ncol(XM))
gamma0 <- rbind(gamma00, gamma01, gamma02, gamma03)

#分散パラメータを設定
Cov0 <- diag(runif(ncol(XM), 0.15, 0.4))

#多変量正規分布からディクレリ多項回帰モデルの回帰パラメータを発生
beta_mu <- ZX %*% gamma0   #平均構造
er <- mvrnorm(hh, rep(0, ncol(XM)), Cov0)   #変量効果の誤差
beta0 <- beta_mu + er   #ディクレリ多項回帰モデルのパラメータ

##ディクレリ分布から確率を発生
Alpha_m <- matrix(0, nrow=hhpt, ncol=select)
Prob <- matrix(0, nrow=hhpt, ncol=select)

for(i in 1:hh){
  Alpha_m[ID$id==i, ] <- matrix(exp(XM[ID$id==i, ] %*% beta0[i, ]), nrow=pt, ncol=select, byrow=T)
  Prob[ID$id==i, ] <- t(apply(Alpha_m[ID$id==i, ], 1, function(x) rdirichlet(1, x)))
}

##多項分布から応答変数を発生
Y <- t(apply(cbind(freq, Prob), 1, function(x) rmultinom(1, x[1], x[-1])))
y <- as.numeric(t(Y))
Pr <- Y / matrix(rowSums(Y), nrow=hhpt, ncol=select)

####マルコフ連鎖モンテカルロ法で階層ベイズディクレリ多項回帰モデルを推定####
##ディクレリ多項回帰モデルの対数尤度関数を設定
fr <- function(theta, y, X, hh, select){
  #パラメータを設定
  beta <- theta
  
  #ディクレリ分布のパラメータの平均構造を設定
  alpha <- matrix(exp(X %*% beta), nrow=hh, ncol=select, byrow=T)
  
  #ディクレリ分布の対数尤度を和
  LLi <- log(ddirichlet(y, alpha))
  LL <- sum(LLi)
  return(LL)
}


##アルゴリズムの設定
R <- 20000   #サンプリング回数
keep <- 4   #2回に1回の割合でサンプリング結果を格納
iter <- 0

##事前分布の設定
Deltabar <- matrix(0, nrow=ncol(ZX), ncol=ncol(XM))
Adelta <- 0.01 * diag(rep(1, ncol(XM)))   #階層モデルの回帰係数の分散の事前分布
nu <- (ncol(XM)) + select   #逆ウィシャート分布の自由度
V <- nu * diag(ncol(XM))

##サンプリング結果の保存用配列
PROB <- array(0, dim=c(hhpt, ncol(XM), R/10))
BETA <- array(0, dim=c(hh, ncol(XM), R/keep))
GAMMA <- matrix(0, nrow=R/keep, ncol=length(gamma0))
SIGMA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
gc(); gc()

#棄却率と対数尤度の保存用配列
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##初期値の設定
alpha_mu <- matrix(0, nrow=hh, ncol=select)
for(j in 1:ncol(Y)){
  alpha_mu[, j] <- tapply(Y[, j] + 1, ID$id, sum)
}
freq_mu <- tapply(rowSums(Y) + 1, ID$id, mean)
Pr <- alpha_mu/rowSums(alpha_mu)

#多変量回帰モデルから回帰パラメータの初期値を設定



#ディクレリ多項回帰のパラメータを推定
#応答変数の設定
Pr <- (Y + 1)/rowSums(Y + 1)

#パラメータの初期値
phi0 <- runif(select-1, -0.5, 0.5)
phi1 <- c(runif(2, 0, 0.5), runif(3, -0.4, 0.4), runif(3, -0.3, 0.3))
phi <- c(phi0, phi1)

#対数尤度を最大化する
res <- optim(phi, fr, gr=NULL, Pr, XM, hhpt, select, method="BFGS", hessian=FALSE, control=list(fnscale=-1, trace=TRUE))
r


