#####欠損データのある混合多変量回帰モデル#####
library(MASS)
library(nlme)
library(glmm)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(5783)

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  X.Sigma <- eigen(Sigma)
  Lambda <- diag(X.Sigma$values)
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda)
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
  return(Sigma)
}

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #固有値分解で強制的に正定値行列に修正する
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

####データの発生####
n <- 1000   #評価対象数
g <- rpois(n, rgamma(n, 13.5, 1.1))
g <- ifelse(g==0, 1, g)
hh <- sum(g)   #評価者数
k <- 8   #応答変数数

##idの設定
c.id <- rep(1:length(g), g)   #評価対象ID

u.id <- c()   #ユーザーID
for(i in 1:length(g)){ 
  u.id <- c(u.id, 1:g[i])
}

ID <- data.frame(no=1:sum(g), c.id=c.id, u.id=u.id)


####説明変数の発生####
##階層モデルの説明変数
cont1 <- 3; bin1 <- 4; multi1 <- 4
X.cont <- matrix(rnorm(hh*cont1), nrow=hh, ncol=cont1)
X.bin <- matrix(0, nrow=hh, ncol=bin1)
X.multi <- matrix(0, nrow=hh, ncol=multi1)

#二値説明変数を設定
for(i in 1:bin1){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(hh, 1, p)
}

#多値説明変数を設定
p <- runif(multi1)
X.multi <- t(rmultinom(hh, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))] #冗長な変数は削除

#データを結合
X <- cbind(1, X.cont, X.bin, X.multi)


##変量効果のデザイン行列を設定
Z <- matrix(0, nrow=sum(g), ncol=n)

for(j in 1:n){
  index <- subset(1:nrow(ID), ID$c.id==j)
  Z[index, j] <- 1 
}

####応答変数の発生####
##パラメータの設定
#分散共分散行列の設定
Cor0 <- corrM(k, -0.6, 0.9, 0.01, 0.2)   #個体内モデルの相関行列
Cov0 <- covmatrix(k, Cor0, 0.6, 0.6)$covariance   #分散共分散行列に変換
CorH <- diag(runif(k, 0.5, 0.75), k)   #階層モデルの分散共分散行列

#変量効果の設定
b.random <- matrix(0, nrow=n, ncol=k)
B.Random <- matrix(0, nrow=sum(g), ncol=k)

for(i in 1:n){
  b.random[i, ] <- mvrnorm(1, rep(0, k), CorH)
  B.Random[ID$c.id==i, ] <- matrix(b.random[i, ], nrow=g[i], ncol=k, byrow=T)
}

#階層モデルのパラメータを設定
mu_score <- rnorm(k, 3.2, 0.25)   #スコアの平均構造
b1 <- matrix(runif(k*cont1, 0, 0.6), nrow=cont1, ncol=k)
b2 <- matrix(runif(k*(bin1+multi1-1), -0.7, 0.7), nrow=bin1+multi1-1, ncol=k)

BETA <- rbind(mu_score, b1, b2)   #パラメータを結合
rownames(BETA) <- c()

##応答変数の発生
Mu <- X %*% BETA + Z %*% b.random   #平均構造
Y.comp <- Mu + mvrnorm(hh, rep(0, k), Cov0)   #平均構造+誤差

##応答変数を欠損させる
#欠損パラメータ
#変量効果のパラメータ
CorA <- diag(runif(k, 0.55, 0.8))
a.random <- matrix(0, nrow=n, ncol=k)
for(i in 1:n){a.random[i, ] <- mvrnorm(1, rep(0, k), CorA)}

#固定効果のパラメータ
a0 <- runif(k, -1.8, -0.9)   #スコアの平均構造
a1 <- matrix(runif(k*cont1, 0, 0.6), nrow=cont1, ncol=k)
a2 <- matrix(runif(k*(bin1+multi1-1), -0.9, 0.5), nrow=bin1+multi1-1, ncol=k)
alpha0 <- rbind(a0, a1, a2)   #パラメータを結合
rownames(alpha0) <- c()

#ロジットと欠損確率を計算
logit <- X %*% alpha0 + Z %*% a.random
Pr <- exp(logit)/(1+exp(logit))

#欠損変数を発生
Z.na <- 1-apply(Pr, 2, function(x) rbinom(length(x), 1, x))
Y <- Y.comp * Z.na
Y <- ifelse(Y==0, NA, Y)

#欠損数を集計
colSums(is.na(Y)); round(colMeans(is.na(Y)), 2)


####マルコフ連鎖モンテカルロ法で混合多変量回帰モデルを推定####
#アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4

##事前分布の設定
#固定効果(多変量回帰)の事前分布
Deltabar <- matrix(0, nrow=ncol(X), ncol=k)   #回帰パラメータの平均の事前分布
Adelta <- 0.01 * diag(1, ncol(X))   #回帰パラメータの分散の事前分布
nu <- (ncol(X)+1)+k   #逆ウィシャート分布の自由度
V <- nu * diag(k)   #逆ウィシャート分布のパラメータ

#変量効果の事前分布
Bbar <- rep(0, n)




