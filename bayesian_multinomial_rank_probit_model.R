#####多項-ランク同時プロビットモデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(gtools)
library(MNP)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
library(lattice)

####多変量正規乱数を発生させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
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
#set.seed(8437)
##データの設定
hh <- 1500   #サンプル数
member <- 10   #選択可能メンバー数
k <- 3   #3位まで選択

##説明変数の発生
cont <- 2
X1.cont <- matrix(rnorm(hh*member*2, 0, 1), nrow=hh, ncol=2*member)

bin <- 2
X1.bin <- matrix(0, nrow=hh, ncol=2*member)
for(i in 1:(2*member)){
  X1.bin[, i]  <- rbinom(hh, 1, runif(1, 0.35, 0.6))
}

##説明変数をベクトル形式のフォーマットに変更
#IDを設定
id <- rep(1:hh, rep(member-1, hh))
m <- rep(1:(member-1), hh)
ID <- data.frame(no=1:length(id), id=id, m=m)

#切片の設定
p <- c(1, rep(0, member-1))
Pop <- matrix(p, nrow=hh*member, ncol=member-1, byrow=T)
POP <- subset(Pop, rowSums(Pop) > 0)

#説明変数の設定
#相対効用に変更
X1r.cont1 <- X1.cont[, 1:(member-1)] - X1.cont[, member]
X1r.cont2 <- X1.cont[, (member+1):(2*member-1)] - X1.cont[, (2*member)]
X1r.bin1 <- X1.bin[, 1:(member-1)] - X1.bin[, member]
X1r.bin2 <- X1.bin[, (member+1):(2*member-1)] - X1.bin[, (2*member)]

#ベクトル形式に変更
X1v.cont1 <- as.numeric(t(X1r.cont1))
X1v.cont2 <- as.numeric(t(X1r.cont2))
X1v.bin1 <- as.numeric(t(X1r.bin1))
X1v.bin2 <- as.numeric(t(X1r.bin2))

##データを結合
X <- data.frame(pop=POP, c1=X1v.cont1, c2=X1v.cont2, b1=X1v.bin1, b2=X1v.bin2)
XM <- as.matrix(X)


####応答変数を発生####
##ランクデータを発生させる
##分散共分散行列の設定
corM <- corrM(col=member-1, lower=-0.6, upper=0.70)   #相関行列を作成
Sigma <- covmatrix(col=member-1, corM=corM, lower=1, upper=1)   #分散共分散行列
Cov <- Sigma$covariance

##妥当なランキングが発生するまで繰り返す
for(i in 1:10000){
  print(i)
  
  #回帰係数のパラメ-タの設定
  a0 <- runif(member-1, 0.2, 3.1)
  a1 <- runif(cont, 0, 1.1)
  a2 <- runif(bin, -0.9, 1.2)
  a <- c(a0, a1, a2)
  alpha.t <- a
  
  ##相対効用を発生させる
  err1 <- mvrnorm(hh, rep(0, member-1), Cov)   #誤差構造
  U1.mean <- matrix(XM %*% a, nrow=hh, ncol=member-1, byrow=T)   #相対効用の平均構造
  U1 <- U1.mean + err1   #誤差構造を加えた構造
  
  ##効用最大化原理に基づき相対順位を決定
  Rank.full <- t(apply(cbind(U1, 0), 1, function(x) order(x, decreasing=TRUE)))
  Rank <- Rank.full[, 1:3]
  
  #基準メンバーが適当な人数に選ばれるまでループさせる
  if(sum(Rank.full[, 1:k]==member) > 30 & sum(Rank.full[, 1:k]==member) < 150) {break} else {next}
}

#発生させたデータの確認
Rank.full
apply(Rank.full, 2, table) #順位ごとの集計


##多項データを発生
##パラメータの設定
##妥当な応答変数が発生するまで繰り返す
for(i in 1:10000){
  print(i)
  
  #回帰係数のパラメ-タの設定
  b <- a + runif(length(a), -0.35, 0.35)
  beta.t <- b
  
  ##相対効用を発生させる
  err2 <- mvrnorm(hh, rep(0, member-1), Cov)   #誤差構造
  U2.mean <- matrix(XM %*% b, nrow=hh, ncol=member-1, byrow=T)   #相対効用の平均構造
  U2 <- U2.mean + err2   #誤差構造を加えた構造
  
  ##効用最大化原理に基づき選択メンバーを決定
  y <- apply(cbind(U2, 0), 1, which.max)
  
  #基準メンバーが適当な人数に選ばれるまでループさせる
  if(sum(y==member) > 20 & sum(y==member) < 100) {break} else {next}
}

#選択メンバーを0、1行列に変換
Y <- matrix(0, nrow=hh, ncol=member)
for(i in 1:hh){
  Y[i, y[i]] <- 1
}

table(y)   #選択メンバーの単純集計


####マルコフ連鎖モンテカルロ法で多項-ランクプロビットモデルを推定####
####MCMC推定のための推定準備####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##多変量正規分布の条件付き期待値と分散を計算する関数
cdMVN <- function(mu, Cov, dependent, U){
  
  #分散共分散行列のブロック行列を定義
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent]
  Cov21 <- Cov[-dependent, dependent]
  Cov22 <- Cov[-dependent, -dependent]
  
  #条件付き分散と条件付き平均を計算
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mu[, dependent] + t(CDinv %*% t(U[, -dependent] - mu[, -dependent]))   #条件付き平均を計算
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #条件付き分散を計算
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- array(0, dim=c(R/keep))   #対数尤度の保存用

##データの設定
#説明変数の多次元配列化
X.array <- array(0, dim=c(member-1, ncol(XM), hh))
for(i in 1:hh){
  X.array[, , i] <- XM[ID$id==i, ]
}
YX.array <- array(0, dim=c(member-1, ncol(XM)+1, hh))

#IDの設定
id_r <- matrix(1:nrow(XM), nrow=hh, ncol=member-1, byrow=T)

##推定プロセスの格納用配列
UM1 <- matrix(0, nrow=hh, ncol=member-1)
UM2 <- matrix(0, nrow=hh, ncol=member-1)
util.M1 <- matrix(0, nrow=hh, ncol=member-1)
util.M2 <- matrix(0, nrow=hh, ncol=member-1)

##事前分布の設定
nu <- member   #逆ウィシャート分布の自由度
V <- solve(0.1*diag(member-1))    #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(XM))  #回帰係数の平均の事前分布
Adelta <- solve(100 * diag(rep(1, ncol(XM)))) #回帰係数の事前分布の分散

##サンプリング結果の保存用配列
Util1 <- array(0, dim=c(hh, member-1, R/keep))
Util2 <- array(0, dim=c(hh, member-1, R/keep))
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
THETA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
SIGMA1 <- matrix(0, nrow=R/keep, ncol=(member-1)^2)
SIGMA2 <- matrix(0, nrow=R/keep, ncol=(member-1)^2)


##初期値の設定
oldalpha <- c(runif(member-1, 0, 3), runif(ncol(XM)-(member-1), 0, 1.0))   

#分散共分散行列の初期値
corM.f <- corrM(col=member-1, lower=0, upper=0)   #相関行列を作成
Sigma.f <- covmatrix(col=member-1, corM=corM.f, lower=1, upper=1)   #分散共分散行列
oldcov <- Sigma.f$covariance

#効用の平均構造の初期値
old.utilm1 <- matrix(XM %*% oldalpha, nrow=hh, ncol=member-1, byrow=T)

#効用の初期値
old.util1 <- old.utilm1 + mvrnorm(nrow(old.utilm1), rep(0, member-1), oldcov)


####マルコフ連鎖モンテカルロ法で多項-ランクプロビットモデルを推定####
