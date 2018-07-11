#####Latent variable bayesian hieralchical probit model#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)


####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
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
##データの設定
s <- 2   #応答変数数
dir <- 150   #ディレクトリ数
item <- 10000   #アイテム数
dir_freq <- rtpois(item, 1.0, a=0, b=5)   #アイテムごとのディレクトリ数
max_dir <- max(dir_freq)   #ディレクトリの最大数
w <- rpois(item, rgamma(item, 10.0, 0.075))   #アイテムあたりのサンプル数
f <- sum(w)   #総サンプル数

#IDの設定
item_id <- rep(1:item, w)
t_id <- as.numeric(unlist(tapply(1:f, item_id, rank)))

#ディレクトリの生成
dir_x <- matrix(0, nrow=item, ncol=dir)
dir_data <- matrix(0, nrow=item, ncol=max(dir_freq))
pr <- runif(dir, 0.1, 3.0)

for(i in 1:item){
  repeat {
    dir_x[i, ] <- rmnom(1, dir_freq[i], pr)
    if(sum(dir_x[i, ] <= 1)==dir) break
  }
  dir_data[i, 1:sum(dir_x[i, ])] <- (dir_x[i, ] * 1:dir)[dir_x[i, ]!=0]
}
dir_vec0 <- as.numeric(t(dir_x * matrix(1:dir, nrow=item, ncol=dir, byrow=T)))
dir_vec <- dir_vec0[dir_vec0!=0]
storage.mode(dir_data) <- "integer"


##説明変数の生成
#アイテムの説明変数
k1 <- 3; k2 <- 4; k3 <- 5
x1 <- matrix(runif(f*k1, 0, 1), nrow=f, ncol=k1)
x2 <- matrix(0, nrow=f, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(f, 1, pr)
}
x3 <- rmnom(f, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #データを結合
k <- ncol(x)

##パラメータを生成
#ディレクトリ生成確率の生成
lambda <- matrix(0, nrow=item, ncol=max(dir_freq))
for(i in 1:item){
  if(dir_freq[i]==1){
    lambda[i, 1] <- 1
  } else {
    lambda[i, 1:dir_freq[i]] <- extraDistr::rdirichlet(1, rep(2.5, dir_freq[i]))
  }
}
#ディレクトリ割当を生成
lambdat <- lambda
dir_z <- rmnom(f, 1, lambda[item_id, ])

#応答変数の誤差相関
Cov <- Covt <- matrix(c(1, 0.5, 0.5, 1.0), nrow=s, ncol=s)  

##応答変数が妥当な数値になるまで繰り返す
for(rp in 1:1000){
  
  #回帰パラメータを生成
  theta1 <- thetat1 <- runif(k, -1.3, 1.1); theta2 <- thetat2 <- runif(k, -1.4, 0.9)   #階層モデルの平均
  tau1 <-taut1 <- runif(k, 0.025, 0.3) * diag(k); tau2 <- taut2 <- runif(k, 0.025, 0.3) * diag(k)   #階層モデルの分散
  beta1 <- betat1 <- mvrnorm(dir, theta1, tau1); beta2 <- betat2 <- mvrnorm(dir, theta2, tau2)   #ディレクトリ別の回帰係数
  
  ##応答変数を生成
  #回帰モデルの平均構造
  z <- rowSums(dir_data[item_id, ] * dir_z)   #ディレクトリ割当
  mu1 <- as.numeric((x * beta1[z, ]) %*% rep(1, k))
  mu2 <- as.numeric((x * beta2[z, ]) %*% rep(1, k))
  mu <- cbind(mu1, mu2)
  
  #多変量正規分布から応答変数を生成
  er <- mvrnorm(f, rep(0, s), Cov)  
  U <- mu + er   #潜在効用を設定
  y <- ifelse(U > 0, 1, 0)   #応答絵変数を生成
  
  #ストップ判定
  if(max(colMeans(y)) < 0.4 & min(colMeans(y)) > 0.15) break
}

####マルコフ連鎖モンテカルロ法でLatent variable bayesian hieralchical probit modelを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##多変量正規分布の条件付き期待値と条件付き分散を計算する関数
cdMVN <- function(mean, Cov, dependent, U){
  
  #分散共分散行列のブロック行列を定義
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent, drop=FALSE]
  Cov21 <- Cov[-dependent, dependent, drop=FALSE]
  Cov22 <- Cov[-dependent, -dependent]

  
  #条件付き分散と条件付き平均を計算
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mean[, dependent] + t(CDinv %*% t(U[, -dependent] - mean[, -dependent]))   #条件付き平均を計算
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #条件付き分散を計算
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##MCMCの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##事前分布の設定
#モデルの事前分布
alpha0 <- 1.0   #ディレクトリ割当の事前分布
nu1 <- s   #逆ウィシャート分布の自由度
V1 <- nu1 * diag(s)   #逆ウィシャート分布の自由度

#階層モデルの事前分布
nu2 <- k   #逆ウィシャート分布の自由度
V2 <- nu * diag(k)   #逆ウィシャート分布の自由度
Deltabar <- rep(0, k)   #回帰係数の事前分布

##パラメータの真値
lambda <- lambdat
theta1 <- thetat1
theta2 <- thetat2
tau1 <- taut1
tau2 <- taut2
beta1 <- betat1
beta2 <- betat2

##切断性分布の切断領域を定義
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)

####ギブスサンプリングでパラメータをサンプリング####
##切断正規分布から潜在効用を生成
#効用の平均構造
mu1 <- as.numeric((x * beta1[dir_data[item_id, 1], ]) %*% rep(1, k))
mu2 <- as.numeric((x * beta2[dir_data[item_id, 1], ]) %*% rep(1, k))
mu <- cbind(mu1, mu2)

cdMVN(mu, Cov, 1, U)


