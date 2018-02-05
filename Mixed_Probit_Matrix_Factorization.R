#####Mixed Probit Matrix Factorization#####
library(MASS)
library(matrixStats)
library(FAdist)
library(mnormt)
library(NMF)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

####多変量正規乱数を発生させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
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
  D <- ifelse(val < 0, val + abs(val) + 0.1, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}


####データの発生####
hh <- 3000
item <- 300
k <- 10   #潜在変数数

##モデルの仮定に従いデータを生成
#相関行列を設定
Cor <- corrM(k, -0.7, 1.0)

#ユーザーアイテム特徴行列のパラメータを生成
W <- WT <- mvrnorm(hh, rep(0, k), Cor)
H <- HT <- matrix(0, nrow=k, ncol=item)
for(j in 1:k){
  par <- rbeta(1, 20.0, 20.0)
  H[j, ] <- HT[j, ] <- rbinom(item, 1, par)
}

#ユーザーとアイテムの変量効果パラメータを生成
alpha <- alphat <- rnorm(hh, -1.25, 1)
beta <- betat <- rnorm(item, -1.25, 1)


#プロビットモデルからアイテム購買行列を生成
alpha_matrix <- matrix(alpha, nrow=hh, ncol=item)
beta_matrix <- matrix(beta, nrow=hh, ncol=item, byrow=T)
Util <- alpha_matrix + beta_matrix + W %*% H   #効用関数
y0 <- rnorm(hh*item, as.numeric(Util), 1)   #正規分布から購買有無を生成
y <- ifelse(y0 > 0, 1, 0)
Data <- matrix(y, nrow=hh, ncol=item)
storage.mode(Data) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")   #スパース行列化


####マルコフ連鎖モンテカルロ法でMixed Probit Matrix Factorizationを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b, hh, item){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(matrix(runif(length(mu)), nrow=hh, ncol=item)*(FB-FA)+FA, mu, sigma))
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


##アルゴリズムの設定
R <- 10000
keep <- 2
disp <- 10
iter <- 0

##事前分布の設定
#変量効果の分散の事前分布
alpha01 <- 0.01
beta01 <- 0.01

#行列分解のパラメータの事前分布
nu <- k   #逆ウィシャーと分布の自由度
V <- nu * diag(rep(1, k))   #逆ウィシャート分布のパラメータ

##初期値の設定
#変量効果の初期値
user_rnorm <- rnorm(hh, 0, 1)
item_rnorm <- rnorm(item, 0, 1)
user_sort <- order(user_rnorm)
item_sort <- order(item_rnorm)
alpha <- user_rnorm[user_sort][ceiling(rank(rowSums(Data)))]
beta <- item_rnorm[item_sort][ceiling(rank(colSums(Data)))]

#行列分解のパラメータの初期値
W <- mvrnorm(hh, rep(0, k), diag(k))
H <- matrix(0, nrow=k, ncol=item)
for(j in 1:k){
  par <- rbeta(1, 30.0, 30.0)
  H[j, ] <- rbinom(item, 1, par)
}

##サンプリング結果の保存用配列
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, item, R/keep))
ALPHA <- matrix(0, nrow=R/keep, ncol=hh)
BETA <- matrix(0, nrow=R/keep, ncol=item)

##切断領域を定義
a <- ifelse(Data==0, -100, 0)
b <- ifelse(Data==1, 100, 0)
a_vec <- as.numeric(a)
b_vec <- as.numeric(b)


#パラメータの真値
W <- WT
H <- HT
r <- rowMeans(H)

alpha <- alphat
alpha_matrix <- matrix(alpha, nrow=hh, ncol=item)
beta <- betat
beta_matrix <- matrix(beta, nrow=hh, ncol=item, byrow=T)
util_mu <- alpha_matrix + beta_matrix + W %*% H   #効用関数
cov_inv <- solve(Cor)

####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####

##切断正規分布よりアイテム購買の効用値をサンプリング
util_mu <- alpha_matrix + beta_matrix + W %*% H   #効用関数
util <- rtnorm(util_mu, 1, a, b, hh, item)

##ユーザー特徴行列をサンプリング
#アイテム特徴行列の誤差を設定
util_fearture <- util - alpha_matrix - beta_matrix

#アイテム特徴行列の事後分布のパラメータ
XX <- H %*% t(H)
XXV <- solve(XX + cov_inv)

for(i in 1:hh){
  XXb <- H %*% util_fearture[i, ]
  beta_mean <- beta_mean <- XXV %*% XXb
  
  #多変量正規分布からbetaをサンプリング
  W[i, ] <- mnormt::rmnorm(1, beta_mean, XXV)
}

##アイテム特徴行列をサンプリング
const <- -1/2*log(2*pi)   #標準正規分布の対数尤度の定数

for(j in 1:k){
  #変数パターンを設定
  H1 <- H0 <- H
  H1[j, ] <- 1
  H0[j, ] <- 0
  
  #パターンごとの対数尤度を計算
  WH0 <- W %*% H0
  WH1 <- W %*% H1
  LL_item0 <- colSums(const -1/2*(util_fearture - WH0)^2)   
  LL_item1 <- colSums(const -1/2*(util_fearture - WH1)^2)
  LLi0 <- cbind(LL_item0, LL_item1)
  LLi <- exp(LLi0 - rowMaxs(LLi0)) * matrix(c(1-r[j], r[j]), nrow=item, ncol=2, byrow=T)   #尤度に変換
  
  #潜在変数の割当確率からHをサンプリング
  z_rate <- LLi[, 2] / rowSums(LLi)  #潜在変数の割当確率 
  H[j, ] <- rbinom(item, 1, z_rate)
  
  #混合率を更新
  r[j] <- mean(H[j, ])
}

##個体間変量効果をサンプリング
util_user <-  util - beta_matrix - W %*% H

#ユーザーごとの平均を推定
mu <- rowMeans(util_user)



