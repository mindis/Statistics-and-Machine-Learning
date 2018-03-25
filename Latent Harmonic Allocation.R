#####Latent Harmonic Allocation(Nested Mixture Gaussian Distribution Model)#####
options(warn=2)
library(MASS)
library(mclust)
library(flexmix)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(mvtnorm)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####多変量正規分布の乱数を発生させる関数を定義####
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

####データの生成####
##データの設定
k1 <- 5
k2 <- 7
d <- 3000   #データ数
w <- rpois(d, rgamma(d, 30, 0.7))
f <- sum(w)
v <- 5   #変数数

#IDの設定
d_id <- rep(1:d, w)
n_id <- c()
for(i in 1:d){
  n_id <- c(n_id, 1:w[i])
}

##パラメータの設定
#ディリクレ分布のパラメータ
alpha1 <- rep(0.4, k1)
alpha2 <- rep(0.5, k2)


#多変量正規分布のパラメータ
Cov <- Covt <- array(0, dim=c(v, v, k1))
Mu <- Mut <- matrix(0, nrow=k1, ncol=v)
for(j in 1:k1){
  corr <- corrM(v, -0.8, 0.9, 0.01, 0.2)
  Cov[, , j] <- Covt[, , j] <- covmatrix(v, corr, 1.0, 5.0)$covariance
  Mu[j, ] <- Mut[j, ] <- runif(v, 10.0, 30.0)
}

#固定パラメータ
beta0 <- betat0 <- mvrnorm(k2, rep(0, v), diag(9.0, v))

#ディリクレ分布のパラメータを生成
theta1 <- thetat1 <- extraDistr::rdirichlet(d, alpha1)
theta2 <- thetat2 <- extraDistr::rdirichlet(k1, alpha2)


##LHAモデルに基づきデータを生成
Data_list <- Z1_list <- Z2_list <- list()

for(i in 1:d){
  #基底を生成
  z1 <- rmnom(w[i], 1, theta1[i, ])
  z1_vec <- as.numeric(z1 %*% 1:k1)
  
  #固定パラメータの割当を生成
  z2 <- rmnom(w[i], 1, theta2[z1_vec, ])
  z2_vec <- as.numeric(z2 %*% 1:k2)
  
  #多変量正規分布から観測データを生成
  mu <- Mu[z1_vec, ] + beta0[z2_vec, ]   #平均パラメータ
  y <- matrix(0, nrow=w[i], ncol=v)
  for(j in 1:w[i]){
    y[j, ] <- mvrnorm(1, mu[j, ], Cov[, , z1_vec[j]])
  }
  
  #データを格納
  Data_list[[i]] <- y
  Z1_list[[i]] <- z1
  Z2_list[[i]] <- z2
}

#リストを変換
Data <- do.call(rbind, Data_list)
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)


####変分ベイズEMアルゴリズムでHLAを推定####
##多変量正規分布の尤度関数
dmv <- function(x, mean.vec, S, S_det, S_inv){
  LLo <- (2*pi)^(-nrow(S)/2) * S_det^(-1/2) *
    exp(-1/2 * (x - mean.vec) %*% S_inv %*% (x - mean.vec))
  return(LLo)
}

##事前分布の設定
#多変量正規分布の事前分布
mu0 <- rep(0, v)
sigma0 <- 100
sigma0_inv <- 1/sigma0

#逆ウィシャート分布の事前分布
nu <- v + 1
V <- nu * diag(v)
inv_V <- solve(V)

#ディリクレ分布の事前分布
alpha <- 1

##パラメータの真値
theta1 <- thetat1
theta2 <- thetat2
mu <- Mu
Cov <- Covt
beta0 <- betat0


##初期値の設定
theta1 <- rdirichlet(d, rep(10.0, k1))
theta2 <- rdirichlet(k1, rep(10.0, k2))
mu <- mvrnorm(k1, rep(20, v), diag(10^2, v))
Cov <- array(diag(3^2, v), dim=c(v, v, k1))
beta0 <- mvrnorm(k2, rep(0, v), diag(3^2, v))


#インデックスを作成
index_k1 <- matrix(1:(k1*k2), nrow=k1, ncol=k2, byrow=T)
index_column <- rep(1:k1, rep(k2, k1))
d_list <- d_vec <- list()
for(i in 1:d){
  d_list[[i]] <- which(d_id==i)
  d_vec[[i]] <- rep(1, length(d_list[[i]]))
}


####変分ベイズEMアルゴリズムでパラメータを推定####
for(rp in 1:1000){
  
  ##Eステップで潜在変数zを推定
  #多変量正規分布の対数尤度を計算
  Li <- matrix(0, nrow=f, ncol=k1*k2)
  for(i in 1:k1){
    for(j in 1:k2){
      Li[, index_k1[i, ][j]] <- dmvnorm(Data, mu[i, ] + beta0[j, ] , Cov[, , i], log=TRUE)
    }
  }
 
  #潜在変数zの事前分布の設定
  theta1_z <- log(theta1)[d_id, index_column]
  theta2_z <- matrix(as.numeric(t(log(theta2))), nrow=f, ncol=k1*k2, byrow=T)
  
  #潜在変数zの計算
  LLi <- theta1_z + theta2_z + Li   #潜在変数zの対数尤度
  z_par <- exp(LLi - rowMaxs(LLi))   #尤度に変換
  z_rate <- z_par / rowSums(z_par)   #潜在変数z
  
  #混合率を更新
  r <- colSums(z_rate) / f
  
  ##Mステップで多変量正規分布のパラメータと固定パラメータを更新
  #平均および分散共分散行列の変分事後平均を推定
  for(i in 1:k1){
    weighted_data <- matrix(0, nrow=f, ncol=v)
    weighted_er <- matrix(0, nrow=f, ncol=v)
    
    for(j in 1:k2){
      #重み付きパラメータを更新
      weighted_data <- weighted_data + z_rate[, index_k1[i, j]] * (Data - matrix(beta0[j, ], nrow=f, ncol=v, byrow=T))
       weighted_er <- weighted_er + (z_rate[, index_k1[i, j]] * Data) -
         (z_rate[, index_k1[i, j]] * matrix(mu[i, ] + beta0[j, ], nrow=f, ncol=v, byrow=T))
    }
    
    #変分事後平均を推定
    n <- sum(z_rate[, index_k1[i, ]])
    mu[i, ] <- colSums(weighted_data) / (f * sum(r[index_k1[i, ]]) + sigma0_inv)   #変分事後平均
    Cov[, , i] <- (inv_V + t(weighted_er) %*% weighted_er) / n   #変分事後分散共分散行列
  }
  
  #固定パラメータの変分事後平均を推定
  for(i in 1:k2){
    weighted_data <- matrix(0, nrow=f, ncol=v)
    for(j in 1:k1){
      #重み付きパラメータを更新
      weighted_data <- weighted_data + z_rate[, index_k1[j, i]] * (Data - matrix(mu[j, ], nrow=f, ncol=v, byrow=T))
    }
    beta0[i, ] <- colSums(weighted_data) / (f * sum(r[index_k1[, i]]) + sigma0_inv)   #変分事後平均を推定
  }
  
  ##潜在変数の割当確率を更新
  #基底の分布の変分事後平均を推定
  Zi_T <- t(z_rate)   #潜在変数zの転置行列
  wsum01 <- matrix(0, nrow=d, ncol=k1*k2)
  wsum1 <- matrix(0, nrow=d, ncol=k1)
  for(i in 1:d){
    wsum01[i, ] <- Zi_T[, d_list[[i]]] %*% d_vec[[i]]
  }
  for(j in 1:k1){
    wsum1[, j] <- rowSums(wsum01[, index_k1[j, ]]) + alpha   #ディリクレ分布のパラメータ
  }
  theta1 <- wsum1 / rowSums(wsum1)   #変分事後平均の推定
  
  
  #固定パラメータの割当分布の変分事後平均を推定
  for(j in 1:k1){
    wsum2 <- colSums(z_rate[, index_k1[j, ]]) + alpha   #ディリクレ分布のパラメータ
    theta2[j, ] <- wsum2 / sum(wsum2)   #変分事後平均の推定
  }
  
  print(sum(log(rowSums(exp(LLi)))))
}

