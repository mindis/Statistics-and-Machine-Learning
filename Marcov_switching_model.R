#####マルコフ切り替えモデル#####
library(MASS)
library(MSwM) 
library(reshape2)
library(gtools)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
n <- 300   #サンプル数
k1 <- 4   #切り替え数
k2 <- 6   #観測確率のパラメータ数


##初期確率の定義
Pf <- c(0.4, 0.1, 0.3, 0.2)

##推移行列の定義
pr1 <- c(0.2, 0.1, 0.4, 0.3)
pr2 <- c(0.1, 0.5, 0.2, 0.2)
pr3 <- c(0.3, 0.1, 0.3, 0.3)
pr4 <- c(0.1, 0.2, 0.3, 0.4)
Pr <- rbind(pr1, pr2, pr3, pr4)

##観測確率の定義
P <- matrix(0, nrow=k1, ncol=k2)
for(i in 1:k1){
  alpha <- runif(1, 0.4, 1)
  P[i, ] <- rdirichlet(1, rep(alpha, k2))
}


##応答変数の発生
Z <- matrix(0, nrow=n, ncol=k1)
Y <- matrix(0, nrow=n, ncol=k2)

#潜在変数の初期値
Z[1, ] <- rmultinom(1, 1, Pf)
Y[1, ] <- rmultinom(1, 1, P[which.max(Z[1, ]), ])

#2回目以降の応答変数を逐次的に発生させる
for(i in 2:n){
  Z[i, ] <- rmultinom(1, 1, Pr[which.max(Z[i-1, ]), ])
  Y[i, ] <- rmultinom(1, 1, P[which.max(Z[i, ]), ])
}


####EMアルゴリズムでマルコフ切り替えモデルを推定####
##初期値の設定
#初期確率の設定
rho <- rep(0.25, k1)   

#マルコフ推移行列の初期値の設定
A <- matrix(0, nrow=k1, ncol=k1)
for(i in 1:k1){
  p_rand <- runif(k1, 0.1, 1)
  A[i, ] <- p_rand / sum(p_rand)
}

#観測モデルのパラメータ
B <- matrix(0, nrow=k1, ncol=k2)
for(i in 1:k1){
  p_rand <- runif(k2, 0.1, 1)
  B[i, ] <- p_rand / sum(p_rand)
}
y <- Y %*% 1:k2

####EMアルゴリズム####
##前向きアルゴリズムでalphaを推定
alpha <- matrix(0, nrow=n, ncol=k1)
alpha1 <- matrix(0, nrow=n, ncol=k1)
alpha_s <- matrix(0, nrow=n, ncol=k1)
alpha_mu <- rep(0, n)
B_vec <- matrix(0, nrow=n, ncol=k1)

B_vec[1, ] <- B[, y[1]]
alpha[1, ] <- rho * B_vec[1, ]
alpha1[1, ] <- alpha[1, ]
alpha_mu[1] <- 1 / sum(alpha[1, ])
alpha_s[1, ] <- alpha_mu[1] * alpha[1, ] 

for(i in 2:nrow(Y)){
 B_vec[i, ] <- B[, y[i]]
 alpha[i, ] <- (matrix(alpha_s[i-1, ], nrow=1, ncol=k1) %*% A) * B_vec[i, ]
 alpha1[i, ] <- alpha1[i-1, ] %*% A * B_vec[i, ]
 alpha_mu[i] <- 1 / sum(alpha[i, ])
 alpha_s[i, ] <- alpha[i, ] * alpha_mu[i]
}

##後ろ向きアルゴリズムでbetaを推定
beta <- matrix(0, nrow=n, ncol=k1)
beta1 <- matrix(0, nrow=n, ncol=k1)
beta_s <- matrix(0, nrow=n, ncol=k1)

beta[n, ] <- 1
beta1[n, ] <- 1
beta_s[n, ] <- alpha_mu[n]

for(i in n:2){
  beta[i-1, ] <- A %*% (B_vec[i, ] * beta_s[i, ])
  beta1[i-1, ] <- A %*% (B_vec[i, ] * beta1[i, ])
  beta_s[i-1, ] <- beta[i-1, ] * alpha_mu[i-1, ] 
}

##パラメータを更新
A_vec <- matrix(A[1, ], nrow=n-1, ncol=k1, byrow=T)

a11 <- alpha1[1:(n-1), ] * A_vec * B_vec[2:n, ] * beta1[2:n, ]
a12 <- alpha1[1:(n-1), ] * beta1[1:(n-1), ] 

a21 <- alpha_s[1:(n-1), ] * A_vec * B_vec[2:n, ] * beta_s[2:n, ]
a22 <- alpha_s[1:(n-1), ] * beta_s[1:(n-1), ] / alpha_mu[1:(n-1), ]

sum(colSums(a11)/colSums(a12))
sum(colSums(a21)/colSums(a22))


P2 <- function(x, A, B, rho) {
  # forward algorithm
  N <- length(x)  # size of data
  alpha <- rho * B[, x[1]]  # alpha1(i) = rho_i * b(i, x_1)
  for (n in 2:N)
    alpha <- (matrix(alpha, nrow=1, ncol=length(alpha)) %*% A) * B[, x[n]]
  return(alpha)
}

P3 <- function(x, A, B, rho) {
  # backward algorithm
  N <- length(x)  # size of data
  beta <- 1       # beta_n
  for (n in N:2) {
    # compute upto beta_1
    beta <- A %*% matrix(B[, x[n]] * beta, ncol=1, nrow=nrow(B))
  }
  AA <- rho * beta * B[, x[1]] 
  return(AA)
}
beta_s[1, ]
rho * beta1[2, ] *B[, Y[1, ]==1]/sum(rho * beta1[2, ] *B[, Y[1, ]==1])
rho * beta1[2, ] * B_vec[1, ] 
P3(Y %*% 1:k2, A, B, rho)
beta1[1, ]

P2(Y %*% 1:k2, A, B, rho)

