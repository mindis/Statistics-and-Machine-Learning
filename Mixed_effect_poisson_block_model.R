#####mixed effect poisson block model#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(flexmix)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(506832)

####データの発生####
##データの設定
d <- 150   #アイテム数
k <- 7   #潜在変数数
N <- d*(d-1)/2   #総サンプル数
vec <- rep(1, k)

##IDを設定
id1 <- id2 <- c()
for(i in 1:(d-1)){
  id1 <- c(id1, rep(i, length((i+1):d)))
  id2 <- c(id2, (i+1):d)
}

##潜在変数の生成
#ディリクレ分布からパラメータを生成
alpha0 <- rep(0.1, k)
theta <- thetat <- extraDistr::rdirichlet(d, alpha0)
Z1 <- rmnom(N, 1, theta[id1, ])
Z2 <- rmnom(N, 1, theta[id2, ])
z1_vec <- as.numeric(Z1 %*% 1:k); z2_vec <- as.numeric(Z2 %*% 1:k)


##応答変数の生成
#パラメータを生成
cov <- covt <- 0.75
beta <- betat <- 0.8
alphat <- alpha <- rnorm(d, beta, cov)   #変量効果のパラメータ
phi0 <- matrix(rnorm(k*k, 0, 0.85), nrow=k, ncol=k)   #潜在変数のパラメータ
phi0[upper.tri(phi0)] <- 0
phi <- phi0 + t(phi0)
diag(phi) <- diag(phi0)
phit <- phi

#ポアソン分布の平均構造
mu <- alpha[id1] + alpha[id2] + (phi[z1_vec, ] * Z2) %*% vec
lambda <- exp(mu)


#ポアソン分布から応答変数を生成
y <- rpois(N, lambda)
sum(y); mean(y)
hist(y, xlab="頻度", main="アイテム間の出現頻度", col="grey", breaks=25)


####マルコフ連鎖モンテカルロ法でmixed effect poisson block modelを推定####
##変量効果ポアソン回帰モデルの対数尤度
loglike <- function(alpha, theta, y, y_factorial, z1, Z2, vec, id1, id2){
  
  #尤度を定義する
  lambda <- exp(alpha[id1] + alpha[id2] + (phi[z1, ] * Z2) %*% vec)   #平均構造
  LLi <- as.numeric(y*log(lambda)-lambda - lfactorial(y))   #対数尤度
  LL <- sum(LLi)   #対数尤度の和
  
  #結果を返す
  LL_value <- list(LLi=LLi, LL=LL)
  return(LL_value)
}


##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##事前分布の設定
#回帰係数の事前分布
alpha01 <- 0
beta01 <- 0
tau01 <- 0.01

#変量効果の事前分布
alpha02 <- 0
s02 <- 0.01
v02 <- 0.01

##パラメータの真値
theta <- thetat
alpha <- alphat
beta <- betat
cov <- covt
Zi1 <- Z1; Zi2 <- Z2
z_vec1 <- as.numeric(Zi1 %*% 1:k)
z_vec2 <- as.numeric(Zi2 %*% 1:k)

##初期値の設定
#変量効果の初期値
cov <- 0.5   #階層モデルの標準偏差の初期値
beta <- 0.8   #階層モデルの平均の初期値
mu <- rep(0, d)
for(i in 1:d){
  mu[i] <- mean(y[which(id1==i | id2==i)])
}
x <- rnorm(d, 0, 0.5)
rank_mu <- ceiling(rank(mu))
alpha <- sort(x, decreasing=TRUE)[rank_mu]   #変量効果の初期値

#潜在変数のパラメータ
phi0 <- matrix(rnorm(k*k, 0, 0.3), nrow=k, ncol=k)   #潜在変数のパラメータ
phi0[upper.tri(phi0)] <- 0
phi <- phi0 + t(phi0)

#潜在変数の初期値
theta <- extraDistr::rdirichlet(d, rep(2.0, k))
Zi1 <- rmnom(N, 1, theta[id1, ])
Zi2 <- rmnom(N, 1, theta[id2, ])
z1_vec <- as.numeric(Zi1 %*% 1:k); z2_vec <- as.numeric(Zi2 %*% 1:k)


##定数を計算
y_factorial <- lfactorial(y)
Y_factorial <- matrix(y_factorial, nrow=N, ncol=k*k)
upper_tri <- matrix(as.logical(upper.tri(phi) + diag(1, k)), nrow=k, ncol=k)
vec <- rep(1, k)

##インデックスを設定
item_list <- list()
seg_list1 <- seg_list2 <- list()
for(i in 1:d){
  item_list[[i]] <- which(id1==i | id2==i)
  seg_list1[[i]] <- as.numeric(id1[item_list[[i]]]!=i)
  seg_list2[[i]] <- as.numeric(id2[item_list[[i]]]!=i)
}


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
##メトロポリスヘイスティング法で変量効果をサンプリング
#新しいパラメータをサンプリング
alphad <- alpha
alphan <- alphad + rnorm(d, 0, 0.2)

#事前分布の誤差
er_new <- alphan - beta 
er_old <- alphad - beta

#対数尤度と対数事前分布を設定
lognew0 <- loglike(alpha, phi, y, y_factorial, z1_vec, Zi2, vec, id1, id2)$LLi
logold0 <- loglike(alpha, phi, y, y_factorial, z1_vec, Zi2, vec, id1, id2)$LLi
logpnew <- -0.5 * (er_new^2 / cov)
logpold <- -0.5 * (er_old^2 / cov)
 
#アイテムごとに対数尤度の和を取る
lognew <- logold <- rep(0, d)
for(i in 1:d){
  lognew[i] <- sum(lognew0[item_list[[i]]])
  logold[i] <- sum(logold0[item_list[[i]]])
}

#MHサンプリング
rand <- runif(d)   #一様分布から乱数を発生
LLind_diff <- exp(lognew + logpnew - logold - logpold)   #採択率を計算
alpha2 <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff

#alphaの値に基づき新しいbetaを採択するかどうかを決定
flag <- (alpha2 >= rand)*1 + (alpha2 < rand)*0
alpha <- flag*alphan + (1-flag)*alphad   #alphaがrandを上回っていたら採択


##ギブスサンプリングで潜在変数をサンプリング
#多項分布から潜在変数をサンプリング
seg_list[[1]]

for(i in 1:d){
  i <- 1
  #インデックスを抽出
  index1 <- item_list[[i]]
  z1_allocation <- seg_list1[[i]]
  z2_allocation <- seg_list2[[i]]
  
  #ポアソン分布の平均構造
  Zi_pairs <- Zi1[index1, ] * z1_allocation + Zi2[index1, ] * z2_allocation   #対となる潜在変数
  Zi_mu <- phi[as.numeric(Zi_pairs %*% 1:k), ]
  lambda <- exp(matrix(alpha[id1[index1]] + alpha[id2[index1]], nrow=length(index1), ncol=k) + Zi_mu)
  
  #対数尤度と潜在変数の割当確率
  LLi <- y[index1]*log(lambda)-lambda - y_factorial[index1]   #対数尤度
  z_par <- exp(LLi - rowMaxs(LLi)) * matrix(theta[i, ], nrow=length(index1), ncol=k, byrow=T)
  z_rate <- z_par / rowSums(z_par)
  
  #多項分布より潜在変数をサンプリング
  Zi[i, ] <- rmnom(1, 1, z_rate)
  id1
  id2
}
colSums(Zi)
colSums(Z)        

