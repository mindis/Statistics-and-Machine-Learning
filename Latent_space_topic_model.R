#####連続空間トピックモデル#####
options(warn=2)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(5723)

####データの発生####
k <- 10   #トピック数
d <- 2000   #ユーザー数
v <- 500   #語彙数
w <- rpois(d, rgamma(d, 45, 0.3))   #1人あたりのページ閲覧数
f <- sum(w)   #総語彙数

#IDを設定
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}


##パラメータの設定
G0 <- GT0 <- extraDistr::rdirichlet(1, rep(2.0, v))   #単語出現率
u <- ut <- mvrnorm(d, rep(0, k), diag(k))
phi <- phit <- mvrnorm(v, rep(0, k), diag(k))

##データの生成
word_list <- list()
word_vec_list <- list()
WX <- matrix(0, nrow=d, ncol=v)

for(i in 1:d){
  #ディクレリ-多項分布から単語を生成
  alpha <- G0 * exp(u[i, ] %*% t(phi))
  words <- extraDistr::rdirmnom(w[i], 1, alpha)
  words_vec <- as.numeric(words %*% 1:v)
  
  #生成した単語を格納
  WX[i, ] <- colSums(words)
  word_list[[i]] <- words
  word_vec_list[[i]] <- words_vec
}

#リストを変換
word_vec <- unlist(word_vec_list)
word_data <- do.call(rbind, word_list)
sparse_data <- as(word_data, "CsparseMatrix")
storage.mode(WX) <- "integer"
rm(word_data); rm(word_list)


##インデックスを作成
dw_list <- list()
for(j in 1:v){
  index <- which(word_vec==j)
  dw_list[[j]] <- d_id[index]
}


####マルコフ連鎖モンテカルロ法で連続空間トピックモデルを推定####
##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##事前分布の設定
sigma <- diag(k)

##データの設定
X <- diag(v)
storage.mode(X) <- "iteger"

##パラメータの真値
G0 <- GT0
u <- ut
phi <- phit
inv_cov <- solve(diag(k))


##MH法でuをサンプリング
#新しいパラメータをサンプリング
u_old <- u
u_new <- u_old + mvrnorm(d, rep(0, k), diag(0.01, k))
alphan <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u_new %*% t(phi))
alphad <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u_old %*% t(phi))

#対数尤度と対数事前分布を計算
lognew1 <- extraDistr::ddirmnom(WX, w, alphan, log=TRUE)
logold1 <- extraDistr::ddirmnom(WX, w, alphad, log=TRUE)
logpnew1 <- -0.5 * rowSums(u_new %*% inv_cov * u_new)
logpold1 <- -0.5 * rowSums(u_old %*% inv_cov * u_old)

#MHサンプリング
rand <- runif(d)   #一様分布から乱数を発生
LLind_diff <- exp(lognew1 + logpnew1 - logold1 - logpold1)   #採択率を計算
tau <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff

#tauの値に基づき新しいbetaを採択するかどうかを決定
flag <- matrix(((tau >= rand)*1 + (tau < rand)*0), nrow=d, ncol=k)
u <- flag*u_new + (1-flag)*u_old   #alphaがrandを上回っていたら採択


##MH法でphiをサンプリング
#新しいパラメータをサンプリング
phid <- phi
phin <- phid + mvrnorm(v, rep(0, k), diag(0.02, k))
alphan <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u %*% t(phin))
alphad <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u %*% t(phid))

#単語ごとに候補分布を推定
lognew2 <- logold2 <- rep(0, v)

for(j in 1:v){
  #対数尤度を設定
  n <- length(dw_list[[j]])
  x <- matrix(X[j, ], nrow=n, ncol=v, byrow=T)
  par_new <- alphan[dw_list[[j]], ]
  par_old <- alphad[dw_list[[j]], ]
  lognew2[j] <- sum(extraDistr::ddirmnom(x, 1, par_new, log=TRUE))
  logold2[j] <- sum(extraDistr::ddirmnom(x, 1, par_old, log=TRUE))
}

#対数事前分布を設定
logpnew2 <- -0.5 * rowSums(phin %*% inv_cov * phin)
logpold2 <- -0.5 * rowSums(phid %*% inv_cov * phid)

