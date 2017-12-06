#####Twitter LDA#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(gtools)
library(bayesm)
library(extraDistr)
library(monomvn)
library(glmnet)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(54876)

####データの発生####
##文書データの設定
hh <- 1000
tweet <- rpois(hh, rgamma(hh, 150.0, 1.0))
d <- sum(tweet)
w <- rpois(d, 11.5)
v1 <- 150   #評価対象の語彙数
v2 <- 150   #一般語の語彙数
v <- v1+v2   #総語彙数
k <- 10   #トピック数

#IDの設定
u_id <- rep(1:hh, tweet)   
t_id <- c()
index_id <- list()
for(i in 1:hh){
  t_id <- c(t_id, 1:tweet[i])
  index_id[[i]] <- which(u_id==i)
}

##パラメータの設定
#ディクレリ事前分布を設定
alpha0 <- rep(0.3, k)   #ユーザー固有のディクレリ分布のパラメータ
alpha1 <- rep(0.4, v1)   #評価対象語のディクレリ分布のパラメータ
alpha2 <- rep(10, v2)   #一般語のディクレリ分布の事前分布のパラメータ
beta0 <- c(4.5, 3.5)   #一般語かどうかのベータ分布のパラメータ


#ディクレリ分布からパラメータを生成
thetat <- theta <- extraDistr::rdirichlet(hh, alpha0)   #ユーザートピックの生成
phit <- phi <- extraDistr::rdirichlet(k, alpha1)   #評価対象語の出現率の生成
lambdat <- lambda <- extraDistr::rdirichlet(1, alpha2)   #一般語の出現率の生成

##多項分布からトピックおよび単語データを生成
WX <- matrix(0, nrow=d, ncol=v)
Z0 <- list()
y0 <- list()
index_word1 <- 1:v1
index_word2 <- (v1+1):v

#tweetごとに1つのトピックを割り当て単語を生成
for(i in 1:hh){
  print(i)
  
  #tweetごとにトピックを生成
  z <- rmnom(tweet[i], 1, theta[i, ])
  z_vec <- z %*% 1:k
  index_hh <- index_id[[i]]
  
  #tweetに割り当てられたトピックから単語を生成
  for(j in 1:nrow(z)){
    
    #トピックに関係あるかどうかの潜在変数
    word <- w[index_hh[j]]
    par <- rbeta(word, beta0[1], beta0[2])
    y <- rbinom(word, 1, par)
    
    #潜在変数に基づいてtweetの単語を生成
    WX[index_hh[j], index_word1] <- rmnom(1, sum(y), phi[z_vec[j, ], ])   #評価対象語の生成
    WX[index_hh[j], index_word2] <- rmnom(1, sum(1-y), lambda)   #一般語の生成
    y0[[index_hh[j]]] <- y
  }
  Z0[[i]] <- as.numeric(z_vec)
}

#データ形式を変換
y_vec <- unlist(y0)
z_vec <- unlist(Z0)
storage.mode(WX) <- "integer"
gc(); gc()


####トピックモデル推定のためのデータの準備####
##データ推定用IDを作成
ID_list <- list()
td_list <- list()
wd_list <- list()

#IDごとにtweet_idおよび単語idを作成
for(i in 1:hh){
  print(i)
  
  #ユーザーIDを記録
  index_hh <- index_id[[i]]
  ID_list[[i]] <- rep(i, sum(w[index_id[[i]]]))
  td_list[[i]] <- rep(1:tweet[i], w[index_hh])
  
  for(j in 1:length(index_hh)){
    num1 <- WX[index_hh[j], ] * 1:v
    num2 <- which(num1 > 0)
    W1 <- WX[index_hh[j], (WX[index_hh[j], ] > 0)]
    wd_list[[index_hh[j]]] <- rep(num2, W1)
  }
}
#リストをベクトルに変換
ID_d <- unlist(ID_list)
td_d <- unlist(td_list)
wd <- unlist(wd_list)

##インデックスを作成
user_list <- list()
word_list <- list()
for(i in 1:length(unique(ID_d))){user_list[[i]] <- which(ID_d==i)}
for(i in 1:length(unique(wd))){word_list[[i]] <- which(wd==i)}

tweet_list <- list()
wsum <- cumsum(w)
for(i in 1:length(w)){
  if(i==1){
    tweet_list[[i]] <- 1:wsum[i]
  } else {
    tweet_list[[i]] <- (wsum[i-1]+1):(wsum[i])
  }
}


####マルコフ連鎖モンテカルロ法で対応トピックモデルを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #負担係数の格納用
  for(j in 1:k){
    #負担係数を計算
    Bi <- rep(theta[, j], w) * phi[j, c(wd)]   #尤度
    Bur[, j] <- Bi   
  }
  
  Br <- Bur / rowSums(Bur)   #負担率の計算
  r <- colSums(Br) / sum(Br)   #混合率の計算
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}



##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep

##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- rep(0.3, k)   
beta01 <- rep(0.4, v1)   
gamma01 <- rep(1, v2)   

##パラメータの初期値
theta <- extraDistr::rdirichlet(hh, rep(1, k))   #ユーザートピックの初期値
phi <- extraDistr::rdirichlet(k, rep(0.5, v))   #評価対象語の出現確率の初期値
lambda <- extraDistr::rdirichlet(k, rep(10, v))   #一般語の出現確率の初期値
r <- c(0.5, 0.5)   #評価対象語と一般語の割当率の混合率

##パラメータの格納用配列
THETA <- array(0, dim=c(hh, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
LAMBDA <- array(0, dim=c(k, v, R/keep))
Zi <- matrix(0, nrow=d, ncol=k)
storage.mode(Zi) <- "integer"


####ギブスサンプリングでパラメータをサンプリング####
phi


