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
tweet <- rpois(hh, rgamma(hh, 140.0, 1.0))
d <- sum(tweet)
w <- rpois(d, 9.5)
f <- sum(w)
v1 <- 150   #評価対象の語彙数
v2 <- 100   #一般語の語彙数
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
phi0 <- cbind(extraDistr::rdirichlet(k, alpha1), matrix(runif(v2*k, 10^-15, 10^-10), nrow=k, ncol=v2)) 
phit <- phi <- phi0 / rowSums(phi0)   #評価対象語の出現率の生成
lambda0 <- c(runif(v1, 10^-15, 10^-10), extraDistr::rdirichlet(1, alpha2))
lambdat <- lambda <- lambda0 / sum(lambda0)


##多項分布からトピックおよび単語データを生成
WX <- matrix(0, nrow=d, ncol=v)
Z0 <- list()
y0 <- list()
index_word1 <- 1:v1
index_word2 <- (v1+1):v

#tweetごとに1つのトピックを割り当て単語を生成
for(i in 1:hh){
  
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
    WX[index_hh[j], ] <- rmnom(1, sum(y), phi[z_vec[j, ], ])   #評価対象語の生成
    WX[index_hh[j], ] <- WX[index_hh[j], ] + rmnom(1, sum(1-y), lambda)   #一般語の生成
    y0[[index_hh[j]]] <- y
  }
  Z0[[i]] <- as.numeric(z_vec)
}
sum(WX)
sum(w)


#データ形式を変換
y_vec <- unlist(y0)
z_vec <- unlist(Z0)
storage.mode(WX) <- "integer"
WX_sparse <- as(WX, "CsparseMatrix")
gc(); gc()


####トピックモデル推定のためのデータの準備####
##データ推定用IDを作成
ID_list <- list()
td_list <- list()
wd_list <- list()

#IDごとにtweet_idおよび単語idを作成
for(i in 1:hh){
  
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
burden_fr <- function(theta, phi, Data, id){
  #尤度を計算
  Bur <- theta[id, ] * exp(Data %*% t(log(phi)))
 
  #負担率を計算
  Br <- Bur / rowSums(Bur)
  r <- colSums(Br) / sum(Br)
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

wd[-index_zeros]

ID_d[index_zeros]
wd[index_zeros]
(count(c(ID_d[index_zeros], 1:hh))-1)


wd[index_zeros]
i <- c(1, 1, 1, 2, 2)
j <- c(1, 3, 5, 4, 2)
x <- c(1, 1, 1, 1, 1)
sparseMatrix(i, j, k, dim=c(5, 5))


##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep

##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- 1  
beta01 <- 0.05  
gamma01 <- 0.1 

##パラメータの初期値
theta <- extraDistr::rdirichlet(hh, rep(1, k))   #ユーザートピックの初期値
phi <- extraDistr::rdirichlet(k, rep(0.5, v))   #評価対象語の出現確率の初期値
lambda <- extraDistr::rdirichlet(k, rep(10, v))   #一般語の出現確率の初期値
y <- rbinom(f, 1, 0.5)
r <- c(0.5, 0.5)   #評価対象語と一般語の割当率の混合率

##パラメータの格納用配列
THETA <- array(0, dim=c(hh, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
LAMBDA <- array(0, dim=c(k, v, R/keep))
Zi <- matrix(0, nrow=d, ncol=k)
storage.mode(Zi) <- "integer"
WX[, 151:ncol(WX)] <- 0


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##tweetトピックをサンプリング
  #tweetごとにトピックの出現確率を計算
  Br <- burden_fr(theta, phi, WX, u_id)
  
  #多項分布から単語トピックをサンプリング
  tweet_rate <- Br$Br
  Z <- rmnom(d, 1, tweet_rate)
  z <- Z %*% 1:k
  
  ##単語がトピックと関係あるかどうかをサンプリング
  #潜在変数のパラメータ
  LH0 <- lambda[wd]
  
  Zi <- matrix(0, nrow=f, ncol=k)
  LH1 <- rep(0, f)
  zw <- rep(z, w)
  pf <- t(phi)[wd, ]
  
  for(j in 1:k){
    index <- which(zw==j)
    Zi[index, j] <- 1
    LH1[index] <- pf[index, j]
  }
  LH <- cbind(r[1]*LH1, r[2]*LH0)
  
  #潜在変数の割当確率から潜在変数を生成
  y_rate <- LH[, 1] / rowSums(LH)
  y <- rbinom(length(y_rate), 1, y_rate)
  cbind(y_rate, wd)
  
  
  #混合率を更新
  r0 <- mean(y)
  r <- c(r0, 1-r0)
  
  ##ディクレリ分布からトピック分布をサンプリング
  wsum0 <- matrix(0, nrow=hh, ncol=k)
  for(i in 1:hh){
    wsum0[i, ] <- colSums(Z[index_id[[i]], ])
  }
  wsum <- wsum0 + alpha01   #ディクレリ分布のパラメータ
  theta <- extraDistr::rdirichlet(hh, wsum)
  
  ##評価対象単語の出現率をサンプリング
  vf <- matrix(0, nrow=k, ncol=v)

  for(j in 1:k){
    index_z <- which(zw==j)
    x0 <- wd[index_z]
    vf[j, ] <- count(c(x0[y[index_z]==1], 1:v))[, 2] - 1 + beta01
  }
  phi <- extraDistr::rdirichlet(k, vf)   #ディクレリ分布からphiをサンプリング
  
  ##一般語の出現率をサンプリング
  lf <- count(c(wd[y==0], 1:v))[, 2] - 1 + gamma01
  lambda <- extraDistr::rdirichlet(1, lf)   #ディクレリ分布からphiをサンプリング

  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    
    #サンプリング結果を確認
    print(rp)
    print(round(cbind(theta[1:10, ], thetat[1:10, ]), 2)) 
    print(round(rbind(lambda[141:160], lambdat[141:160]), 3))
    print(round(cbind(phi[, 146:154], phit[, 146:154]), 3))
  }
}

cbind(wd, y_rate, y)

round(cbind(t(phi), t(phit)), 3)

round(phi[, 140:160], 3)
colMeans(phi)

lambda
lambdat

