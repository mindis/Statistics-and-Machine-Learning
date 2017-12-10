#####Twitter LDA#####
library(MASS)
library(lda)
library(Matrix)
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
w <- rpois(d, 10.0)
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
alpha0 <- rep(0.25, k)   #ユーザー固有のディクレリ分布のパラメータ
alpha1 <- c(rep(0.4, v1), rep(0.005, v2))   #評価対象語のディクレリ分布のパラメータ
alpha2 <- c(rep(0.1, v1), rep(15, v2))   #一般語のディクレリ分布の事前分布のパラメータ
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
  
  for(j in 1:length(index_hh)){
    num1 <- WX[index_hh[j], ] * 1:v
    num2 <- which(num1 > 0)
    W1 <- WX[index_hh[j], (WX[index_hh[j], ] > 0)]
    wd_list[[index_hh[j]]] <- rep(num2, W1)
  }
}

#リストをベクトルに変換
ID_d <- unlist(ID_list)
td_d <- rep(1:length(w), w)
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

##スパース行列用のインデックスを作成
#複数観測された単語を記録
word_freq <- as.matrix(data.frame(td_d, wd) %>%  
                         dplyr::group_by(td_d) %>%
                         dplyr::count(wd))
dup <- word_freq[which(word_freq[, 3] > 1), ]

#複数観測された単語のインデックスを作成
multi_list <- list()
for(i in 1:nrow(dup)){
  print(i)
  multi_list[[i]] <- which(td_d==dup[i, 1] & wd==dup[i, 2])
}
index_multi <- unlist(multi_list)


##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 2000/keep

##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- 1  
beta01 <- 0.1  
gamma01 <- 0.1 

##パラメータの初期値
#tfidfで初期値を設定
#カスタムアイテム
tf <- WX/rowSums(WX)
idf1 <- log(nrow(WX)/colSums(WX > 0))
idf2 <- log(nrow(WX)/colSums(WX == 0))
theta <- extraDistr::rdirichlet(hh, rep(1, k))   #ユーザートピックの初期値
phi <- extraDistr::rdirichlet(k, idf1)   #評価対象語の出現確率の初期値
lambda <- extraDistr::rdirichlet(1, idf2*50)   #一般語の出現確率の初期値
y <- rbinom(f, 1, 0.5)
r <- c(0.5, 0.5)   #評価対象語と一般語の割当率の混合率

##パラメータの格納用配列
THETA <- array(0, dim=c(hh, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
LAMBDA <- matrix(0, nrow=R/keep, ncol=v)
SEG <- matrix(0, nrow=d, ncol=k)
Y <- matrix(0, nrow=f, ncol=2)
storage.mode(SEG) <- "integer"

##MCMC推定用配列
Zi <- matrix(0, nrow=f, ncol=k)
LH1 <- rep(0, f)
wsum0 <- matrix(0, nrow=hh, ncol=k)
vf <- matrix(0, nrow=k, ncol=v)


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##一般語のスパース行列を定義
  #複数データのスパース行列インデックス
  index_y0 <- index_multi[y[index_multi]==0]
  multi_word <- data.frame(td_d=td_d[index_y0], wd=wd[index_y0]) %>%  
    dplyr::group_by(td_d) %>%
    dplyr::count(wd)
  
  #単一データのスパース行列インデックス
  index_y0 <- which(y[-index_multi]==0)
  td_ones <- td_d[-index_multi][index_y0]
  wd_ones <- wd[-index_multi][index_y0]
  
  #ベクトルを結合
  td_vec <- c(td_ones, multi_word$td_d)
  wd_vec <- c(wd_ones, multi_word$wd)
  freq_vec <- c(rep(1, length(td_ones)), multi_word$n)
  
  #スパース行列を作成
  sparse_gword <- sparseMatrix(td_vec, wd_vec, x=freq_vec, dims=c(d, v))
  
  ##tweetトピックをサンプリング
  #tweetごとにトピックの出現確率を計算
  Data_sparse <- WX_sparse - sparse_gword   #一般語を除く
  Br <- burden_fr(theta, phi, Data_sparse, u_id)   #尤度を計算
  
  #多項分布から単語トピックをサンプリング
  tweet_rate <- Br$Br
  Z <- rmnom(d, 1, tweet_rate)
  z <- Z %*% 1:k
  zw <- rep(z, w)
  
  ##単語がトピックと関係あるかどうかをサンプリング
  #潜在変数のパラメータ
  LH0 <- lambda[wd]
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
  
  #混合率を更新
  r0 <- mean(y)
  r <- c(r0, 1-r0)
  
  ##ディクレリ分布からトピック分布をサンプリング
  for(i in 1:hh){
    wsum0[i, ] <- colSums(Z[index_id[[i]], ])
  }
  wsum <- wsum0 + alpha01   #ディクレリ分布のパラメータ
  theta <- extraDistr::rdirichlet(hh, wsum)
  
  ##評価対象単語の出現率をサンプリング
  for(j in 1:k){
    index_z <- which(zw==j)
    x0 <- wd[index_z]
    vf[j, ] <- count(c(x0[y[index_z]==1], 1:v))[, 2] - 1 + beta01
  }
  phi <- extraDistr::rdirichlet(k, vf)   #ディクレリ分布からphiをサンプリング
  
  ##一般語の出現率をサンプリング
  lf <- count(c(wd[y==0], 1:v))[, 2] - 1 + gamma01
  lambda <- extraDistr::rdirichlet(1, lf)   #ディクレリ分布からlambdaをサンプリング
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    LAMBDA[mkeep, ] <- lambda
    if(burnin > rp){
      SEG <- SEG + Z
      Y <- Y + cbind(y, 1-y)
    }
    
    #サンプリング結果を確認
    print(rp)
    print(round(cbind(theta[1:10, ], thetat[1:10, ]), 2)) 
    print(round(rbind(lambda[141:160], lambdat[141:160]), 3))
    print(round(cbind(phi[, 146:154], phit[, 146:154]), 3))
  }
}

####サンプリング結果の可視化と要約####
burnin <- 2000/keep
RS <- R/keep

##サンプリング結果をプロット
matplot(t(THETA[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA[100, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA[1000, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI[, 1, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI[, 150, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI[, 151, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(LAMBDA[, 146:150], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(LAMBDA[, 151:155], type="l", xlab="サンプリング回数", ylab="パラメータ")

##パラメータのサンプリング結果の要約
round(cbind(apply(THETA[, , burnin:RS], c(1, 2), mean), thetat), 3)   #ユーザーのトピック割当確率
round(cbind(t(apply(PHI[, , burnin:RS], c(1, 2), mean)), t(phit)), 3)   #評価対象語のトピック別の出現確率
round(cbind(colMeans(LAMBDA[burnin:RS, ]), as.numeric(lambdat)), 3)   #一般語の出現確率

##トピックのサンプリング結果の要約
round(cbind(SEG/rowSums(SEG), z_vec), 3)   #tweetごとのトピック割当の要約
round(cbind(Y/rowSums(Y), wd), 3)   #一般語かどうかの潜在変数の要約


