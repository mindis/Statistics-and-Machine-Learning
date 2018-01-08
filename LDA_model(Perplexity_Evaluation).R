#####Latent Dirichlet Allocationモデル(Perplexity)#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(21437)

####データの発生####
#set.seed(423943)
#データの設定
k <- 15   #トピック数
d <- 3000   #文書数
v <- 500   #語彙数
w <- rpois(d, rgamma(d, 70, 0.40))   #1文書あたりの単語数
f <- sum(w)


#パラメータの設定
alpha0 <- rep(0.25, k)   #文書のディレクリ事前分布のパラメータ
alpha1 <- rep(0.1, v)   #単語のディレクリ事前分布のパラメータ

#ディレクリ乱数の発生
thetat <- theta <- rdirichlet(d, alpha0)   #文書のトピック分布をディレクリ乱数から発生
phit <- phi <- rdirichlet(k, alpha1)   #単語のトピック分布をディレクリ乱数から発生

#多項分布の乱数からデータを発生
WX <- matrix(0, nrow=d, ncol=v)
Z <- list()

for(i in 1:d){
  z <- rmnom(w[i], 1, theta[i, ])   #文書のトピック分布を発生
  z_vec <- z %*% c(1:k)   #トピック割当をベクトル化
  
  wn <- rmnom(w[i], 1, phi[z_vec, ])   #文書のトピックカラ単語を生成
  wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
  WX[i, ] <- wdn
  Z[[i]] <- z
  print(i)
}

####トピックモデル推定のためのデータと関数の準備####
##それぞれの文書中の単語の出現および補助情報の出現をベクトルに並べる
##データ推定用IDを作成
n1 <- 2000   #学習用サンプル数
n2 <- 1000   #検証用サンプル数
WX1 <- WX[1:n1, ]   #学習用データ
w1 <- w[1:n1]
WX2 <- WX[(n1+1):d, ]   #検証用データ
w2 <- w[(n1+1):d]
ID_list <- list()
wd_list <- list()

#文書ごとに文書IDおよび単語IDを作成
for(i in 1:nrow(WX1)){
  print(i)
  
  #単語のIDベクトルを作成
  ID_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
  num2 <- which(num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
}

#リストをベクトルに変換
ID_d <- unlist(ID_list)
wd <- unlist(wd_list)

##インデックスを作成
doc_list <- list()
word_list <- list()
for(i in 1:length(unique(ID_d))) {doc_list[[i]] <- which(ID_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- which(wd==i)}
gc(); gc()


####マルコフ連鎖モンテカルロ法で対応トピックモデルを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #負担係数の格納用
  for(j in 1:k){
    #負担係数を計算
    Bi <- rep(theta[, j], w) * phi[j, wd]   #尤度
    Bur[, j] <- Bi   
  }
  
  Br <- Bur / rowSums(Bur)   #負担率の計算
  r <- colSums(Br) / sum(Br)   #混合率の計算
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##アルゴリズムの設定
R <- 2500   #サンプリング回数
keep <- 2   #2回に1回の割合でサンプリング結果を格納
iter <- 0
burnin <- 500/keep
disp <- 10

##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- 1.0
beta0 <- 0.5


##パラメータの初期値
theta <- rdirichlet(n1, rep(1, k))   #文書トピックのパラメータの初期値
phi <- rdirichlet(k, colSums(WX1)/sum(WX1)*v)   #単語トピックのパラメータの初期値

##パラメータの格納用配列
f1 <- length(ID_d)
THETA <- array(0, dim=c(n1, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
SEG <- matrix(0, nrow=f1, ncol=k)
storage.mode(SEG) <- "integer"
gc(); gc()

##MCMC推定用配列
wsum0 <- matrix(0, nrow=n1, ncol=k)
vf0 <- matrix(0, nrow=k, ncol=v)


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##単語トピックをサンプリング
  #単語ごとにトピックの出現率を計算
  word_rate <- burden_fr(theta, phi, wd, w1, k)$Br
  
  #多項分布から単語トピックをサンプリング
  Zi <- rmnom(f1, 1, word_rate)   
  z_vec <- Zi %*% 1:k
  
  ##単語トピックのパラメータを更新
  #ディクレリ分布からthetaをサンプリング
  for(i in 1:n1){
    wsum0[i, ] <- colSums(Zi[doc_list[[i]], ])
  }
  wsum <- wsum0 + alpha01 
  theta <- extraDistr::rdirichlet(n1, wsum)
  
  #ディクレリ分布からphiをサンプリング
  for(j in 1:v){
    vf0[, j] <- colSums(Zi[word_list[[j]], ])
  }
  vf <- vf0 + beta0
  phi <- extraDistr::rdirichlet(k, vf)
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(rp%%keep==0 & rp >= burnin){
      SEG <- SEG + Zi
    }
    
    #サンプリング結果を確認
    if(rp%%disp==0){
      print(rp)
      #print(round(cbind(theta[1:10, ], thetat[1:10, ]), 3))
      print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
    }
  }
}


####サンプリング結果の可視化と要約####
burnin <- 1000/keep   #バーンイン期間
RS <- R/keep

##サンプリング結果の可視化
#文書のトピック分布のサンプリング結果
matplot(t(THETA[1, , ]), type="l", ylab="パラメータ", main="文書1のトピック分布のサンプリング結果")
matplot(t(THETA[2, , ]), type="l", ylab="パラメータ", main="文書2のトピック分布のサンプリング結果")
matplot(t(THETA[3, , ]), type="l", ylab="パラメータ", main="文書3のトピック分布のサンプリング結果")
matplot(t(THETA[4, , ]), type="l", ylab="パラメータ", main="文書4のトピック分布のサンプリング結果")

#単語の出現確率のサンプリング結果
matplot(t(PHI[1, 1:10, ]), type="l", ylab="パラメータ", main="トピック1の単語の出現率のサンプリング結果")
matplot(t(PHI[2, 11:20, ]), type="l", ylab="パラメータ", main="トピック2の単語の出現率のサンプリング結果")
matplot(t(PHI[3, 21:30, ]), type="l", ylab="パラメータ", main="トピック3の単語の出現率のサンプリング結果")
matplot(t(PHI[4, 31:40, ]), type="l", ylab="パラメータ", main="トピック4の単語の出現率のサンプリング結果")
matplot(t(PHI[5, 41:50, ]), type="l", ylab="パラメータ", main="トピック5の単語の出現率のサンプリング結果")
matplot(t(PHI[6, 51:60, ]), type="l", ylab="パラメータ", main="トピック6の単語の出現率のサンプリング結果")
matplot(t(PHI[7, 61:70, ]), type="l", ylab="パラメータ", main="トピック7の単語の出現率のサンプリング結果")
matplot(t(PHI[8, 71:80, ]), type="l", ylab="パラメータ", main="トピック8の単語の出現率のサンプリング結果")

##サンプリング結果の要約推定量
#トピック分布の事後推定量
topic_mu <- apply(THETA[, , burnin:(R/keep)], c(1, 2), mean)   #トピック分布の事後平均
round(cbind(topic_mu, thetat[1:n1, ]), 3)
round(topic_sd <- apply(THETA[, , burnin:(R/keep)], c(1, 2), sd), 3)   #トピック分布の事後標準偏差

#単語出現確率の事後推定量
word_mu <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #単語の出現率の事後平均
round(rbind(word_mu, phit)[, 1:50], 3)


####テストデータを用いてPerplexityを評価####
##インデックスの作成
ID2_list <- list()
wd2_list <- list()

#文書ごとに文書IDおよび単語IDを作成
for(i in 1:nrow(WX2)){
  print(i)
  
  #単語のIDベクトルを作成
  ID2_list[[i]] <- rep(i, w2[i])
  num1 <- (WX2[i, ] > 0) * (1:v)
  num2 <- which(num1 > 0)
  W1 <- WX2[i, (WX2[i, ] > 0)]
  number <- rep(num2, W1)
  wd2_list[[i]] <- number
}

#リストをベクトルに変換
ID2_d <- unlist(ID2_list)
wd2 <- unlist(wd2_list)

##インデックスを作成
doc2_list <- list()
word2_list <- list()
for(i in 1:length(unique(ID2_d))) {doc2_list[[i]] <- which(ID2_d==i)}
for(i in 1:length(unique(wd2))) {word2_list[[i]] <- which(wd2==i)}
gc(); gc()

#アルゴリズムの設定
R <- 2500   #サンプリング回数
keep <- 2   #2回に1回の割合でサンプリング結果を格納
iter <- 0
burnin <- 500/keep
disp <- 10
f2 <- length(wd2)

#初期値を設定
theta <- extraDistr::rdirichlet(n2, rep(1, k))

#パラメータの格納用配列
THETA2 <- array(0, dim=c(n2, k, R/keep))
SEG2 <- matrix(0, nrow=f2, k)

for(rp in 1:R){
  
  ##単語トピックをサンプリング
  #単語ごとにトピックの出現率を計算
  word_rate <- burden_fr(theta, word_mu, wd2, w2, k)$Br
  
  #多項分布から単語トピックをサンプリング
  Zi <- rmnom(f2, 1, word_rate)   
  z_vec <- Zi %*% 1:k
  
  ##単語トピックのパラメータを更新
  #ディクレリ分布からthetaをサンプリング
  wsum0 <- matrix(0, nrow=n2, ncol=k)
  for(i in 1:n2){
    wsum0[i, ] <- colSums(Zi[doc2_list[[i]], ])
  }
  wsum <- wsum0 + alpha01 
  theta <- extraDistr::rdirichlet(n2, wsum)

  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA2[, , mkeep] <- theta
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(rp%%keep==0 & rp >= burnin){
      SEG2 <- SEG2 + Zi
    }
    
    #サンプリング結果を確認
    if(rp%%disp==0){
      print(rp)
      print(round(cbind(theta[1:10, ], thetat[(n1+1):(n1+10), ]), 3))
    }
  }
}

##Perplexityを計算
#ユニグラムモデルのPerplexityを計算
par <- colSums(WX1) / sum(WX1)
LL1 <- sum(WX2 %*% log(par))
Perx1 <- exp(-LL1 / f2)   #Perplexity

#トピックモデルのPerplexityを計算
theta_mu <- apply(THETA2[, , burnin:(R/keep)], c(1, 2), mean)   #thetaの事後平均
LL2 <- sum(log(rowSums(burden_fr(theta_mu, word_mu, wd2, w2, k)$Bur)))   #対数尤度
Perx2 <- exp(-LL2 / f2)   #Perplexity
cbind(Perx1, Perx2)

