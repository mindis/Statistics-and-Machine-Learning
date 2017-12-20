#####教師ありトピックモデル####
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


####データの発生####
#set.seed(423943)
#データの設定
k <- 6   #教師数
d <- 3000   #文書数
v1 <- 300   #教師に関係のある語彙数
v2 <- 100   #教師に関係のない語彙数
v <- v1 + v2   #語彙数
w <- rpois(d, rgamma(d, 55, 0.50))   #1文書あたりの単語数
f <- sum(w)

#インデックスを作成
id_vec <- rep(1:d, w)
index_hh <- list()
for(i in 1:d){
  index_hh[[i]] <- which(id_vec==i)
}

#教師データを生成
pr <- runif(k, 1, 5)
Y <- rmnom(d, 1, pr)
y <- as.numeric(Y %*% 1:k)
y_vec <- as.numeric(y)[id_vec] 


#パラメータの設定
alpha0 <- alpha0 <- diag(1.5, k) + 0.15   #文書のディレクリ事前分布のパラメータ
alpha1 <- c(rep(0.4, v1), rep(0.005, v2))   #評点に関係のある単語のディレクリ事前分布のパラメータ
alpha2 <- c(rep(0.05, v1), rep(5, v2))   #教師に関係のない単語のディレクリ事前分布のパラメータ

#ディレクリ乱数の発生
thetat <- theta <- rdirichlet(d, alpha0[y, ])   #文書のトピック分布をディレクリ乱数から発生
phit <- phi <- rdirichlet(k, alpha1)   #評点に関係のある単語分布をディレクリ乱数から発生
gammat <- gamma <- rdirichlet(1, alpha2)   #評点に関係のない単語分布をディレクリ乱数から発生
betat <- beta <- rbeta(sum(f), 15, 15)   #単語が評点と関連するかどうかのパラメータ


##多項分布の乱数からデータを発生
WX <- matrix(0, nrow=d, ncol=v)
x_list <- list()
x <- rep(0, f)
Z <- list()
index_v1 <- 1:v1
index_v2 <- (v1+1):v

for(i in 1:d){
  #文書のトピックを生成
  z <- rmnom(w[i], 1, theta[i, ])   #文書のトピック分布を発生
  z_vec <- z %*% c(1:k)   #トピック割当をベクトル化
  
  #一般語かどうかを生成
  x_list[[i]] <- rbinom(w[i], 1, beta[index_hh[[i]]])
  
  #生成したトピックから単語を生成
  wn <- rmnom(sum(x_list[[i]]), 1, phi[z_vec[x_list[[i]]==1], ])   #文書のトピックから単語を生成
  an <- rmnom(sum(1-x_list[[i]]), 1, gammat)
  wdn <- colSums(wn) + colSums(an)   #単語ごとに合計して1行にまとめる
  WX[i, ] <- wdn
  Z[[i]] <- z
  x[index_hh[[i]]] <- x_list[[i]]
  print(i)
}

####トピックモデル推定のためのデータと関数の準備####
##それぞれの文書中の単語の出現および補助情報の出現をベクトルに並べる
##データ推定用IDを作成
ID_list <- list()
wd_list <- list()

#求人ごとに求人IDおよび単語IDを作成
for(i in 1:nrow(WX)){
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
R <- 10000   #サンプリング回数
keep <- 2   #2回に1回の割合でサンプリング結果を格納
iter <- 0
burnin <- 1000/keep

##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- 1.0
alpha02 <- 0.5


##パラメータの初期値
#トピック分布の初期値
theta <- extraDistr::rdirichlet(d, (diag(1.0, k) + 0.25)[y, ])   #文書トピックのパラメータの初期値

#一般語の単語分布の初期値
inv_idf <- colSums(WX > 0)/d
gamma <- inv_idf / sum(inv_idf)

#教師関連単語分布の初期値
phi <- matrix(0, nrow=k, ncol=v)
for(j in 1:k) {
  M <- colSums(WX[y==j, ])
  phi[j, ] <- M*log(1/inv_idf) / sum(M*log(1/inv_idf))
}


##パラメータの格納用配列
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
GAMMA <- matrix(0, nrow=R/keep, ncol=v)
SEG1 <- matrix(0, nrow=f, ncol=k)
SEG2 <- rep(0, nrow=f)
storage.mode(SEG1) <- "integer"
gc(); gc()

##MCMC推定用配列
wsum0 <- matrix(0, nrow=d, ncol=k)
vf0 <- matrix(0, nrow=k, ncol=v)


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

rep(theta[, j], w) 
phi[y_vec, ]


y_vec
M <- matrix(0, nrow=f, ncol=k+1)
for(j in 1:k){
  M[, j] <- phit[j, wd]
}
M[, k+1] <- gammat[wd] 

round(cbind(y_vec, wd, M/rowSums(M)), 3)

####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##単語トピックをサンプリング
  #単語ごとにトピックの出現率を計算
  word_rate <- burden_fr(theta, phi, wd, w, k)$Br
  
  #多項分布から単語トピックをサンプリング
  Zi <- rmnom(f, 1, word_rate)   
  z_vec <- Zi %*% 1:k
  
  ##単語トピックのパラメータを更新
  #ディクレリ分布からthetaをサンプリング
  for(i in 1:d){
    wsum0[i, ] <- colSums(Zi[doc_list[[i]], ])
  }
  wsum <- wsum0 + alpha01m 
  theta <- extraDistr::rdirichlet(d, wsum)
  
  #ディクレリ分布からphiをサンプリング
  for(j in 1:v){
    vf0[, j] <- colSums(Zi[word_list[[j]], ])
  }
  vf <- vf0 + beta0m
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
    print(rp)
    print(round(cbind(theta[1:10, ], thetat[1:10, ]), 3))
    print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
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
round(cbind(topic_mu, thetat), 3)
round(topic_sd <- apply(THETA[, , burnin:(R/keep)], c(1, 2), sd), 3)   #トピック分布の事後標準偏差

#単語出現確率の事後推定量
word_mu <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #単語の出現率の事後平均
round(rbind(word_mu, phit)[, 1:50], 3)
