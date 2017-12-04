#####トピック追跡モデル####
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
#set.seed(423943)
##データの設定
hh <- 250   #観測人数
pt <- 12   #観測期間
d <- hh*pt   #総文書数
k <- 10   #トピック数
v <- 250   #語彙数
w <- rpois(d, rgamma(d, 100, 0.8))   #1文書あたりの単語数

##IDの設定
id <- rep(1:hh, rep(pt, hh))
time <- rep(1:pt, hh)
ID <- data.frame(no=1:d, id, time)

##パラメータの設定
#ディクレリ事前分布を設定
alpha01 <- rep(0.4, k)   #1期目の文書のディクレリ事前分布のパラメータ
alpha02 <- rep(0.75, v)   #1期目の単語のディクレリ事前分布のパラメータ
alpha11 <- rep(10, k)   #2期目以降の文書のディクレリ事前分布のパラメータ
alpha12 <- rep(1000, v)   #2期目以降の単語のディクレリ事前分布のパラメータ


##時間ごとにディクレリ分布から文書と単語のトピック分布を生成
#1期目のトピック分布を生成
thetat <- theta <- matrix(0, nrow=d, ncol=k)
phit <- phi <- array(0, dim=c(k, v, pt))

index_time <- which(ID$time==1)
thetat[index_time, ] <- theta[index_time, ] <- extraDistr::rdirichlet(hh, alpha01)
phit[, , 1] <- phi[, , 1] <- extraDistr::rdirichlet(k, alpha02)

#2期目以降のトピック分布を逐次生成
for(j in 2:pt){
  index_time <- which(ID$time==j)
  
  #文書トピックを生成
  thetat[index_time, ] <- theta[index_time, ] <- 
    extraDistr::rdirichlet(hh, matrix(alpha11, nrow=hh, ncol=k, byrow=T) * thetat[index_time-1, ])
  thetat[index_time, ] <- theta[index_time, ]<- 
    (thetat[index_time, ] < 0.0001)*0.0001 + (thetat[index_time, ] >= 0.0001)*thetat[index_time, ]
  
  #単語トピックを生成
  phit[, , j] <- phi[, , j] <- 
    extraDistr::rdirichlet(k, matrix(alpha12, nrow=k, ncol=v, byrow=T) * phit[, , j-1])
  phit[, , j] <- phit[, , j]<- 
    (phit[, , j] < 0.00001)*0.00001 + (phit[, , j] >= 0.00001)*phit[, , j]
}

##多項分布から文書データを生成
WX <- matrix(0, nrow=d, ncol=v)
Z1 <- list()

for(i in 1:d){
  print(i)
  
  #文書のトピック分布を発生
  z1 <- extraDistr::rmnom(w[i], 1, theta[i, ])   #文書のトピック分布を発生
  
  #文書のトピック分布から単語を発生
  zn <- z1 %*% c(1:k)   #0,1を数値に置き換える
  zdn <- cbind(zn, z1)   #apply関数で使えるように行列にしておく
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi[x[1], , ID$time[i]])))   #文書のトピックから単語を生成
  wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
  WX[i, ] <- wdn  
  
  #発生させたトピックを格納
  Z1[[i]] <- zdn[, 1]
}
colSums(WX)

#データ行列を整数型行列に変更
storage.mode(WX) <- "integer"


####トピックモデルのためのデータと関数の準備####
##それぞれの文書中の単語の出現および補助情報の出現をベクトルに並べる
##データ推定用IDを作成
ID_list <- list()
hh_list <- list()
pt_list <- list()
wd_list <- list()

#求人ごとに求人IDおよび単語IDを作成
for(i in 1:nrow(WX)){
  print(i)
  
  #単語のIDベクトルを作成
  ID_list[[i]] <- rep(i, w[i])
  hh_list[[i]] <- rep(ID$id[i], w[i])
  pt_list[[i]] <- rep(ID$time[i], w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
  num2 <- subset(num1, num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
}


#リストをベクトルに変換
ID_d <- unlist(ID_list)
hh_d <- unlist(hh_list)
pt_d <- unlist(pt_list)
wd <- unlist(wd_list)

##インデックスを作成
doc1_list <- list()
doc2_list <- list()
doc3_list <- list()
word_list <- list()

for(i in 1:length(unique(ID_d))) {doc1_list[[i]] <- which(ID_d==i)}
for(i in 1:length(unique(hh_d))) {doc2_list[[i]] <- which(hh_d==i)}
for(i in 1:length(unique(pt_d))) {doc3_list[[i]] <- which(pt_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- which(wd==i)}

#時間ごとにインデックスを作成
ptdoc_list <- list()
for(j in 1:pt){
  pd <- list()
  for(i in 1:hh){
    pd[[i]] <- subset(doc3_list[[j]], doc3_list[[j]] %in% doc2_list[[i]])
  }
  ptdoc_list[[j]] <- pd
}

ptwd_list <- list()
for(j in 1:pt){
  pd <- list()
  for(i in 1:length(unique(wd))){
    pd[[i]] <- subset(doc3_list[[j]], doc3_list[[j]] %in% word_list[[i]])
  }
  ptwd_list[[j]] <- pd
}

gc(); gc()


####マルコフ連鎖モンテカルロ法で対応トピックモデルを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #負担係数の格納用
  for(kk in 1:k){
    #負担係数を計算
    Bi <- rep(theta[, kk], w) * phi[kk, c(wd)]   #尤度
    Bur[, kk] <- Bi   
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

##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- matrix(rep(1, k), nrow=hh, ncol=k, byrow=T)   
alpha02 <- rep(0.5, v)  
alpha11 <- matrix(rep(10, k), nrow=hh, ncol=k, byrow=T)   
alpha12 <- matrix(rep(10, v), nrow=k, ncol=v, byrow=T)

##パラメータの初期値
theta.ini <- runif(k, 0.5, 2)
phi.ini <- runif(v, 0.5, 1)
theta <- matrix(0, nrow=d, ncol=k)
phi <- array(0, dim=c(k, v, pt))
theta[ID$time==1, ] <- extraDistr::rdirichlet(hh, theta.ini)   #文書トピックのパラメータの初期値
phi[, , 1] <- extraDistr::rdirichlet(k, phi.ini)   #単語トピックのパラメータの初期値


##パラメータの格納用配列
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, pt, R/keep))
W_SEG <- matrix(0, nrow=R/(keep*10), ncol=sum(w))
storage.mode(W_SEG) <- "integer"
gc(); gc()

#インデックスを作成
index_id <- list()
index_pt <- list()
for(i in 1:hh){index_id[[i]] <- which(ID$id==i)}
for(i in 1:pt){index_pt[[i]] <- which(ID$time==i)}
vec <- 1/1:k


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##1期目の単語トピックをサンプリング
  #単語ごとにトピックの出現率を計算
  word_rate <- burden_fr(theta[index_pt[[1]], ], phi[, , 1], wd[doc3_list[[1]]], w[index_pt[[1]]], k)$Br
  
  #多項分布から単語トピックをサンプリング
  Zi <- matrix(0, nrow=length(wd), ncol=k)
  word_cumsums <- rowCumsums(word_rate)
  rand <- matrix(runif(nrow(word_rate)), nrow=nrow(word_rate), ncol=k)   #一様乱数
  Zi0 <- ((k+1) - (word_cumsums > rand) %*% rep(1, k)) %*% vec   #トピックをサンプリング
  Zi0[Zi0!=1] <- 0
  Zi[doc3_list[[1]], ] <- Zi0
  
  ##単語トピックのパラメータを更新
  #ディクレリ分布からthetaをサンプリング
  wsum0 <- matrix(0, nrow=hh, ncol=k)
  for(i in 1:hh){
    wsum0[i, ] <- colSums(Zi[doc2_list[[i]], ])
  }
  wsum <- wsum0 + alpha01   #ディクレリ分布のパラメータ
  theta[index_pt[[1]], ] <- extraDistr::rdirichlet(hh, wsum)   #ディクレリ分布からトピック割当をサンプリング
  
  system.time(
  wsum <- data.frame(id=ID_d, pt=pt_d, Zi=Zi) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(pt==1) %>%
    dplyr::summarize_all(funs(sum))
  )
  
  #ディクレリ分布からphiをサンプリング
  for(i in 1:v1){
    vf01[i, ] <- colSums(Zi1[word1_list[[i]], ])
  }
  vf <- t(vf01 + beta01m)   #ディクレリ分布のパラメータ
  phi1 <- extraDistr::rdirichlet(k, vf)   #ディクレリ分布からphiをサンプリング
  
  
  ##補助情報トピックをサンプリング
  #発生させた単語トピックからトピック抽出確率を計算
  aux_sums <- wsum - alpha01m
  theta_aux <- aux_sums/rowSums(aux_sums)
  
  #補助情報ごとにトピックの出現率を計算
  aux_rate <- burden_fr(theta_aux, omega, ad, x, k)$Br
  
  #多項分布から補助情報トピックをサンプリング
  aux_cumsums <- rowCumsums(aux_rate)
  rand <- matrix(runif(nrow(aux_rate)), nrow=nrow(aux_rate), ncol=k)   #一様乱数
  aux_z <- (k+1) - (aux_cumsums > rand) %*% rep(1, k)   #トピックをサンプリング
  
  #発生させたトピックをダミー変数化
  Zi2 <- matrix(0, nrow(aux_z), k)
  for(i in 1:k){
    index <- subset(1:nrow(aux_z), aux_z==i)
    Zi2[index, i] <- 1
  }
  
  ##補助情報トピックのパラメータを更新
  af <- as.matrix(data.frame(id=ad, Br=Zi2) %>%
                    dplyr::group_by(id) %>%
                    dplyr::summarize_all(funs(sum)))[, 2:(k+1)] + gamma0m
  omega <- t(apply(t(af), 1, function(x) rdirichlet(1, x)))
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    mkeep1 <- rp/keep
    THETA[, , mkeep1] <- theta
    PHI[, , mkeep1] <- phi
    OMEGA[, , mkeep1] <- omega
    if(rp%%(keep*10)==0){
      mkeep2 <- rp/(keep*5)
      W_SEG[mkeep2, ] <- word_z
      A_SEG[mkeep2, ] <- aux_z
    }
    
    #サンプリング結果を確認
    print(rp)
    print(round(cbind(theta[1:5, ], thetat[1:5, ]), 3))
    #print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
    #print(round(cbind(omega[, 1:10], omegat[, 1:10]), 3))
  }
}





