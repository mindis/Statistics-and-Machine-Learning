#####対応トピックモデル#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
detach("package:gtools", unload=TRUE)
detach("package:bayesm", unload=TRUE)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

set.seed(8079)

####データの発生####
#set.seed(423943)
#文書データの設定
k <- 8   #トピック数
d <- 5000   #文書数
v <- 300   #語彙数
w <- rpois(d, 200)   #1文書あたりの単語数
a <- 25   #補助変数数
x0 <- rpois(d, 10)
x <- ifelse(x0 < 1, 1, x0)

#パラメータの設定
alpha0 <- round(runif(k, 0.1, 1.25), 3)   #文書のディレクリ事前分布のパラメータ
alpha1 <- rep(0.25, v)   #単語のディレクリ事前分布のパラメータ
alpha2 <- rep(0.3, a)   #補助情報のディクレリ事前分布のパラメータ

#ディレクリ乱数の発生
thetat <- theta <- rdirichlet(d, alpha0)   #文書のトピック分布をディレクリ乱数から発生
phit <- phi <- rdirichlet(k, alpha1)   #単語のトピック分布をディレクリ乱数から発生
lambda <- matrix(0, nrow=d, ncol=k)   #文書に含むトピックだけを補助情報のトピックにするための確率を格納する行列
omegat <- omega <- rdirichlet(k, alpha2)   #補助情報のトピック分布をディクレリ乱数から発生

#多項分布の乱数からデータを発生
WX <- matrix(0, nrow=d, ncol=v)
AX <- matrix(0, nrow=d, ncol=a)
Z1 <- list()
Z2 <- list()

for(i in 1:d){
  print(i)
  
  #文書のトピック分布を発生
  z1 <- t(rmultinom(w[i], 1, theta[i, ]))   #文書のトピック分布を発生
  
  #文書のトピック分布から単語を発生
  zn <- z1 %*% c(1:k)   #0,1を数値に置き換える
  zdn <- cbind(zn, z1)   #apply関数で使えるように行列にしておく
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi[x[1], ])))   #文書のトピックから単語を生成
  wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
  WX[i, ] <- wdn  
  
  #文書のトピック分布から補助変数を発生
  #文書で発生させたトピックのみを補助情報のトピック分布とする
  rate <- rep(0, k)
  z_table <- table(zn)
  lambda[i, as.numeric(names(z_table))] <- z_table/sum(z_table)
  
  z2 <- t(rmultinom(x[i], 1, lambda[i, ]))  
  zx <- z2 %*% 1:k
  zax <- cbind(zx, z2)
  an <- t(apply(zax, 1, function(x) rmultinom(1, 1, omega[x[1], ])))
  adn <- colSums(an)
  AX[i, ] <- adn
  
  #文書トピックおよび補助情報トピックを格納
  Z1[[i]] <- zdn[, 1]
  Z2[[i]] <- zax[, 1]
}

#データ行列を整数型行列に変更
storage.mode(WX) <- "integer"
storage.mode(AX) <- "integer"


####EMアルゴリズムでトピックモデルを推定####
####トピックモデルのためのデータと関数の準備####
##それぞれの文書中の単語の出現および補助情報の出現をベクトルに並べる
##データ推定用IDを作成
ID1_list <- list()
wd_list <- list()
ID2_list <- list()
ad_list <- list()

#求人ごとに求人IDおよび単語IDを作成
for(i in 1:nrow(WX)){
  print(i)
  
  #単語のIDベクトルを作成
  ID1_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
  num2 <- subset(num1, num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
  
  #補助情報のIDベクトルを作成
  ID2_list[[i]] <- rep(i, x[i])
  num1 <- (AX[i, ] > 0) * (1:a)
  num2 <- subset(num1, num1 > 0)
  A1 <- AX[i, (AX[i, ] > 0)]
  number <- rep(num2, A1)
  ad_list[[i]] <- number
}

#リストをベクトルに変換
ID1_d <- unlist(ID1_list)
ID2_d <- unlist(ID2_list)
wd <- unlist(wd_list)
ad <- unlist(ad_list)

##インデックスを作成
doc1_list <- list()
word_list <- list()
doc2_list <- list()
aux_list <- list()
for(i in 1:length(unique(ID1_d))) {doc1_list[[i]] <- subset(1:length(ID1_d), ID1_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- subset(1:length(wd), wd==i)}
for(i in 1:length(unique(ID2_d))) {doc2_list[[i]] <- subset(1:length(ID2_d), ID2_d==i)}
for(i in 1:length(unique(ad))) {aux_list[[i]] <- subset(1:length(ad), ad==i)}
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
alpha01 <- rep(1.0, k)
beta0 <- rep(0.5, v)
gamma0 <- rep(0.5, a)
alpha01m <- matrix(alpha01, nrow=d, ncol=k, byrow=T)
beta0m <- matrix(beta0, nrow=v, ncol=k)
gamma0m <- matrix(gamma0, nrow=a, ncol=k)

##パラメータの初期値
theta.ini <- runif(k, 0.5, 2)
phi.ini <- runif(v, 0.5, 1)
omega.ini <- runif(a, 0.5, 1)
theta <- rdirichlet(d, theta.ini)   #文書トピックのパラメータの初期値
phi <- rdirichlet(k, phi.ini)   #単語トピックのパラメータの初期値
omega <- rdirichlet(k, omega.ini)   #補助情報トピックのパラメータの初期値

##パラメータの格納用配列
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
OMEGA <- array(0, dim=c(k, a, R/keep))
W_SEG <- matrix(0, nrow=R/(keep*10), ncol=sum(w))
A_SEG <- matrix(0, nrow=R/(keep*10), ncol=sum(x))
storage.mode(W_SEG) <- "integer"
storage.mode(A_SEG) <- "integer"
gc(); gc()


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##単語トピックをサンプリング
  #単語ごとにトピックの出現率を計算
  word_rate <- burden_fr(theta, phi, wd, w, k)$Br
  
  #多項分布から単語トピックをサンプリング
  word_cumsums <- rowCumsums(word_rate)
  rand <- matrix(runif(nrow(word_rate)), nrow=nrow(word_rate), ncol=k)   #一様乱数
  word_z <- (k+1) - (word_cumsums > rand) %*% rep(1, k)   #トピックをサンプリング
  cbind(word_cumsums > rand, (k+1) - (word_cumsums > rand) %*% rep(1, k))
  
  
  #発生させたトピックをダミー変数化
  Zi1 <- matrix(0, nrow(word_z), k)
  for(i in 1:k){
    index <- subset(1:nrow(word_z), word_z==i)
    Zi1[index, i] <- 1
  }
  
  ##単語トピックのパラメータを更新
  #ディクレリ分布からthetaをサンプリング
  wsum <- as.matrix(data.frame(id=ID1_d, Br=Zi1) %>%
                      dplyr::group_by(id) %>%
                      dplyr::summarize_all(funs(sum)))[, 2:(k+1)] + alpha01m
  theta <- t(apply(wsum, 1, function(x) rdirichlet(1, x)))

  #ディクレリ分布からphiをサンプリング
  vf <- as.matrix(data.frame(id=wd, Br=Zi1) %>%
                    dplyr::group_by(id) %>%
                    dplyr::summarize_all(funs(sum)))[, 2:(k+1)] + beta0m
  phi <- t(apply(t(vf), 1, function(x) rdirichlet(1, x)))
  
  
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

