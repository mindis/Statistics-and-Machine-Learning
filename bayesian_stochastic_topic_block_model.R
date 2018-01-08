#####トピックモデルによる確率的ブロックモデル#####
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

#set.seed(2578)

####データの発生####
##データの設定
d <- 10000   #ユーザー数
k1 <- 10   #ユーザーのセグメント数
k2 <- 12   #トピック数
v <- 1000   #アイテム数
w <- rpois(d, rgamma(d, 15, 0.5))   #購買アイテム数
w[w < 5] <- ceiling(runif(sum(w < 5), 5, 15))
f <- sum(w)

##パラメータの設定
#ディクレリ分布のパラメータを設定
alpha01 <- rep(0.25, k2)   #トピック分布のパラメータ
alpha11 <- rep(0.1, v)   #アイテム購買確率のパラメータ

#ディクレリ分布のパラメータを生成
omega <- omegat <- rep(1/k1, k1)   #混合率
theta <- thetat <- extraDistr::rdirichlet(k1, alpha01)   #トピック分布のパラメータ
phi <- phit <- extraDistr::rdirichlet(k2, alpha11)   #単語分布のパラメータ


##多項分布よりアイテム購買行列を生成
WX <- matrix(0, nrow=d, ncol=v)
Z1_list <- list()
Z2_list <- list()

for(i in 1:d){
  #ユーザーのセグメントを生成
  z1 <- rmnom(1, 1, omega)
  z1_vec <- as.numeric(z1 %*% 1:k1)
  
  #生成したユーザーセグメントからトピックを生成
  z2 <- rmnom(w[i], 1, theta[z1_vec, ])
  z2_vec <- as.numeric(z2 %*% 1:k2)
  
  #トピックからアイテム購買を生成
  wd <- rmnom(w[i], 1, phi[z2_vec, ])
  WX[i, ] <- colSums(wd)
  
  #パラメータを格納
  Z1_list[[i]] <- z1
  Z2_list[[i]] <- z2
}


#リスト形式を変換
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)
v <- length(which(colSums(WX) > 0))
index <- which(colSums(WX) > 0)
WX <- WX[, index]
phit <- phit[, index]
const <- lfactorial(w) - rowSums(lfactorial(WX))   #多項分布の密度関数の対数尤度の定数

##データ推定用IDを作成
ID_list <- list()
wd_list <- list()

#求人ごとに求人IDおよび単語IDを作成
for(i in 1:nrow(WX)){
  
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
for(i in 1:d) {doc_list[[i]] <- which(ID_d==i)}
for(i in 1:v) {word_list[[i]] <- which(wd==i)}
gc(); gc()


####マルコフ連鎖モンテカルロ法で確率的トピックブロックモデルを推定####
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

##潜在変数zを計算する関数
LLobz <- function(WX, theta, r, const, d, k){
  
  #多項分布の対数尤度
  log_theta <- log(t(theta))
  LLi <- const + WX %*% log_theta
  
  #logsumexpの尤度
  LLi_max <- matrix(apply(LLi, 1, max), nrow=d, ncol=k)
  r_matrix <- matrix(r, nrow=d, ncol=k, byrow=T)
  
  #割当確率のパラメータを設定
  expl <- r_matrix * exp(LLi - LLi_max)
  expl_log <- log(expl)
  expl_max <- matrix(log(max(expl[1, ])), nrow=d, ncol=k)
  z <- exp(expl_log - (log(rowSums(exp(expl_log - expl_max))) + expl_max))   #セグメント割当確率
}

##アルゴリズムの設定
R <- 10000
keep <- 2
burnin <- 1000/keep
iter <- 0
disp <- 10
LLt <- sum(log(rowSums(burden_fr(thetat[as.numeric(Z1 %*% 1:k1), ], phit, wd, w, k2)$Bur)))

##事前分布の設定
alpha1 <- 1
beta1 <- 1

##初期値の設定
r <- rep(1/k1, k1)   #混合率の初期値
Zi1 <- rmnom(d, 1, rep(1/k1, k1))   #セグメント割当の初期値
theta <- extraDistr::rdirichlet(k1, rep(100, k2))   #トピック分布の初期値
phi <- extraDistr::rdirichlet(k2, colSums(WX)/sum(WX)*100)   #単語分布の初期値

##パラメータの格納用配列
THETA <- array(0, dim=c(k1, k2, R/keep))
PHI <- array(0, dim=c(k2, v, R/keep))
SEG1 <- matrix(0, nrow=d, ncol=k1)
SEG2 <- matrix(0, nrow=f, ncol=k2)
Zi1 <- Z1


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##アイテムごとにトピックをサンプリング
  #トピック分布の割当確率を計算
  zi1 <- as.numeric(Zi1 %*% 1:k1)
  thetan <- theta[zi1, ]
  out <- burden_fr(thetan, phi, wd, w, k2)
  word_rate <- out$Br
  
  #多項分布からトピックをサンプリング
  Zi2 <- rmnom(f, 1, word_rate)
  z2_vec <- as.numeric(Zi2 %*% 1:k2)
  
  
  ##ユーザーごとにセグメントをサンプリング
  #セグメント割当確率を計算
  ZX <- matrix(0, nrow=d, ncol=k2)
  for(i in 1:d){
   ZX[i, ] <- colSums(Zi2[doc_list[[i]], , drop=FALSE])
  }
  LLi0 <- ZX %*% t(log(theta))
  LLho <- exp(LLi0 - apply(LLi0, 1, max))
  z1_rate <- r*LLho / rowSums(r*LLho)
  
  #多項分布からセグメントをサンプリング
  Zi1 <- rmnom(d, 1, z1_rate)
  z1_vec <- as.numeric(Zi1 %*% 1:k1)
  
  #混合率を更新
  #r <- matrix(colMeans(Zi1), nrow=d, ncol=k1, byrow=T)
  r <- matrix(1/k1, nrow=d, ncol=k1)
  
  ##パラメータを更新
  #トピック分布thetaをサンプリング
  wsum0 <- matrix(0, nrow=k1, ncol=k2)
  for(j in 1:k1){
    wsum0[j, ] <- colSums(ZX[Zi1[, j]==1, , drop=FALSE])
  }
  wsum <- wsum0 + alpha1
  theta <- extraDistr::rdirichlet(k1, wsum)
  
  #単語分布phiをサンプリング
  vf0 <- matrix(0, nrow=k2, ncol=v)
  for(j in 1:v){
    vf0[, j] <- colSums(Zi2[word_list[[j]], , drop=FALSE])
  }
  vf <- vf0 + beta1
  phi <- extraDistr::rdirichlet(k2, vf)
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(rp%%keep==0 & rp >= burnin){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
    }
    
    #サンプリング結果を確認
    if(rp%%disp==0){
      print(rp)
      print(c(sum(log(rowSums(out$Bur))), LLt))
      print(round(rbind(r=colMeans(Zi1), omega), 3))
      #print(round(cbind(theta, thetat), 3))
      print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
    }
  }
}

matplot(t(THETA[1, , ]), type="l")
matplot(t(THETA[2, , ]), type="l")
matplot(t(THETA[3, , ]), type="l")
matplot(t(PHI[, 1, ]), type="l")
matplot(t(PHI[, 2, ]), type="l")
matplot(t(PHI[, 3, ]), type="l")
