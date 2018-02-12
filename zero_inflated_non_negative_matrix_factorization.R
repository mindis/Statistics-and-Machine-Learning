#####ゼロ過剰非負値行列因子分解#####
options(warn=2)
library(MASS)
library(NMF)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(data.table)
library(ggplot2)

#set.seed(319547)

####データの発生####
##データの設定
hh <- 3000   #ユーザー数
item <- 500   #アイテム数
k <- 10   #潜在変数数

##ゼロ過剰NMFの仮定に基づきデータを生成
#ガンマ分布よりパラメータを設定
alpha01 <- 0.25; beta01 <- 1.0
alpha02 <- 0.25; beta02 <- 0.8
W <- WT <- matrix(rgamma(hh*k, alpha01, beta01), nrow=hh, ncol=k)
H <- HT <- matrix(rgamma(item*k, alpha02, beta02), nrow=k, ncol=item)
WH0 <- W %*% H

#潜在購買行列zを生成
tau <- rbeta(hh, 3, 4.5)   #ベータ分布のパラメータ
Z <- matrix(0, nrow=hh, ncol=item)
for(i in 1:hh){
  Z[i, ] <- rbinom(item, 1, tau[i])
}
mean(tau)

#ポアソン分布から購買行列を生成
WHT <- WH <- WH0 * Z
Data <- matrix(0, nrow=hh, ncol=item)
for(j in 1:item){
  Data[, j] <- rpois(hh, WH[, j])
}


####マルコフ連鎖モンテカルロ法でゼロ過剰NMFを推定####
##ポアソン分布の対数尤度関数
pois <- function(Data, lambda, const){
  LLi <- Data*log(lambda) - lambda - const
  return(LLi)
}

##アルゴリズムの設定
R <- 5000
keep <- 2
disp <- 10
iter <- 0

##事前分布の設定
alpha1 <- 0.01; beta1 <- 0.01
alpha2 <- 0.01; beta2 <- 0.01

##パラメータの真値
W <- WT
H <- HT
r <- rowMeans(Z)

##初期値の設定
W <- matrix(rgamma(hh*k, 0.1, 0.1), nrow=hh, ncol=k)
H <- matrix(rgamma(item*k, 0.1, 0.1), nrow=k, ncol=item)
r <- rowMeans(Data > 0)

##パラメータの格納用配列
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, item, R/keep))
Z_data <- matrix(0, nrow=hh, ncol=item)


##データの設定
z <- as.numeric(Data > 0)
const <- lfactorial(Data)   #ポアソン分布の対数尤度の定数
const_vec <- as.numeric(const)
data_vec <- as.numeric(Data)   #データをベクトル化
u_id <- rep(1:hh, item)
t_id <- rep(1:item, rep(hh, item))

#インデックスを作成
index_zeros <- which(data_vec==0)   #潜在変数の推定対象のベクトル


##対数尤度の基準値
LLbest <- sum(dpois(Data, WT %*% HT, log=TRUE) * Z)
LLst <- sum(dpois(Data, mean(Data), log=TRUE))


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  ##潜在変数Zをサンプリング
  #潜在変数zの割当確率のパラメータ
  r_vec <- r[u_id][index_zeros]   #混合率のベクトル
  Li_zeros <- exp(-as.numeric(W %*% H)[index_zeros])   #データがゼロの時の尤度
  z_rate <- r_vec*Li_zeros / (r_vec*Li_zeros + (1-r_vec))   #潜在変数の割当確率
  
  #二項分布から潜在変数zをサンプリング
  z[index_zeros] <- rbinom(length(index_zeros), 1, z_rate)
  Zi <- matrix(z, nrow=hh, ncol=item)
  r <- rowMeans(Zi)   #混合率を更新
  
  
  ##補助変数lambdaを更新
  lambda <- array(0, dim=c(hh, item, k))
  WH <-  W %*% H
  for(j in 1:k){
    lambda[, , j] <- (W[, j] %*% t(H[j, ]) * Zi / WH) 
  }

  ##ガンマ分布よりWをサンプリング
  weights <- colMeans(Zi)
  for(j in 1:k){
    w1 <- alpha1 + rowSums(lambda[, , j] * Data)
    w2 <- beta1 + sum(H[j, ] * weights)
    W[, j] <- rgamma(hh, w1, w2)   
  }
  W <- W / matrix(colSums(W), nrow=hh, ncol=k, byrow=T) * hh/5   #各列ベクトルを正規化
  
  ##補助変数lambdaを更新
  lambda <- array(0, dim=c(hh, item, k))
  WH <-  W %*% H
  for(j in 1:k){
    lambda[, , j] <- (W[, j] %*% t(H[j, ]) * Zi / WH)
  }
  
  ##ガンマ分布よりHをサンプリング
  weights <- rowMeans(Zi)
  for(j in 1:k){
    h1 <- alpha2 + colSums(lambda[, , j] * Data)
    h2 <- beta2 + sum(W[, j] * weights)
    H[j, ] <- rgamma(item, h1, h2)  
  }
  
  ##サンプリング結果の保存と表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    W_array[, , mkeep] <- W
    H_array[, , mkeep ] <- H
    if(rp > 1000){
      Z_data <- Z_data + Zi
    }
    
    if(rp%%disp==0){
      print(rp)
      print(c(sum(dpois(Data, W %*% H, log=T) * Zi), LLbest, LLst))
      print(round(c(mean(Z), mean(Zi)), 3))
      print(round(cbind(W[1:10, 1:k], WT[1:10, 1:k]), 3))
      print(round(cbind(H[1:k, 1:7], HT[, 1:7]), 3))
    }
  }
}

