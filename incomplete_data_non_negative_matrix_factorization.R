#####欠損のある非負値行列因子分解#####
library(MASS)
library(matrixStats)
library(FAdist)
library(NMF)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

####データの発生####
#データの設定
hh <- 5000   #ユーザー数
item <- 1500  #カテゴリー数
k <- 10   #潜在変数数

##IDの設定
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)


##非負値行列因子分解の仮定に従いデータを生成
#ガンマ分布よりパラメータを設定
alpha01 <- 0.25; beta01 <- 1.0
alpha02 <- 0.15; beta02 <- 0.85
W0 <- matrix(rgamma(hh*k, alpha01, beta01), nrow=hh, ncol=k)
H0 <- matrix(rgamma(item*k, alpha02, beta02), nrow=k, ncol=item)
WH <- W0 %*% H0

#欠損有無のベータ分布のパラメータを設定
beta1 <- rbeta(hh, 9.5, 10.0)   #ユーザ-購買確率
beta2 <- rbeta(item, 7.5, 8.0)   #アイテム購買確率

#ポアソン分布よりデータを生成
Data0 <- matrix(0, nrow=hh, ncol=item)
for(j in 1:item){
  Data0[, j] <- rpois(hh, WH[, j])
}

##欠損がある購買データを生成
#欠損行列を生成
Z0 <- matrix(0, nrow=hh, ncol=item)
for(j in 1:item){
  deficit <- rbinom(hh, 1, beta1 * beta2[j])
  Z0[, j] <- deficit   #欠損を代入
}

#欠損インデックス
Z <- Z0 * Data0 > 0
z_vec <- as.numeric(t(Z))
index_z <- which(z_vec==1)
N <- length(index_z)

#欠損のある購買ベクトル
user_id <- user_id0[index_z]
item_id <- item_id0[index_z]
y <- as.numeric(t(Data0))[index_z]


##ベストなパラメータに対する対数尤度
LLc1 <- sum(dpois(Data0, W0 %*% H0, log=TRUE))   #完全データに対する対数尤度
LLc2 <- sum(as.numeric(t(dpois(Data0, W0 %*% H0, log=TRUE)))[index_z])   #不完全データに対する対数尤度
sparse_data <- as(Data0, "CsparseMatrix")   #スパース行列の作成


####マルコフ連鎖モンテカルロ法でNMFを推定####
##アルゴリズムの設定
R <- 5000
keep <- 4
iter <- 0
disp <- 10

##事前分布の設定
alpha1 <- 0.1; beta1 <- 1
alpha2 <- 0.1; beta2 <- 1

##初期値の設定
W <- matrix(rgamma(hh*k, 0.1, 0.25), nrow=hh, ncol=k)
H <- matrix(rgamma(item*k, 0.1, 0.25), nrow=k, ncol=item)

##サンプリング結果の保存用配列
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, item, R/keep))
lambda <- array(0, dim=c(hh, item, k))

##ユーザーおよびアイテムのインデックスを作成
user_list <- user_vec <- list()
item_list <- item_vec <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
  user_vec[[i]] <- rep(1, length(user_list[[i]]))
}
for(j in 1:item){
  item_list[[j]] <- which(item_id==j)
  item_vec[[j]] <- rep(1, length(item_list[[j]]))
}


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##補助変数lambdaを更新
  lambda <- matrix(0, nrow=N, ncol=k)
  WH <- as.numeric(t(W %*% H))[index_z]
  for(j in 1:k){
    lambda[, j] <- as.numeric(t(W[, j] %*% t(H[j, ])))[index_z] / WH
  }
  
  ##ガンマ分布よりユーザー特徴行列Wをサンプリング
  #ユーザーごとのガンマ分布のパラメータを設定
  W1 <- W2 <- matrix(0, nrow=hh, ncol=k)
  W1_T <- t(lambda * y)   #要素ごとの期待値
  for(i in 1:hh){
    W1[i, ] <- alpha1 + W1_T[, user_list[[i]], drop=FALSE] %*% user_vec[[i]]
    W2[i, ] <- beta1 + H[, item_id[user_list[[i]]], drop=FALSE] %*% user_vec[[i]]
  }
  
  #ガンマ分布よりユーザー特徴行列Wをサンプリング
  W <- matrix(rgamma(hh*k, W1, W2), nrow=hh, ncol=k)
  W <- W / matrix(colSums(W), nrow=hh, ncol=k, byrow=T) * hh/5   #各列ベクトルを正規化
  
  
  ##補助変数lambdaを更新
  lambda <- matrix(0, nrow=N, ncol=k)
  WH <- as.numeric(t(W %*% H))[index_z]
  for(j in 1:k){
    lambda[, j] <- as.numeric(t(W[, j] %*% t(H[j, ])))[index_z] / WH
  }
  
  ##ガンマ分布よりアイテム特徴行列Hをサンプリング
  #アイテムごとのガンマ分布のパラメータを設定
  H1 <- H2 <- matrix(0, nrow=item, ncol=k)
  H1_T <- t(lambda * y)   #要素ごとの期待値
  for(i in 1:item){
    H1[i, ] <- alpha1 + H1_T[, item_list[[i]], drop=FALSE] %*% item_vec[[i]]
    H2[i, ] <- beta1 + t(W[user_id[item_list[[i]]], , drop=FALSE]) %*% item_vec[[i]]
  }
  #ガンマ分布よりユーザー特徴行列Wをサンプリング
  H <- t(matrix(rgamma(item*k, H1, H2), nrow=item, ncol=k))
  
  
  ##サンプリング結果の保存と表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    W_array[, , mkeep] <- W[, 1:k]
    H_array[, , mkeep] <- H[1:k, ]
  }
  
  if(rp%%disp==0){
    print(rp)
    print(c(sum(dpois(y, as.numeric(t(W %*% H))[index_z], log=TRUE)), LLc2))
    print(round(W[1:10, ], 3))
  }
}

####サンプリング結果の確認####
##サンプリング結果をプロット
matplot(t(W_array[100, , ]), type="l")
matplot(t(H_array[, 1, ]), type="l")

##事後平均を計算


##テストデータに対する対数尤度
sum(dpois(as.numeric(t(Data0))[-index_z], as.numeric(t(W %*% H))[-index_z], log=TRUE))
sum(dpois(as.numeric(t(Data0))[-index_z], as.numeric(t(W0 %*% H0))[-index_z], log=TRUE))
round(cbind(as.numeric(t(Data0))[-index_z], as.numeric(t(W %*% H))[-index_z], as.numeric(t(W0 %*% H0))[-index_z]), 3)
round(cbind(y, as.numeric(t(W %*% H))[index_z], as.numeric(t(W0 %*% H0))[index_z]), 3)


