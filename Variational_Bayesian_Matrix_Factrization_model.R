#####欠損値のある変分ベイズ行列因子分解#####
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

#set.seed(25897)

####データの発生####
##データの設定
k <- 10   #基底数
hh <- 2000   #ユーザー数
item <- 1000   #アイテム数

##パラメータの設定
sigma <- 1
A <- A_T <- mvrnorm(hh, rep(0, k), diag(1, k))   #ユーザーの特徴行列
B <- B_T <- mvrnorm(item, rep(0, k), diag(1, k))   #アイテムの特徴行列
beta1 <- rbeta(hh, 8.5, 10.0)   #ユーザ-購買確率
beta2 <- rbeta(item, 5.0, 6.0)   #アイテム購買確率


##モデルに基づき応答変数を生成
AB <- A %*% t(B)   #期待値
Y0 <- Y <- matrix(0, nrow=hh, ncol=item)

for(j in 1:item){
  #評価ベクトルを生成
  y_vec <- rnorm(hh, AB[, j], 1)   #正規分布から評価ベクトルを生成
  y_vec0 <- y_vec
  
  #欠損を生成
  deficit <- rbinom(hh, 1, beta1 * beta2[j])
  y_vec[deficit==0] <- NA   #欠損セルにNAを代入
  
  #評価ベクトルを代入
  Y[, j] <- y_vec
  Y0[, j] <- y_vec0
}

##IDと評価ベクトルを設定
N <- length(Y[is.na(Y)==FALSE])
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)

#評価がある要素のみ抽出
index_user <- which(is.na(as.numeric(t(Y)))==FALSE)
user_id <- user_id0[index_user]
item_id <- item_id0[index_user]
y0 <- as.numeric(t(Y))[index_user]

#評価ベクトルを1～5の間の離散値に収める
score_mu <- 3   #平均スコア
y <- as.numeric(round(scale(y0) + score_mu))   #平均3の整数値評価ベクトル
y[y < 1] <- 1
y[y > 5] <- 5


##インデックスの作成
index_user <- list()
index_item <- list()
for(i in 1:hh){
  index_user[[i]] <- which(user_id==i)
}
for(j in 1:item){
  index_item[[j]] <- which(item_id==j)
}

####変分ベイズ法でパラメータを推定####
##事前分布の設定
sigma <- 1
Ca <- diag(1, k)
Cb <- diag(1, k)

##初期値の設定
Cov_A <- array(diag(1, k), dim=c(k, k, hh))
Cov_B <- array(diag(1, k), dim=c(k, k, item))
A <- mvrnorm(hh, rep(0.0, k), diag(0.5, k))
B <- mvrnorm(item, rep(0.0, k), diag(0.5, k))

####変分ベイズ法でパラメータを更新####
##ユーザー特徴行列のパラメータを更新
A <- matrix(0, nrow=hh, ncol=k)
Cov_A <- array(0, dim=c(k, k, hh))
index_item <- item_id[index_user[[i]]]

#分散成分を更新
for(i in 1:hh){
  Cov_sum <- apply(Cov_B[, , index_item], c(1, 2), sum)
  Cov_A[, , i] <- sigma^2 * (solve((t(B[index_item, ]) %*% B[index_item, ] + Cov_sum) + sigma^2 * solve(Ca)))
}
#ユーザーごとに特徴ベクトルを更新
for(i in 1:hh){
  A[i, ] <- sigma^-2 * Cov_A[, , i] %*% colSums(y[index_user[[i]]] * B[item_id[index_user[[i]]], ])
}

##アイテム特徴行列のパラメータを更新
#分散成分を更新
Cov_B <- sigma^2 * solve((t(A) %*% A + Cov_A) + sigma^2 + solve(Cb))

#アイテムごとに特徴ベクトルを更新
B <- matrix(0, nrow=item, ncol=k)
for(j in 1:item){
  B[j, ] <- sigma^-2 * Cov_B %*% colSums(y[index_item[[j]]] * A[user_id[index_item[[j]]], ]) 
}

B[, ]

##ハイパーパラメータを更新

Cov_A
colSums(A^2)
Cov_A


t(A[1, , drop=FALSE]) %*% A[1, ]

