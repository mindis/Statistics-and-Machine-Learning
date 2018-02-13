#####Hidden Marcov Stochastic Block Multinomial model#####
library(MASS)
library(Matrix)
library(bayesm)
library(MCMCpack)
library(gtools)
library(extraDistr)
library(matrixStats)
library(reshape2)
library(qrmtools)
library(slfm)
library(HMM)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

#set.seed(318)

####データの発生####
#データの設定
hh <- 3000   #ユーザー数
pt <- rpois(hh, 20)   #観測期間
max_pt <- max(pt)
hhpt <- sum(pt)   #総レコード数
item <- 2000   #アイテム数
seg_u <- 8   #ユーザーのセグメント数
seg_i <- 7   #アイテムのセグメント数

#IDを設定
u_id <- rep(1:hh, pt)
t_id <- c()
for(i in 1:hh){
  t_id <- c(t_id, 1:pt[i])
}

#インデックスを作成
u_index <- list()
for(i in 1:hh){u_index[[i]] <- which(u_id==i)}

##パラメータとセグメントを生成
#パラメータを生成
alpha01 <- rep(10, seg_u)
alpha02 <- matrix(2.0, nrow=seg_u, ncol=seg_u)
diag(alpha02) <- runif(seg_u, 10, 16)
gamma01 <- extraDistr::rdirichlet(1, alpha01)
gamma02 <- extraDistr::rdirichlet(seg_u, alpha02)

##ユーザーセグメントを生成
z01 <- matrix(0, nrow=hhpt, ncol=seg_u)
z01_vec <- rep(0, hhpt)
for(i in 1:hh){
  for(j in 1:pt[i]){
    if(j==1){
      #1期目のセグメントを生成
      z01[u_index[[i]][j], ] <- rmnom(1, 1, gamma01)
      z01_vec[u_index[[i]][j]] <- as.numeric(z01[u_index[[i]][j], ] %*% 1:seg_u)
      
    } else {
      
      #2期目以降のセグメントを生成
      z01[u_index[[i]][j], ] <- rmnom(1, 1, gamma02[z01_vec[u_index[[i]][j-1]], ])
      z01_vec[u_index[[i]][j]] <- as.numeric(z01[u_index[[i]][j], ] %*% 1:seg_u)
    }
  }
}
z1_cnt <- colSums(z01)

##アイテムセグメントを生成
alpha03 <- rep(20, seg_i)
omega01 <- extraDistr::rdirichlet(1, alpha03)
z02 <- rmnom(item, 1, omega01)
z02_vec <- as.numeric(z02 %*% 1:seg_i)
z2_cnt <- colSums(z02)


##観測モデルのパラメータを生成
#ユーザーセグメント×アイテムセグメントのベータ事前分布のパラメータを発生
hist(rbeta(10000, 0.25, 4.5), col="grey", breaks=25, main="ベータ分布からの乱数", xlab="パラメータ")
theta0 <- matrix(rbeta(seg_u*seg_i, 0.4, 4.5), nrow=seg_u, ncol=seg_i)
round(theta0, 3)


##ベルヌーイ分布から観測行列を発生させる
#ベルヌーイ分布から共起行列を生成
Data <- matrix(0, nrow=hhpt, ncol=item)

for(i in 1:seg_u){
  print(i)
  for(j in 1:seg_i){
    n <- z1_cnt[i] * z2_cnt[j]
    z1_cnt[i]
    Data[z01_vec==i, z02_vec==j] <- matrix(rbinom(n, 1, theta0[i, j]), nrow=z1_cnt[i], ncol=z2_cnt[j])
  }
}

Data_T <- t(Data)
storage.mode(Data) <- "integer"
storage.mode(Data_T) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")
sparse_data_T <- as(Data_T, "CsparseMatrix")  
gc(); gc()

#多項分布のパラメータを設定
thetat <- matrix(0, nrow=seg_u, ncol=seg_i)
for(i in 1:seg_u){
  n <- sum(sparse_data[z01_vec==i, ])
  for(j in 1:seg_i){
    thetat[i, j] <- sum(sparse_data[z01_vec==i, z02_vec==j]) / n
  }
}

phit <- matrix(0, nrow=seg_i, nco=seg_u)
for(i in 1:seg_i){
  n <- sum(sparse_data[, z02_vec==i])
  for(j in 1:seg_u){
    phit[i, j] <- sum(sparse_data[z01_vec==j, z02_vec==i]) / n
  }
}


####マルコフ連鎖モンテカルロ法でHM Sparse Stochastic Block modelを推定####
##アルゴリズムの設定
R <- 5000 
keep <- 2
disp <- 10
iter <- 0
sbeta <- 1.5

##事前分布の設定
alpha1 <- matrix(1, nrow=seg_u, ncol=seg_i)   #ユーザーセグメントのディクレリ事前分布
alpha2 <- matrix(1, nrow=seg_i, ncol=seg_u)   #アイテムセグメントのディクレリ事前分布

##初期値の設定
#混合率の初期値
alpha1 <- rep(20, seg_u)
alpha2 <- matrix(4, nrow=seg_u, ncol=seg_u)
alpha3 <- rep(20, seg_i)
diag(alpha2) <- 20
gamma1 <- extraDistr::rdirichlet(1, alpha1)
gamma2 <- extraDistr::rdirichlet(seg_u, alpha2)
omega1 <- extraDistr::rdirichlet(1, alpha3)

#ブロックごとのパラメータの初期値
index_u <- floor(seq(1, hhpt, length=seg_u+1))
index_i <- floor(seq(1, item, length=seg_i+1))
sortlist1 <- order(rowSums(sparse_data))
sortlist2 <- order(colSums(sparse_data), decreasing=TRUE)

#クラス割当の初期値を設定
z1 <- rep(0, hhpt)
z2 <- rep(0, item)

for(i in 1:(length(index_u)-1)){
  #ユーザーのクラス割当のインデックスを設定
  index1 <- sortlist1[index_u[i]:index_u[i+1]]
  z1[index1] <- i
  
  for(j in 1:(length(index_i)-1)){
    #アイテムのクラス割当のインデックスを設定
    index2 <- sortlist2[index_i[j]:index_i[j+1]]
    z2[index2] <- j
  }
}

#セグメントインデックスを作成
index1 <- list()
index2 <- list()
user_vec <- matrix(0, nrow=hhpt, ncol=seg_u)
item_vec <- matrix(0, nrow=item, ncol=seg_i)
n1 <- c()
n2 <- c()

for(i in 1:seg_u){
  index1[[i]] <- which(z1==i)
  user_vec[index1[[i]], i] <- 1
  n1 <- c(n1, sum(sparse_data[index1[[i]], ]))
}
for(j in 1:seg_i){
  index2[[j]] <- which(z2==j)
  item_vec[index2[[j]], j] <- 1
  n2 <- c(n2, sum(sparse_data[, index2[[j]]]))
}

#パラメータの初期値を設定
#アイテムセグメントの初期パラメータ
oldtheta <- matrix(0, nrow=seg_u, ncol=seg_i) 
for(i in 1:seg_u){
  for(j in 1:(seg_i-1)){
    freq <- sum(sparse_data[index1[[i]], index2[[j]]])
    oldtheta[i, j] <- freq / n1[i]
  }
}
oldtheta[, seg_i] <- 1 - rowSums(oldtheta)
oldtheta <- (oldtheta + 0.00001) / rowSums(oldtheta + 0.00001)

#アイテムセグメントの初期パラメータ
oldphi <- matrix(0, nrow=seg_i, ncol=seg_u)
for(i in 1:seg_i){
  for(j in 1:(seg_u-1)){
    freq <- sum(sparse_data[index1[[j]], index2[[i]]])
    oldphi[i, j] <- freq / n2[i]
  }
}
oldphi[, seg_u] <- 1 - rowSums(oldphi)
oldphi <- (oldphi + 0.00001) / rowSums(oldphi + 0.00001)

#ユーザー、アイテムごとの購買数
n_user <- rowSums(sparse_data)
n_item <- colSums(sparse_data)

##パラメータの格納用配列
THETA <- array(0, dim=c(seg_u, seg_i, R/keep))
PHI <- array(0, dim=c(seg_i, seg_u, R/keep))
SEG1 <- matrix(0, hhpt, ncol=seg_u)
SEG2 <- matrix(0, item, ncol=seg_i)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"
gc(); gc()


##インデックスを作成
max_pt <- max(pt)
index_t11 <- which(t_id==1)
index_t21 <- list()
index_t22 <- list()
n_time <- rep(0, max_pt)
n_time[1] <- length(index_t11)
for(j in 2:max_pt){
  index_t21[[j]] <- which(t_id==j)-1
  index_t22[[j]] <- which(t_id==j)
  n_time[j] <- length(index_t22[[j]])
}


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
##ユーザーごとのセグメント割当を生成
#ユーザーのセグメント割当ごとの尤度を推定
y1 <- as.matrix(sparse_data %*% item_vec)   #アイテムセグメントごとの購入頻度
LLi0 <- y1 %*% t(log(oldtheta))   #セグメントごとの多項分布の対数尤度
LLi_max <- rowMaxs(LLi0)
LLi1 <- exp(LLi0 - LLi_max)   #尤度に変換

#セグメント割当確率の推定とセグメントの生成
z_rate1 <- matrix(0, nrow=hhpt, ncol=seg_u)
Zi1 <- matrix(0, nrow=hhpt, ncol=seg_u)
z1_vec <- rep(0, hhpt)
rf02 <- matrix(0, nrow=seg_u, ncol=seg_u)

for(j in 1:max_pt){
  if(j==1){
    #セグメント割当確率を推定
    n <- n_time[j]
    LLs <- matrix(gamma1, nrow=n, ncol=seg_u, byrow=T) * LLi1[index_t11, ]   #重み付き尤度
    matrix(gamma1, nrow=n, ncol=seg_u, byrow=T)
    z_rate1[index_t11, ] <- LLs / rowSums(LLs)   #割当確率
    
    #多項分布からセグメントを生成
    Zi1[index_t11, ] <- rmnom(n, 1, z_rate1[index_t11, ])
    z1_vec[index_t11] <- as.numeric(Zi1[index_t11, ] %*% 1:seg_u)
    
    #混合率のパラメータを更新
    rf01 <- colSums(Zi1[index_t11, ])
    
  } else {
    
    #セグメントの割当確率
    index <- index_t22[[j]]
    n <- n_time[j]
    LLs <- gamma2[z1_vec[index_t21[[j]]], , drop=FALSE] * LLi1[index, , drop=FALSE]   #重み付き尤度
    z_rate1[index, ] <- LLs / rowSums(LLs)   #割当確率
    
    #多項分布よりセグメントを生成
    Zi1[index, ] <- rmnom(n, 1, z_rate1[index, ])
    z1_vec[index] <- as.numeric(Zi1[index, ] %*% 1:seg_u)
    
    #混合率のパラメータを更新
    rf02 <- rf02 + t(Zi1[index_t21[[j]], , drop=FALSE]) %*% Zi1[index, , drop=FALSE]   #マルコフ推移
  }
}


