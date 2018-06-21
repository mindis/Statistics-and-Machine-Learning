#####Hierarchical Hidden Marcov Mixture Model#####
options(warn=2)
library(MASS)
library(lda)
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

#set.seed(5723)

####データの発生####
##データの設定
k1 <- 4   #文脈数
k2 <- 7   #状態数
k3 <- 10   #トピック数
d <- 3000   #文書数
s <- rep(300, k1)   #文脈ごとの語彙数
v <- sum(s)   #総語彙数
v1 <- c(1, (cumsum(s)+1)[-k1])
v2 <- cumsum(s)
w <- rpois(d, rgamma(d, 70, 0.5))   #文書ごとの単語数
f <- sum(w)   #総語彙数

##IDの設定
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){t_id <- c(t_id, 1:w[i])}

##パラメータの事前分布の設定
#マルコフ推移確率のディレクリ分布の事前分布を設定
alpha11 <- rep(5.0, k1)
alpha12 <- matrix(1.5, nrow=k1, ncol=k1)
diag(alpha12) <- 10^-300
alpha21 <- rep(1.5, k2)
alpha22 <- matrix(0.2, nrow=k2, ncol=k2)
diag(alpha22) <- 1.0

#切り替え確率のベータ分布の事前分布を設定
s0 <- 7.5
v0 <- 50.0

#トピック分布のディリクレ分布の事前分布を設定
alpha3 <- array(0.15, dim=c(k2, k3, k1))

#トピックごとの単語分布のディリクレ分布の事前分布を設定
gamma <- matrix(0.001, nrow=k1, ncol=v)
for(j in 1:k1){
  gamma[j, v1[j]:v2[j]] <- 0.125
}

##すべての単語が生成されるまでループ
for(rp in 1:1000){
  print(rp)
  
  ##事前分布からパラメータを生成
  #ディリクレ分布からマルコフ推移確率を生成
  theta11 <- thetat11 <- as.numeric(extraDistr::rdirichlet(1, alpha11))
  theta12 <- thetat12 <- extraDistr::rdirichlet(k1, alpha12)
  theta21 <- thetat21 <- extraDistr::rdirichlet(k1, alpha21)
  theta22 <- array(0, dim=c(k2, k2, k1))
  for(j in 1:k1){
    theta22[, , j] <- extraDistr::rdirichlet(k2, alpha22)
  }
  thetat22 <- theta22
  
  #ベータ分布から切り替え変数を生成
  beta <- betat <- rbeta(d, s0, v0)
  
  #ディリクレ分布からトピック分布を生成
  theta3 <- array(0, dim=c(k2, k3, k1))
  for(j in 1:k1){
    theta3[, , j] <- extraDistr::rdirichlet(k2, alpha3[, , j])
  }
  thetat3 <- theta3
  
  #ディリクレ分布から単語分布を生成
  phi <- array(0, dim=c(k3, v, k1))
  for(j in 1:k1){
    phi[, , j] <- extraDistr::rdirichlet(k3, gamma[j, ])
  }
  phit <- phi
  
  ##HHMMMモデルの仮定からデータを生成
  #データの格納用
  z1_list <- list()
  z21_list <- list()
  z22_list <- list()
  z3_list <- list()
  wd_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  
  for(i in 1:d){
    #切換え変数を生成
    z1_vec <- rbinom(w[i], 1, beta[i])
    z1_vec[1] <- 1
    
    ##多項分布よりマルコフ状態推移を生成
    z21_vec <- rep(0, w[i])
    z22_vec <- rep(0, w[i])
    z3_vec <- rep(0, w[i])
    word_vec <- rep(0, w[i])
    words <- matrix(0, nrow=w[i], ncol=v)
    
    for(j in 1:w[i]){
      if(j==1){
        #上位階層の状態数位を生成
        z21 <- rmnom(1, 1, theta11)
        z21_vec[j] <- as.numeric(z21 %*% 1:k1)
        
        #下位階層の状態推移を生成
        z22 <- rmnom(1, 1, theta21[z21_vec[j], ])
        z22_vec[j] <- as.numeric(z22 %*% 1:k2)
        
      } else {
        
        if(z1_vec[j]==1){
          #上位階層の状態推移を生成
          z21 <- rmnom(1, 1, theta12[z21_vec[j-1], ])
          z21_vec[j] <- as.numeric(z21 %*% 1:k1)
          
          #下位階層の状態推移を生成
          z22 <- rmnom(1, 1, theta21[z21_vec[j], ])
          z22_vec[j] <- as.numeric(z22 %*% 1:k2)
          
        } else {
          
          #上位階層の状態推移を生成
          z21_vec[j] <- z21_vec[j-1]
          
          #下位階層の状態推移を生成
          z22 <- rmnom(1, 1, theta22[z22_vec[j-1], , z21_vec[j]])
          z22_vec[j] <- as.numeric(z22 %*% 1:k2)
        }
      }
      ##トピックと単語を生成
      #多項分布からトピックを生成
      z3 <- rmnom(1, 1, theta3[z22_vec[j], , z21_vec[j]])
      z3_vec[j] <- as.numeric(z3 %*% 1:k3)
      
      #トピックから単語を生成
      word <- rmnom(1, 1, phi[z3_vec[j], , z21_vec[j]])
      words[j, ] <- word
      word_vec[j] <- as.numeric(word %*% 1:v)
    }
    
    #生成したデータを格納
    z1_list[[i]] <- z1_vec
    z21_list[[i]] <- z21_vec
    z22_list[[i]] <- z22_vec
    z3_list[[i]] <- z3_vec
    wd_list[[i]] <- word_vec
    WX[i, ] <- colSums(words)
  }
  if(min(colSums(WX)) > 0) break 
}

##データを変換
Z1 <- unlist(z1_list)
Z21 <- unlist(z21_list)
Z22 <- unlist(z22_list)
Z3 <- unlist(z3_list)
wd <- unlist(wd_list)
sparse_data <- sparseMatrix(1:f, wd, dims=c(f, v))
sparse_data_T <- t(sparse_data)


####マルコフ連鎖モンテカルロ法でHHMMMモデルを推定####
##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10


##初期値を設定
#多項分布の密度関数の対数尤度の定数
const <- lfactorial(w) - rowSums(lfactorial(WX))

#パラメータの初期値の事前分布
alpha11 <- rep(5.0, k1)
alpha12 <- matrix(2.5, nrow=k1, ncol=k1)
diag(alpha12) <- 5.0
alpha21 <- rep(5.0, k2)
alpha22 <- matrix(2.5, nrow=k2, ncol=k2)
diag(alpha22) <- 5.0
alpha3 <- rep(5.0, k3)
gamma <- rep(5.0, v)

#パラメータの初期値を生成
beta <- rep(0.2, d)
theta11 <- rep(1/k1, k1)
theta12 <- extraDistr::rdirichlet(k1, alpha12)
theta21 <- extraDistr::rdirichlet(k2, alpha21)
theta22 <- array(0, dim=c(k2, k2, k1))
theta3 <- array(0, dim=c(k2, k3, k1))
phi <- array(0, dim=c(k3, v, k1))
for(j in 1:k1){
  theta22[, , j] <- extraDistr::rdirichlet(k2, alpha22)
  theta3[, , j] <- extraDistr::rdirichlet(k2, alpha3)
  phi[, , j] <- extraDistr::rdirichlet(k3, gamma)
}

##パラメータの真値
beta <- betat
theta11 <- thetat11
theta12 <- thetat12
theta21 <- thetat21
theta22 <- thetat22
theta3 <- thetat3
phi <- phit

##事前分布の設定
#ハイパーパラメータの事前分布
s0 <- 1
v0 <- 1
alpha01 <- 0.1
alpha02 <- 0.1


##MCMC用インデックスを作成
max_word <- max(t_id)
index_t11 <- which(t_id==1)
index_t21 <- index_t22 <- index_t23 <- list()
n_vec <- c(length(index_t11), rep(0, max_word-1))
for(j in 1:max_word){
  index_t21[[j]] <- which(t_id==j)-1
  index_t22[[j]] <- which(t_id==j)
  index <- index_t22[[j]]+1
  index_t23[[j]] <- index[d_id[index_t22[[j]]]==d_id[index]]
  n_vec[j] <- length(index_t22[[j]])
}
vec_k1 <- rep(1, k1)
vec_k2 <- rep(1, k2)
vec_k3 <- rep(1, k3)


#対数尤度の基準値
LLst <- sum(dmnom(WX, w, colSums(WX)/sum(WX), log=TRUE))


####MCMCでパラメータをサンプリング####

##配列の設定
rf21 <- matrix(0, nrow=k1, ncol=k1)
rf22 <- array(0, dim=c(k2, k2, k1))
Zi1 <- rep(0, f)
Zi21 <- matrix(0, nrow=f, ncol=k1)
z21_vec <- rep(0, f)
Zi22 <-  matrix(0, nrow=f, ncol=k2)
z22_vec <- rep(0, f)


##トピックごとの期待尤度を推定
LLt <- array(0, dim=c(f, k2, k1))
phi_wd <- array(0, dim=c(f, k3, k1)) 
for(j in 1:k1){
  phi_wd[, , j] <- t(phi[, , j])[wd, ]
  for(k in 1:k2){
    LLt[, k, j] <- (matrix(theta3[k, , j], nrow=f, ncol=k3, byrow=T) * phi_wd[, , j]) %*% vec_k3 
  }
}

for(pd in 1:max_word){
  if(pd==1){
    ##1単語目の潜在状態を生成
    ##上位階層の潜在状態を生成
    #対数尤度を上位階層のの期待尤度に変換
    LLt21 <- LLt[index_t11, , ]
    n <- length(index_t11)
    par <- matrix(0, nrow=n[pd], ncol=k1)
    for(j in 1:k1){
      par[, j] <- (matrix(theta21[j, ], nrow=n, ncol=k2, byrow=T) * LLt21[, , j]) %*% vec_k2
    }
    
    #多項分布から上位階層の状態を生成
    r <- matrix(theta11, nrow=n, ncol=k1, byrow=T)
    par_rate <- r*par / rowSums(r*par)   #潜在変数の割当確率
    Zi21[index_t11, ] <- rmnom(n, 1, par_rate)   #多項分布から状態を生成
    z21_vec[index_t11] <- as.numeric(Zi21[index_t11, ] %*% 1:k1)
    
    
    ##下位階層の潜在状態を生成
    #上位階層の状態に応じて尤度計算
    zi21 <- Zi21[index_t11, ]
    par <- matrix(0, nrow=n, ncol=k2)
    for(j in 1:k1){
      par <- par + LLt21[, , j] * zi21[, j]
    }
    
    #多項分布から下位階層の状態を生成
    r <- theta21[z21_vec[index_t11], ]
    par_rate <- r*par / (rowSums(r*par))   #潜在変数の割当確率
    Zi22[index_t11, ] <- rmnom(n, 1, par_rate)   #多項分布から状態を生成
    z22_vec[index_t11] <- as.numeric(Zi22[index_t11, ] %*% 1:k2)
    
    } else {
      
    ##2単語目以降のマルコフ推移に基づき潜在状態を生成
    ##潜在状態が切り替わるかどうかを生成
    #データの設定
    index1 <- index_t21[[pd]]
    index2 <- index_t22[[pd]]
    n <- length(index2)
    zi21_j <- Zi21[index1, , drop=FALSE]
    zi22_j <- Zi22[index1, , drop=FALSE]
    
    #期待尤度を計算
    LLt21 <- LLt[index2, , , drop=FALSE]
    LLt22 <- matrix(0, nrow=n, ncol=k2)
    par <- matrix(0, nrow=n, ncol=k1)
    for(j in 1:k1){
      theta22_r <- theta22[zi22_j, , j]
      par[, j] <- (LLt21[, , j]*theta22_r) %*% vec_k2
      LLt22 <- LLt22 + LLt21[, , j]*zi21_j[, j]*theta22_r
    }
    
    #上位階層の潜在状態切換え変数を生成
    par1 <- par0 <- rep(0, n)
    for(j in 1:k1){
      par1 <- par1 + par[, j]*zi21_j[, j]
      par0 <- par0 + matrix(theta12[j, -j], nrow=n, ncol=k1-1, byrow=T)*par[, -j]*zi21_j[, j]
    }
    r <- beta[d_id[index2]]   #混合率
    beta_par0 <- r * rowSums(par0)
    beta_rate <- beta_par0 / ((1-r)*par1+beta_par0)   #切換え変数の割当確率
    Zi1[index2] <- rbinom(n, 1, beta_rate)   #ベルヌーイ分布から切換え変数を生成
    index_z1 <- which(Zi1[index2]==1)
    
    
    ##切換え変数を条件として上位階層を生成
    if(length(index_z1)==0){
      Zi21[index2, ] <- Zi21[index1, ]
      z21_vec[index2] <- z21_vec[index1]
    } else {
      r <- theta12[z21_vec[index1], ]
      Zi21[index2, ] <- Zi21[index1, ]; z21_vec[index2] <- z21_vec[index1]   #1期前の潜在状態を繰り返す
      par_rate <- r*par / rowSums(r*par)   #潜在状態の割当確率
      
      if(nrow(par_rate)==1){
        Zi21[index2, ] <- rmnom(length(index_z1), 1, par_rate[index_z1, ])   #多項分布から潜在状態を生成
        z21_vec[index2] <- as.numeric(Zi21[index2, ] %*% 1:k1)
      } else {
        Zi21[index2, ][index_z1, ] <- rmnom(length(index_z1), 1, par_rate[index_z1, ])   #多項分布から潜在状態を生成
        z21_vec[index2][index_z1] <- as.numeric(Zi21[index2, ][index_z1, ] %*% 1:k1)
      }
      #上位階層のマルコフ推移行列を更新
      rf21 <- rf21 + t(Zi21[index1, , drop=FALSE]) %*% (Zi21[index2, , drop=FALSE] * Zi1[index2])
    }
    
    ##下位階層の潜在状態を生成
    #上位階層の状態に応じて尤度計算
    par <- matrix(0, nrow=n, ncol=k2)
    for(j in 1:k1){
      index <- which(Zi21[index2, j]==1)
      if(length(index)==0) next
      LLt22 <- LLt[index2, , , drop=FALSE]
      par[index, ] <- theta22[z22_vec[index1][index], , j] * LLt22[, , j, drop=FALSE][index, , 1]
    }
    
    #多項分布から下位階層の状態を生成
    par_rate <- par / rowSums(par)   #潜在変数の割当確率
    Zi22[index2, ] <- rmnom(n, 1, par_rate)   #多項分布から潜在状態を生成 
    z22_vec[index2] <- as.numeric(Zi22[index2, ] %*% 1:k2)
    
    #マルコフ推移行列を更新
    for(j in 1:k1){
    rf22[, , j] <- rf22[, , j] + t(Zi22[index1, , drop=FALSE]) %*% 
      (Zi22[index2, , drop=FALSE] * Zi21[index2, j] * (1-Zi1[index2]))
    }
  }
}

rf22
rf21
