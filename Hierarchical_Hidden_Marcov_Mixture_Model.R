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
alpha12 <- matrix(1.0, nrow=k1, ncol=k1)
diag(alpha12) <- 2.5
alpha21 <- rep(1.5, k2)
alpha22 <- matrix(0.2, nrow=k2, ncol=k2)
diag(alpha22) <- 1.0

#切り替え確率のベータ分布の事前分布を設定
s0 <- 7.5
v0 <- 50.0

#トピック分布のディリクレ分布の事前分布を設定
alpha3 <- array(0.15, dim=c(k2, k3, k1))

#トピックごとの単語分布のディリクレ分布の事前分布を設定
gamma <- matrix(0.003, nrow=k1, ncol=v)
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

##パラメータの真値
beta <- betat
theta11 <- thetat11
theta12 <- thetat12
theta21 <- thetat21
theta22 <- thetat22
theta3 <- thetat3
phi <- phit

##初期値を設定
#多項分布の密度関数の対数尤度の定数
const <- lfactorial(w) - rowSums(lfactorial(WX))

#パラメータの初期値を設定


