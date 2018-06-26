#####Tree Structured Mixture Multinomial Model#####
options(warn=0)
library(stringr)
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
library(ggplot2)
#set.seed(2506787)

####データの発生####
##データの設定
#木構造のトピックの設定
m <- 3   #木の最大深さ
k1 <- 3   #1階層目のトピック数
k2 <- rtpois(k1, a=1, b=6, 2.5)   #2階層目のトピック数
k3 <- c()
for(j in 1:k1){
  k3 <- c(k3, rtpois(k2[j], a=0, b=6, 2.0))   #3階層目のトピック数
}
index_k3 <- cbind(c(1, cumsum(k2)[-length(k2)]+1), cumsum(k2))
k <- 1 + sum(c(k1, k2, k3))   #総トピック数
max_k <- max(c(k1, k2, k3))

#文書の設定
d <- 3000   #文書数
w <- rpois(d, rgamma(d, 75, 0.5))   #単語数
f <- sum(w)   #総単語数
v0 <- 200
v1 <- rep(300, k1)
v <- sum(c(v0, v1))   #総語彙数
index_v0 <- 1:v0
index_v1 <- cbind(c(v0+1, (v0+cumsum(v1)[-k1])+1), (v0+cumsum(v1)))

#IDを設定
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}

##パラメータを生成
#単語分布の事前分布
alpha1 <- rep(0.01, v)
alpha21 <- alpha22 <- alpha23 <- rep(0.0025, v)
alpha1[index_v0] <- 5.0
alpha21[index_v1[1, 1]:index_v1[1, 2]] <- 0.1
alpha22[index_v1[2, 1]:index_v1[2, 2]] <- 0.1
alpha23[index_v1[3, 1]:index_v1[3, 2]] <- 0.1
alpha2 <- rbind(alpha21, alpha22, alpha23)

#トピック割当の事前分布
beta1 <- c(30.0, 20.0)   #木の深さの事前分布
beta2 <- rep(5.0, max(c(k1, k2, k3)))   #分岐の事前分布


##すべての単語が出現するまでデータの生成を続ける
for(rp in 1:1000){
  print(rp)
  
  ##木構造のパラメータを生成
  #ベータ分布から停止確率を生成
  gamma1 <- 0.775
  gamma2 <- rbeta(k1, beta1[1], beta1[2])
  gamma3 <- list()
  for(j in 1:k1){
    gamma3[[j]] <- rbeta(k2[j], beta1[1], beta1[2])
  }
  
  #ディリクリ分布から木構造のノード選択確率を生成
  theta1 <- as.numeric(extraDistr::rdirichlet(1, beta2[1:k1]))
  theta2 <- theta3 <- list()
  for(i in 1:k1){
    theta_list <- list()
    theta2[[i]] <- as.numeric(extraDistr::rdirichlet(1, beta2[1:k2[i]]))
    
    for(j in 1:k2[i]){
      x <- k3[(index_k3[i, 1]:index_k3[i, 2])][j]
      if(x==1){
        theta_list[[j]] <- 1
      } else {
        theta_list[[j]] <- as.numeric(extraDistr::rdirichlet(1, beta2[1:x]))
      }
    }
    theta3[[i]] <- theta_list
  }
  
  ##単語分布を生成
  phi0 <- as.numeric(extraDistr::rdirichlet(1, alpha1))
  phi1 <- extraDistr::rdirichlet(k1, alpha2)
  phi2 <- phi3 <- list()
  
  for(i in 1:k1){
    phi_list <- list()
    phi2[[i]] <- extraDistr::rdirichlet(k2[i], alpha2[i, ])
    
    for(j in 1:k2[i]){
      x <- k3[(index_k3[i, 1]:index_k3[i, 2])][j]
      phi_list[[j]] <- extraDistr::rdirichlet(x, alpha2[i, ])
    }
    phi3[[i]] <- phi_list
  }
  
  ##モデルに基づきデータを生成
  #データの格納用配列
  Z1_list <- Z2_list <- word_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  
  for(i in 1:d){
    ##木構造のノードと通過を生成
    #初期ノードの通過を生成
    z11_vec <- rbinom(w[i], 1, gamma1)   
    
    #1階層目のトピックを生成
    z21 <- rmnom(w[i], 1, theta1) * z11_vec
    z21_vec <- as.numeric(z21 %*% 1:k1)
    
    #1階層目のノードの通過を生成
    z12_vec <- rep(0, w[i])
    z12_vec[z11_vec==1] <- rbinom(sum(z11_vec), 1, gamma2[z21_vec])
    
    #2階層目のトピックを生成
    n <- sum(z12_vec)
    x21 <- z21_vec[z12_vec==1]
    theta_par <- matrix(0, nrow=n, ncol=max_k)
    for(j in 1:n){
      theta_par[j, 1:k2[x21[j]]] <- theta2[[x21[j]]]
    }
    z22 <- matrix(0, nrow=w[i], ncol=max_k)
    z22[z12_vec==1, ] <- rmnom(n, 1, theta_par)   #多項分布からトピックを生成
    z22_vec <- as.numeric(z22 %*% 1:max_k)
    
    #2階層目のノードの通過を生成
    gamma_par <- rep(0, n)
    z22_vec1 <- z22_vec[z22_vec > 0]
    for(j in 1:n){
      gamma_par[j] <- gamma3[[x21[j]]][z22_vec1[j]]
    }
    z13_vec <- rep(0, w[i])
    z13_vec[z12_vec==1] <- rbinom(n, 1, gamma_par)
    
    #3階層目のトピックを生成
    n <- sum(z13_vec)
    x22 <- cbind(x21[z13_vec[z12_vec==1]==1], z22_vec[z13_vec==1])
    theta_par <- matrix(0, nrow=n, ncol=max_k)
    
    for(j in 1:n){
      k0 <- k3[(index_k3[x22[j, 1], 1]:index_k3[x22[j, 1], 2])[x22[j, 2]]]
      if(k0==1){
        theta_par[j, 1] <- 1
      } else {
        theta_par[j, 1:k0] <- theta3[[x22[j, 1]]][[x22[j, 2]]]
      }
    }
    z23 <- matrix(0, nrow=w[i], ncol=max_k)
    z23[z13_vec==1, ] <- rmnom(n, 1, theta_par)   #多項分布からトピックを生成
    z23_vec <- as.numeric(z23 %*% 1:max_k)
    
    #トピックを結合
    z2 <- cbind(z0=1, z21_vec, z22_vec, z23_vec)
    depth <- rowSums(z2 > 0)
    
    
    ##単語を生成
    #木構造に応じて単語分布のパラメータを格納
    phi <- matrix(0, nrow=w[i], ncol=v)
    phi[depth==1, ] <- matrix(phi0, nrow=sum(1-z11_vec), ncol=v, byrow=T)
    phi[depth==2, ] <- phi1[z2[depth==2, 2], ]
    
    for(j in 1:w[i]){
      if(depth[j] < 3) next
      if(depth[j]==3){
        phi[j, ]  <- phi2[[z2[j, 2]]][z2[j, 3], ]
      }
      if(depth[j]==4){
        phi[j, ] <- phi3[[z2[j, 2]]][[z2[j, 3]]][z2[j, 4], ]
      }
    }
    
    #多項分布から単語を生成
    word <- rmnom(w[i], 1, phi)
    word_vec <- as.numeric(word %*% 1:v)
    
    
    ##データの格納
    Z1_list[[i]] <- cbind(z11_vec, z12_vec, z13_vec)
    Z2_list[[i]] <- z2
    word_list[[i]] <- word_vec
    WX[i, ] <- colSums(word)
  }
  if(min(colSums(WX)) > 0) next
}
min(colSums(WX))

sum(colSums(WX)==0)

do.call(rbind, Z2_list)
