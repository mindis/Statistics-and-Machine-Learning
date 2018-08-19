#####Joint Geometric Topic Model####
options(warn=0)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(HMM)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####データの発生####
##データの発生
k <- 20   #トピック数
hh <- 5000   #ユーザー数
item <- 1500   #場所数
w <- rtpois(hh, rgamma(item, 25, 0.25), a=1, b=Inf)   #訪問数
f <- sum(w)   #総訪問数

##IDとインデックスの設定
#IDの設定
d_id <- rep(1:hh, w)
t_id <- as.numeric(unlist(tapply(1:f, d_id, rank)))
geo_id01 <- rep(1:hh, rep(item, hh))
geo_id02 <- rep(1:item, hh)

#インデックスの設定
d_index <- geo_index <- list()
for(i in 1:hh){
  d_index[[i]] <- which(d_id==i)
}
for(i in 1:hh){
  geo_index[[i]] <- which(geo_id01==i)
}

##すべてのアイテムが生成されるまで繰り返す
for(rp in 1:1000){
  print(rp)
  
  ##ユーザーとアイテムの経緯度を生成
  #ユーザーの場所集合を生成
  s <- 30 
  rate <- extraDistr::rdirichlet(1, rep(2.0, s))
  point <- as.numeric(rmnom(hh, 1, rate) %*% 1:s)
  
  #経緯度を生成
  longitude <- c(0, 5); latitude <- c(0, 5)
  geo_user0 <- matrix(0, nrow=hh, ncol=2)
  for(j in 1:s){
    index <- which(point==j)
    cov <- runif(2, 0.01, 0.15) * diag(2)
    cov[1, 2] <- cov[2, 1] <- runif(1, -0.6, 0.6) * prod(sqrt(diag(cov)))
    geo_user0[index, ] <- mvrnorm(length(index), c(runif(1, longitude[1], longitude[2]), runif(1, latitude[1], latitude[2])), cov)
  }
  geo_user <- min(geo_user0) + geo_user0
  plot(geo_user, xlab="経度", ylab="緯度", main="ユーザーの場所集合の分布") 
  
  
  #スポットの場所集合を生成
  s <- 25
  rate <- extraDistr::rdirichlet(1, rep(2.0, s))
  point <- as.numeric(rmnom(item, 1, rate) %*% 1:s)
  
  #経緯度を生成
  longitude <- c(0, 5); latitude <- c(0, 5)
  geo_item0 <- matrix(0, nrow=item, ncol=2)
  for(j in 1:s){
    index <- which(point==j)
    cov <- runif(2, 0.005, 0.125) * diag(2)
    cov[1, 2] <- cov[2, 1] <- runif(1, -0.6, 0.6) * prod(sqrt(diag(cov)))
    geo_item0[index, ] <- mvrnorm(length(index), c(runif(1, longitude[1], longitude[2]), runif(1, latitude[1], latitude[2])), cov)
  }
  geo_item <- min(geo_item0) + geo_item0
  plot(geo_item, xlab="経度", ylab="緯度", main="スポットの分布") 
  
  
  #ユーザーと場所のユークリッド距離
  d0 <- sqrt(rowSums((geo_user[geo_id01, ] - geo_item[geo_id02, ])^2))
  hist(d0, breaks=50, xlab="ユークリッド距離", main="2地点間のユークリッド距離の分布", col="grey")
  
  
  ##パラメータを生成
  #トピック分布を生成
  alpha1 <- rep(0.1, k)
  theta <- thetat <- extraDistr::rdirichlet(hh, alpha1)
  
  #場所分布の生成
  alpha2 <- 1.75
  beta <- betat <- 1.0   #バンド幅のパラメータ
  phi <- phit <- mvrnorm(k, rep(0, item), alpha2^2*diag(item))
  
  
  ##応答変数を生成
  Z_list <- V_list <- d_list <- list()
  VX <- matrix(0, nrow=hh, ncol=item); storage.mode(VX) <- "integer"
  
  for(i in 1:hh){
    #トピックを生成
    z <- rmnom(w[i], 1, theta[i, ])
    z_vec <- as.numeric(z %*% 1:k)
    
    #訪問確率を決定
    par <- exp(phi[z_vec, ]) * matrix(exp(-beta/2 * d0[geo_index[[i]]]), nrow=w[i], ncol=item)
    prob <- par / rowSums(par)
    
    #訪問した場所を生成
    v <- rmnom(w[i], 1, prob)
    v_vec <- as.numeric(v %*% 1:item)
     
    #データを格納
    d_list[[i]] <- d0[geo_index[[i]]][v_vec]
    Z_list[[i]] <- z
    V_list[[i]] <- v_vec
    VX[i, ] <- colSums(v)
  }
  #break条件
  if(min(colSums(VX)) > 0) break
}

#リストを変換
d <- unlist(d_list)
Z <- do.call(rbind, Z_list); storage.mode(Z) <- "integer"
v <- unlist(V_list)
sparse_data <- sparseMatrix(i=1:f, j=v, x=rep(1, f), dims=c(f, item))
sparse_data_T <- t(sparse_data)

#データの可視化
plot(geo_user, xlab="経度", ylab="緯度", main="ユーザーの場所集合の分布") 
plot(geo_item, xlab="経度", ylab="緯度", main="スポットの分布")
hist(d0, breaks=50, xlab="ユークリッド距離", main="2地点間のユークリッド距離の分布", col="grey")
