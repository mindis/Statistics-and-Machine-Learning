#####Collaborative Topic Model#####
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

#set.seed(93441)

####データの発生####
##データの設定
k <- 10
hh <- 1500   #レビュー人数
item <- 500   #アイテム数
s <- rpois(hh, rgamma(hh, 6.5, 0.7))   #1人あたりのレビュー数
s[s==0] <- ceiling(runif(length(s[s==0]), 1, 5))
d <- sum(s)   #総文書数
v <- 1000   #語彙数
w <- rpois(d, rgamma(d, 40, 0.8))   #文書あたりの単語数
f <- sum(w)   #総単語数

##アイテム購買を設定
u_list <- list()
par <- as.numeric(extraDistr::rdirichlet(1, rep(3.0, item)))   #アイテム購買確率
for(rp in 1:1000){
  for(i in 1:hh){
    for(j in 1:1000){
      pi <- rmnom(s[i], 1, par)
      if(max(colSums(pi))==1){
        break
      }
    }
    u_list[[i]] <- pi
  }
  if(min(colSums(do.call(rbind, u_list)))==0){
    break
  }
}
U <- do.call(rbind, u_list)
u_vec <- as.numeric(U %*% 1:item)
colSums(U)

##IDを設定
u_id <- rep(1:hh, s)
w_id <- rep(u_vec, w)


##パラメータの設定
#トピックモデルのパラメータ
alpha01 <- rep(0.2, k)
alpha02 <- rep(0.15, v)
theta1 <- thetat1 <- extraDistr::rdirichlet(item, alpha01)
phi <- phit <- extraDistr::rdirichlet(k, alpha02)

#潜在因子モデルのパラメータ
tau1 <- 0.5
tau2 <- 0.025
sigma <- 0.2
theta2 <- thetat2 <- mvrnorm(hh, rep(0, k), diag(tau1, k))
psi <- psit <- t(theta + mvrnorm(item, rep(0, k), diag(tau2, k)))


##Collavorative topic modelのデータを生成
Z_list <- list()
WX <- matrix(0, nrow=d, ncol=v)
word_list <- list()
score <- rep(0, d)

for(i in 1:d){
  #アイテムトピックと単語を生成
  z <- rmnom(w[i], 1, theta1[u_vec[i], ])   #トピックを生成
  z_vec <- as.numeric(z %*% 1:k)
  words <- rmnom(w[i], 1, phi[z_vec, ])   #単語を生成
  words_vec <- colSums(words)
  
  #スコアを生成
  r <- as.numeric(theta2[u_id[i], ] %*% psi[, u_vec[i]]) + rnorm(1, 0, sigma)
  
  #パラメータを格納
  Z_list[[i]] <- z
  WX[i, ] <- words_vec
  word_list[[i]] <- as.numeric(words %*% 1:v)
  score[i] <- r
}

#リストを変換
Z <- do.call(rbind, Z_list)
words_vec <- unlist(word_list)
storage.mode(WX) <- "integer"

