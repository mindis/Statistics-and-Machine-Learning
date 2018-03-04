#####対応言語トピックモデル#####
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
##文書の設定
k1 <- 15   #名詞のトピック数
k2 <- 15   #動詞と形容詞のトピック数
d <- 2000   #文書数
v1 <- 1000   #名詞の語彙数
v2 <- 500   #動詞(形容詞)の語彙数
s0 <- rpois(d, rgamma(d, 22.5, 0.5))   #動詞と名詞のペア数
s1 <- extraDistr::rtpois(sum(s0), 4.0, a=0, b=Inf)   #ペアごとの名詞数
s2 <- extraDistr::rtpois(sum(s0), 0.7, a=0, b=3)   #ペアごとの動詞数
f0 <- sum(s0)   #総文章数
f1 <- sum(s1)   #名詞の総単語数
f2 <- sum(s2)   #動詞の総単語数


##IDを設定
p_id <- rep(1:d, s0)
w <- as.numeric(tapply(s1, p_id, sum))   #単語数
w_id <- rep(1:d, w)
w_id1 <- rep(1:f0, s1)
w_id2 <- rep(1:f0, s2)


##モデルに基づき単語を生成
#ディレクリ分布のパラメータ
alpha01 <- rep(0.1, k1)   #トピック分布のパラメータ
alpha02 <- rep(0.05, k2)   #混合モデルのパラメータ
beta01 <- rep(0.1, v1)   #名詞の単語分布のパラメータ
beta02 <- rep(0.1, v2)   #動詞の単語分布のパラメータ

#パラメータを生成
theta1 <- thetat1 <- extraDistr::rdirichlet(d, alpha01)
beta <- betat <- matrix(rgamma(k1*k2, 0.4, 0.8), nrow=k2, ncol=k1)
phi <- phit <- extraDistr::rdirichlet(k1, beta01)
psi <- psit <- extraDistr::rdirichlet(k2, beta02)

#逐次的にトピックと単語を生成
WX1 <- matrix(0, nrow=f0, ncol=v1)
WX2 <- matrix(0, nrow=f0, ncol=v2)
word_list1 <- list()
word_list2 <- list()
Z_list1 <- list()
Z_list2 <- list()
Z_sums <- matrix(0, nrow=k2, ncol=k1)
Pr <- matrix(0, nrow=f0, ncol=k2)

for(i in 1:f0){
  if(i%%1000==0){
    print(i)
  }
  #名詞トピックを生成
  pr1 <- theta1[p_id[i], ]
  z1 <- rmnom(s1[i], 1, pr1)
  z1_vec <- as.numeric(z1 %*% 1:k1)

  #動詞トピックを生成
  U <- exp(colSums(z1) %*% beta)
  pr2 <- U / sum(U)
  z2 <- rmnom(s2[i], 1, pr2)
  z2_vec <- as.numeric(z2 %*% 1:k2)
  
  #トピックにもとづき単語を生成
  word1 <- rmnom(s1[i], 1, phi[z1_vec, ])
  word2 <- rmnom(s2[i], 1, psi[z2_vec, ])

  #データを格納
  Z_list1[[i]] <- z1
  Z_list2[[i]] <- z2
  WX1[i, ] <- colSums(word1)
  WX2[i, ] <- colSums(word2)
  word_list1[[i]] <- as.numeric(word1 %*% 1:v1)
  word_list2[[i]] <- as.numeric(word2 %*% 1:v2)
  
  #theta2を推定のためにトピックを格納
  for(j in 1:s2[i]){
    Z_sums[z2_vec[j], ] <- Z_sums[z2_vec[j], ] + colSums(z1)
  }
  Pr[i, ] <- pr2
}

#リストを変換
Z1 <- do.call(rbind, Z_list1)
Z2 <- do.call(rbind, Z_list2)
word_vec1 <- unlist(word_list1)
word_vec2 <- unlist(word_list2)

#動詞のトピック分布を生成
theta2 <- thetat2 <- extraDistr::rdirichlet(k2, Z_sums + 0.01)


####マルコフ連鎖モンテカルロ法で対応言語トピックモデルを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k){
  #負担係数を計算
  Bur <- theta1[w, ] * t(phi)[wd, ]   #尤度
  Br <- Bur / rowSums(Bur)   #負担率
  r <- colSums(Br) / sum(Br)   #混合率
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##パラメータの事前分布
alpha01 <- 10
alpha02 <- 1
beta01 <- 0.1
beta02 <- 0.0001

##パラメータの真値
theta1 <- thetat1
theta2 <- thetat2
phi <- phit
psi <- psit


par <- burden_fr(theta1, phi, word_vec1, w_id, k1)
round(par$Br, 2)
