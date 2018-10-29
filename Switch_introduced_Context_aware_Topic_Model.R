#####Switch introduced Context aware Topic Model#####
options(warn=0)
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

####データの発生####
##データの設定
L <- 2
k1 <- 15   #同行者トピック数
k2 <- 12   #同行者単語トピック数
k3 <- 15   #単語トピック数
d <- 5000   #文書数
m <- 50   #同行者種類数
v1 <- 500   #同行者依存の語彙数
v2 <- 700   #トピック依存の語彙数
v <- v1 + v2   #総語彙数
index_v1 <- 1:v1   #同行者依存の語彙のインデックス
index_v2 <- (max(index_v1)+1):v   #トピック依存の語彙のインデックス
w <- rpois(d, rgamma(d, 45, 0.3))   #文書あたりの単語数
f <- sum(w)   #総単語数
vec_k1 <- rep(1, k1)
vec_k2 <- rep(1, k2)
vec_k3 <- rep(1, k3)

##文書IDを設定
d_id <- rep(1:d, w)
a_id <- as.numeric(unlist(tapply(1:f, d_id, rank)))

##パラメータの設定
#事前分布のパラメータを設定
alpha1 <- rep(5.0, k1)
alpha2 <- rep(0.2, k2)
alpha3 <- rep(0.15, k3)
beta1 <- rep(0.15, m)
beta2 <- c(rep(0.05, v1), rep(0.001, v2))
beta3 <- c(rep(0.001, v1), rep(0.05, v2))
tau <- c(20.0, 40.0)

##モデルに基づきデータを生成
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  
  ##パラメータを生成
  #ディレクリ分布からパラメータを生成
  theta1 <- thetat1 <- as.numeric(extraDistr::rdirichlet(1, alpha1))
  theta2 <- thetat2 <- extraDistr::rdirichlet(k1, alpha2)
  theta3 <- thetat3 <- extraDistr::rdirichlet(d, alpha3)
  gamma <- gammat <- extraDistr::rdirichlet(k1, beta1)
  phi1 <- extraDistr::rdirichlet(k2, beta2)
  phi2 <- extraDistr::rdirichlet(k3, beta3)

  #単語出現確率が低いトピックを入れ替える
  index_v1 <- which(colMaxs(phi1) < (k3*k3)/f)[which(colMaxs(phi1) < (k3*k3)/f) %in% index_v1]
  index_v2 <- which(colMaxs(phi2) < (k3*k3)/f)[which(colMaxs(phi2) < (k3*k3)/f) %in% index_v2]
  for(j in 1:length(index_v1)){
    phi1[as.numeric(rmnom(1, 1, extraDistr::rdirichlet(1, rep(2.0, k2))) %*% 1:k2), index_v1[j]] <- (k3*k3)/f
  }
  for(j in 1:length(index_v2)){
    phi2[as.numeric(rmnom(1, 1, extraDistr::rdirichlet(1, rep(2.0, k3))) %*% 1:k3), index_v2[j]] <- (k3*k3)/f
  }
  phit1 <- phi1; phit2 <- phi2
  
  #スイッチングパラメータを生成
  pi <- rbeta(d, tau[1], tau[2])
  
  ##文書ごとにトピックと単語を生成
  Z1_list <- Z2_list <- list()
  word_list <- list()
  switch_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  
  #同行者トピックと同行者を生成
  M <- extraDistr::rmnom(d, 1, theta1)
  m_vec <- as.numeric(M %*% 1:k1)
  context <- rmnom(d, 1, gamma[m_vec, ])
  context_vec <- as.numeric(context %*% 1:m)
  
  #単語単位でトピックと単語を生成
  for(i in 1:d){
    
    #スイッチング変数を生成
    s <- rbinom(w[i], 1, pi[i])
    index_s <- which(s==1)
    n <- length(index_s)
    
    #同行者トピックを生成
    z1 <- rmnom(w[i], 1, theta2[m_vec[i], ])
    z1_vec <- as.numeric(z1 %*% 1:k2)
    
    #トピックを生成
    z2 <- rmnom(w[i], 1, theta3[i, ])
    z2_vec <- as.numeric(z2 %*% 1:k3)
    
    
    #スイッチング変数とトピックから単語を生成
    word <- matrix(0, nrow=w[i], ncol=v)
    word[index_s, ] <- rmnom(n, 1, phi1[z1_vec[index_s], ])
    word[-index_s, ] <- rmnom(w[i]-n, 1, phi2[z2_vec[-index_s], ])
    
    #データを格納
    word_list[[i]] <- as.numeric(word %*% 1:v)
    WX[i, ] <- colSums2(word)
    switch_list[[i]] <- s
    Z1_list[[i]] <- z1
    Z2_list[[i]] <- z2
  }
  
  #break条件
  if(min(colSums2(WX)) > 0){
    break
  }
}

#データを変換
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)
storage.mode(WX) <- "integer"
wd <- unlist(word_list)
word_data <- sparseMatrix(1:f, wd, x=rep(1, f), dims=c(f, v))
word_data_T <- t(word_data)
rm(word_list); rm(Z1_list); rm(Z2_list)
gc(); gc()


####マルコフ連鎖モンテカルロ法でSwitch introduced Context aware Topic Modelを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k, vec_k){
  #負担係数を計算
  Bur <- theta[w, ] * t(phi)[wd, ]   #尤度
  Br <- Bur / as.numeric(Bur %*% vec_k)   #負担率
  bval <- list(Br=Br, Bur=Bur)
  return(bval)
}

##アルゴリズムの設定
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10


