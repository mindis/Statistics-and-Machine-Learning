#####確率的潜在要素解析(PLCAモデル)#####
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

#set.seed(86751)

####データの発生####
k <- 10   #トピック数
d <- 2000   #文書数
v <- 500   #語彙数
f <- 250000   #要素数

##パラメータの設定
#ディレクリ分布のパラメータ
alpha01 <- rep(3.0, k)
alpha11 <- rep(1.0, d)
alpha12 <- rep(0.25, v)

#パラメータを生成
pi <- pit <- extraDistr::rdirichlet(1, alpha01)
theta <- thetat <- extraDistr::rdirichlet(k, alpha11)
phi <- phit <- extraDistr::rdirichlet(k, alpha12)
Mu <- array(0, dim=c(d, v, k))
Par <- matrix(0, nrow=k, ncol=d*v)


##モデルに基づきデータを生成
#要素ごとに多項分布からトピックを生成
Z <- rmnom(f, 1, pi)
z <- as.numeric(Z %*% 1:k) 

#トピックに基づき単語を生成
ID_d <- rep(0, f)
wd <- rep(0, f)

for(j in 1:k){
  print(j)
  index <- which(z==j)
  #文書の割当を生成
  Z1 <- rmnom(length(index), 1, theta[j, ])
  ID_d[index] <- as.numeric(Z1 %*% 1:d)
  
  #単語の割当を生成
  Z2 <- rmnom(length(index), 1, phi[j, ])
  wd[index] <- as.numeric(Z2 %*% 1:v)
}

#文書行列を作成
WX <- matrix(0, nrow=d, ncol=v)
for(i in 1:f){
  WX[ID_d[i], wd[i]] <- WX[ID_d[i], wd[i]] + 1
}


####マルコフ連鎖モンテカルロ法で


##インデックスを作成
doc_list <- list()
word_list <- list()
for(i in 1:d){doc_list[[i]] <- which(ID_d==i)}
for(i in 1:v){word_list[[i]] <- which(wd==i)}



