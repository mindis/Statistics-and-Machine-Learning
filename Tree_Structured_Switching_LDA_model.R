#####Tree_Structured_Switching_LDA_model#####
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
s <- 3   #データタイプ数
k0 <- 10   #一般語のトピック数
k11 <- 3   #ユーザーの上位トピック数
k12 <- rtpois(k11, a=1, b=5, 2.5)   #ユーザーの下位トピックス
k21 <- 3   #アイテムの上位トピック数
k22 <- rtpois(k21, a=1, b=5, 2.5)   #アイテムの下位トピック数
hh <- 2000   #レビュアー数
item <- 500   #アイテム数
v1 <- 300   #評価スコアの語彙数
v2 <- 400   #ユーザートピックの語彙数
v3 <- 400   #アイテムトピックの語彙数
v <- v1 + v2 + v3   #総語彙数
spl <- matrix(1:v1, nrow=s, ncol=v1/s, byrow=T)
v1_index <- 1:v1
v2_index <- (v1+1):(v1+v2)
v3_index <- (v2+1):v

##IDと欠損ベクトルの作成
#IDを仮設定
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)

#欠損ベクトルを作成
for(rp in 1:100){
  m_vec <- rep(0, hh*item)
  for(i in 1:item){
    prob <- runif(1, 0.025, 0.16)
    m_vec[item_id0==i] <- rbinom(hh, 1, prob)
  }
  m_index <- which(m_vec==1)
  
  #完全なIDを設定
  user_id <- user_id0[m_index]
  item_id <- item_id0[m_index]
  d <- length(user_id)   #総レビュー数
  
  #すべてのパターンが生成されればbreak
  if(length(unique(user_id))==hh & length(unique(item_id))==item) break
}

#単語数を設定
w <- rpois(d, rgamma(d, 25, 0.5))   #文書あたりの単語数
f <- sum(w)   #総単語数
n_user <- plyr::count(user_id)$freq
n_item <- plyr::count(item_id)$freq

#単語IDを設定
u_id <- rep(user_id, w)
i_id <- rep(item_id, w)
d_id <- rep(1:d, w)

#インデックスを設定
user_index <- list()
item_index <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
}
