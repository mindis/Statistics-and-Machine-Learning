#####Switching Multinom LDA model#####
options(warn=0)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
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
r <- 5   #評価スコア数
s <- 3   #極性値数
a <- 3   #分岐数
k11 <- 5   #ユーザーの評価スコアのトピック数
k12 <- 5   #アイテムの評価スコアのトピック数
K1 <- matrix(1:(k11*k12), nrow=k11, ncol=k12, byrow=T)   #トピックの配列
k21 <- 10   #ユーザーのテキストのトピック数
k22 <- 15   #アイテムのテキストのトピック数
hh <- 1000   #レビュアー数
item <- 200   #アイテム数
v1 <- 300   #評価スコアの語彙数
v2 <- 350   #ユーザートピックの語彙数
v3 <- 350   #アイテムトピックの語彙数
v <- v1 + v2 + v3   #総語彙数
spl <- matrix(1:v1, nrow=s, ncol=v1/s, byrow=T)
v1_index <- 1:v1
v2_index <- (v1+1):v2
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

#インデックスを設定
user_index <- list()
item_index <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
}

##パラメータの設定
#ディリクレ分布の事前分布の設定
alpha11 <- rep(0.2, k11)
alpha12 <- rep(0.2, k12)
alpha21 <- rep(0.15, k21)
alpha22 <- rep(0.15, k22)
alpha3 <- c(0.1, 0.225, 0.3, 0.25, 0.125) * r
alpha41 <- c(rep(0.5, v1/s), rep(0.025, v1/s), rep(0.0025, v1/s), rep(0.0025, v2), rep(0.0025, v3))
alpha42 <- c(rep(0.3, v1/s), rep(0.1, v1/s), rep(0.025, v1/s), rep(0.0025, v2), rep(0.0025, v3))
alpha43 <- c(rep(0.2, v1/s), rep(1.0, v1/s), rep(0.2, v1/s), rep(0.0025, v2), rep(0.0025, v3))
alpha44 <- c(rep(0.025, v1/s), rep(0.1, v1/s), rep(0.3, v1/s), rep(0.0025, v2), rep(0.0025, v3))
alpha45 <- c(rep(0.0025, v1/s), rep(0.025, v1/s), rep(0.5, v1/s), rep(0.0025, v2), rep(0.0025, v3))
alpha4 <- rbind(alpha31, alpha32, alpha33, alpha34, alpha35)
alpha51 <- c(rep(0.002, v1/s), rep(0.002, v1/s), rep(0.002, v1/s), rep(0.1, v2), rep(0.002, v3))
alpha52 <- c(rep(0.002, v1/s), rep(0.002, v1/s), rep(0.002, v1/s), rep(0.002, v2), rep(0.1, v3))
beta1 <- c(1.6, 4.8, 5.6)

#事前分布からパラメータを生成
theta11 <- thetat11 <- extraDistr::rdirichlet(hh, alpha11)
theta12 <- thetat12 <- extraDistr::rdirichlet(item, alpha12)
theta21 <- thetat21 <- extraDistr::rdirichlet(hh, alpha21)
theta22 <- thetat22 <- extraDistr::rdirichlet(item, alpha22)
eta <- etat <- extraDistr::rdirichlet(k11*k12, alpha3)
omega <- omegat <- extraDistr::rdirichlet(r, alpha4)
phi <- phit <- extraDistr::rdirichlet(k21, alpha51)
gamma <- gammat <- extraDistr::rdirichlet(k22, alpha52)
lambda <- lambdat <- extraDistr::rdirichlet(hh, beta1)


##モデルに基づきデータを生成
WX <- matrix(0, nrow=d, ncol=v)
y <- rep(0, d)
U1 <- U2 <- list()
Z1_list <- Z21_list <- Z22_list <- wd_list <- list()

i <- 1

#ユーザーとアイテムを抽出
u_index <- user_id[i]
i_index <- item_id[i]

#評価スコアのトピックを生成
u1 <- as.numeric(rmnom(1, 1, theta11[u_index, ]))
u2 <- as.numeric(rmnom(1, 1, theta12[i_index, ]))

#評価スコアのトピックからスコアを生成
y[i] <- as.numeric(rmnom(1, 1, eta[K1[which.max(u1), which.max(u2)], ]) %*% 1:r)


#多項分布からスイッチング変数を生成
z1 <- rmnom(w[i], 1, lambda[u_index, ])
z1_vec <- as.numeric(z1 %*% a)


#ユーザートピックを生成
z21 <- matrix(0, nrow=w[i], ncol=k21)
if(sum(z1[, 2]) > 0){
  z21[z1[, 2]==1, ] <- rmnom(sum(z1[, 2]), 1, theta21[u_index, ])
}
z21_vec <- as.numeric(z21 %*% 1:k21)

#アイテムトピックを生成
z22 <- matrix(0, nrow=w[i], ncol=k22)
if(sum(z1[, 3]) > 0){
  z22[z1[, 3]==1, ] <- rmnom(sum(z1[, 3]), 1, theta22[i_index, ])
}
z22_vec <- as.numeric(z22 %*% 1:k22)

#トピックから単語を生成
words <- matrix(0, nrow=w[i], ncol=v)
if(sum(z1) > 0){
  words[index_z1, ] <- rmnom(sum(z1), 1, phi[z21_vec[index_z1], ])
}
if(sum(1-z1) > 0){
  words[-index_z1, ] <- rmnom(sum(1-z1), 1, gamma[z22_vec[-index_z1], ])
}
word_vec <- as.numeric(words %*% 1:v)
WX[i, ] <- colSums(words)

#データを格納
wd_list[[i]] <- word_vec
Z1_list[[i]] <- z1
Z21_list[[i]] <- z21
Z22_list[[i]] <- z22



