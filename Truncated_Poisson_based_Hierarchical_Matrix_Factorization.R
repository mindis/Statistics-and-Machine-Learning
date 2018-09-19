#####Truncated Poisson based Hierarchical Matrix Factorization#####
options(warn=0)
library(MASS)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####データの発生####
##データの設定
k <- 10   #基底数
hh <- 10000   #ユーザー数
item <- 3000   #アイテム数
pt <- rtpois(hh, rgamma(hh, 27.5, 0.25), a=1, b=Inf)   #購買接触数
hhpt <- sum(pt)
vec <- rep(1, k)

#IDを設定
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
ID <- data.frame(no=1:hhpt, id=user_id, t=pt_id)   #データの結合
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}


##素性ベクトルを生成
k1 <- 2; k2 <- 3; k3 <- 4
x1 <- matrix(runif(hhpt*k1, 0, 1), nrow=hhpt, ncol=k1)
x2 <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #データを結合


##階層モデルの説明変数を設定
#ユーザーの説明変数
k1 <- 3; k2 <- 3; k3 <- 4
u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #データを結合

#アイテムの説明変数
k1 <- 2; k2 <- 3; k3 <- 4
v1 <- matrix(runif(item*k1, 0, 1), nrow=item, ncol=k1)
v2 <- matrix(0, nrow=item, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  v2[, j] <- rbinom(item, 1, pr)
}
v3 <- rmnom(item, 1, runif(k3, 0.2, 1.25)); v3 <- v3[, -which.min(colSums(v3))]
v <- cbind(1, v1, v2, v3)   #データを結合

##アイテムの割当を生成
#セグメント割当を生成
topic <- 25
phi <- extraDistr::rdirichlet(topic, rep(0.5, item))
z <- as.numeric(rmnom(hh, 1,  extraDistr::rdirichlet(hh, rep(2.5, topic))) %*% 1:topic)

#多項分布からアイテムを生成
item_id_list <- list()
for(i in 1:hh){
  if(i%%100==0){
    print(i)
  }
  item_id_list[[i]] <- as.numeric(rmnom(pt[i], 1, phi[z[user_id[user_list[[i]]]], ]) %*% 1:item)
}
item_id <- unlist(item_id_list)
item_list <- list()
for(j in 1:item){
  item_list[[j]] <- which(item_id==j)
}

#個別に和を取るためのスパース行列
user_vec <- sparseMatrix(sort(user_id), unlist(user_list), x=rep(1, hhpt), dims=c(hh, hhpt))
item_vec <- sparseMatrix(sort(item_id), unlist(item_list), x=rep(1, hhpt), dims=c(item, hhpt))

#生成したデータを可視化
freq_item <- plyr::count(item_id); freq_item$x <- as.character(freq_item$x)
hist(freq_item$freq, breaks=25, col="grey", xlab="アイテムの購買頻度", main="アイテムの購買頻度分布")
gc(); gc()


####応答変数を生成####
for(rp in 1:1000){
  print(rp)
  
  ##素性ベクトルのパラメータ
  beta <- betat <- c(-0.7, runif(ncol(x)-1, -1.25, 1.25))
  
  ##階層モデルのパラメータを生成
  #階層モデルの分散パラメータ
  Cov_u <- Cov_ut <- diag(runif(k, 0.01, 0.25), k)   #ユーザー-アイテムの階層モデルの分散
  Cov_v <- Cov_vt <- diag(runif(k, 0.01, 0.25), k)   #アイテムの階層モデルの分散
  
  #階層モデルの回帰係数を設定
  alpha_u <- alpha_ut <- matrix(rnorm(k*ncol(u), 0, 0.75), nrow=ncol(u), ncol=k)
  alpha_v <- alpha_vt <- matrix(rnorm(k*ncol(v), 0, 0.75), nrow=ncol(v), ncol=k)
  
  ##行列分解のパラメータを生成
  theta_u <- theta_ut <- u %*% alpha_u + mvrnorm(hh, rep(0, k), Cov_u)
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, k), Cov_v)
  
  ##ロジットモデルから購買ベクトルを生成
  #ロジットと購買確率を設定
  uv <- uvt <- rowSums(theta_u[user_id, ] * theta_v[item_id, ])
  mu <- as.numeric(x %*% beta)
  logit <- mu + uv
  prob <- exp(logit) / (1 + exp(logit))
  
  #ベルヌーイ分布から購買ベクトルを生成
  y <- rbinom(hhpt, 1, prob)
  print(mean(y))
  if(mean(y) > 0.25 & mean(y) < 0.4) break   #break条件
}

#購買数を確認
sum(y); mean(y)
mean(prob[y==0]); mean(prob[y==1])
hist(prob, main="購買確率の真値の分布", xlab="購買確率", col="grey", breaks=25)

