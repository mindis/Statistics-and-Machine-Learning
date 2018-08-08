#####Bayesian Hierarchical Factorization Machines#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(FAdist)
library(bayesm)
library(extraDistr)
library(condMVNorm)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #固有値分解で強制的に正定値行列に修正する
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

####データの発生####
##データの設定
s <- 5   #基底数
hh <- 5000   #ユーザー数
item <- 2000   #アイテム数
N0 <- hh*item

#IDを設定
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)

##欠損ベクトルを生成
#欠損確率を生成
user_prob <- rbeta(hh, 8.5, 10.0)   #ユーザ-購買確率
item_prob <- rbeta(item, 6.5, 8.0)   #アイテム購買確率
prob <- user_prob[user_id0]*item_prob[item_id0]

#ベルヌーイ分布から欠損ベクトルを生成
z_vec <- rbinom(N0, 1, prob)
N <- sum(z_vec)

#欠損ベクトルからidを再構成
user_id <- user_id0[z_vec==1]
item_id <- item_id0[z_vec==1]
rm(user_id0); rm(item_id0); rm(z_vec); rm(prob)
gc(); gc()

#購買数をカウント
freq <- plyr::count(user_id); user_freq <- freq$freq[order(freq$x)]
freq <- plyr::count(item_id); item_freq <- freq$freq[order(freq$x)]
hist(user_freq, col="grey", breaks=25, main="ユーザーごとの購買数", xlab="購買数")
hist(item_freq, col="grey", breaks=25, main="アイテムごとの購買数", xlab="購買数")


##説明変数の生成
#モデルの説明変数を生成
k1 <- 4; k2 <- 5; k3 <- 5
k <- k1 + k2 + k3
x1 <- matrix(0, nrow=N, ncol=k1)
x2 <- matrix(0, nrow=N, ncol=k2)
for(j in 1:k1){
  par <- runif(2, 1.0, 2.5)
  x1[, j] <- rbeta(N, par[1], par[2])
}
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(N, 1, pr)
}
x3 <- rmnom(N, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #データを結合

#ユーザーの説明変数を生成
k1 <- 3; k2 <- 3; k3 <- 4
u1 <- matrix(0, nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k1){
  par <- runif(2, 1.0, 2.5)
  u1[, j] <- rbeta(hh, par[1], par[2])
}
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #データを結合


#アイテムの説明変数を生成
k1 <- 3; k2 <- 2; k3 <- 4
v1 <- matrix(0, nrow=item, ncol=k1)
v2 <- matrix(0, nrow=item, ncol=k2)
for(j in 1:k1){
  par <- runif(2, 1.0, 2.5)
  v1[, j] <- rbeta(item, par[1], par[2])
}
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  v2[, j] <- rbinom(item, 1, pr)
}
v3 <- rmnom(item, 1, runif(k3, 0.2, 1.25)); v3 <- v3[, -which.min(colSums(v3))]
v <- cbind(1, v1, v2, v3)   #データを結合


##交差項を設定
#交差項のインデックスを作成
v_matrix <- matrix(0, nrow=k-1, ncol=k-1)
v_matrix[upper.tri(v_matrix)] <- 1
v_list <- list()
for(j in 1:(k-1)){
  index <- which(v_matrix[j, ]==1)
  if(length(index) > 0){
    v_list[[j]] <- cbind(j1=j, j2=which(v_matrix[j, ]==1))
  }
}
v_index <- do.call(rbind, v_list)



