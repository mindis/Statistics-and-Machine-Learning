#####Bias Smoothed Hierarchical Probit Regression model#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)


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
hh <- 5000   #ユーザー数
item <- 2000   #アイテム数
context <- 200   #コンテキスト数
N0 <- hh*item

##IDの設定
#ユーザーとアイテムIDを設定
user_id0 <- rep(1:hh, rep(item, hh))   #ユーザーID
item_id0 <- rep(1:item, hh)   #アイテムID

#コンテキストIDを設定
alpha <- 2.5
G0 <- as.numeric(extraDistr::rdirichlet(1, rep(alpha, context)))*context/(alpha*2)
context_list <- list()
for(i in 1:hh){
  if(i%%100==0){
    print(i)
  }
  prob <- as.numeric(extraDistr::rdirichlet(1, G0))
  context_list[[i]] <- as.numeric(rmnom(item, 1, prob) %*% 1:context)
}
context_id0 <- unlist(context_list)
storage.mode(context_id0) <- "integer"
plyr::count(context_id0)
rm(context_list)

#ユーザー×コンテキストのID
uw_index <- paste(user_id0, context_id0, sep="")
uw_id0 <- left_join(data.frame(id=uw_index, no_vec=1:length(uw_index)),
                    data.frame(id=unique(uw_index), no=1:length(unique(uw_index))), by="id")$no

#アイテム×コンテキストのID
vw_index <- paste(item_id0, context_id0, sep="")
vw_id0 <- left_join(data.frame(id=vw_index, no_vec=1:length(vw_index)),
                    data.frame(id=unique(vw_index), no=1:length(unique(vw_index))), by="id")$no
table(uw_id0[user_id0==1])

##応答変数が妥当になるまでパラメータの生成を繰り返す
for(rp in 1:1000){
  print(rp)

  ##素性ベクトルを生成
  k1 <- 2; k2 <- 3; k3 <- 4
  x1 <- matrix(runif(hh*item*k1, 0, 1), nrow=hh*item, ncol=k1)
  x2 <- matrix(0, nrow=hh*item, ncol=k2)
  for(j in 1:k2){
    pr <- runif(1, 0.25, 0.55)
    x2[, j] <- rbinom(hh*item, 1, pr)
  }
  x3 <- rmnom(hh*item, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
  x0 <- cbind(x1, x2, x3)   #データを結合
  
  
  ##階層モデルの説明変数を生成
  #ユーザーの説明変数
  k1 <- 1; k2 <- 3; k3 <- 5
  u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
  u2 <- matrix(0, nrow=hh, ncol=k2)
  for(j in 1:k2){
    pr <- runif(1, 0.25, 0.55)
    u2[, j] <- rbinom(hh, 1, pr)
  }
  u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
  u <- cbind(1, u1, u2, u3)   #データを結合
  
  #アイテムの説明変数
  k1 <- 2; k2 <- 2; k3 <- 4
  v1 <- matrix(runif(item*k1, 0, 1), nrow=item, ncol=k1)
  v2 <- matrix(0, nrow=item, ncol=k2)
  for(j in 1:k2){
    pr <- runif(1, 0.25, 0.55)
    v2[, j] <- rbinom(item, 1, pr)
  }
  v3 <- rmnom(item, 1, runif(k3, 0.2, 1.25)); v3 <- v3[, -which.min(colSums(v3))]
  v <- cbind(1, v1, v2, v3)   #データを結合
  
  #コンテキストの説明変数
  k1 <- 3; k2 <- 4
  w1 <- matrix(0, nrow=context, ncol=k1)
  for(j in 1:k1){
    pr <- runif(1, 0.25, 0.55)
    w1[, j] <- rbinom(context, 1, pr)
  }
  w2 <- rmnom(context, 1, runif(k2, 0.2, 1.25)); w2 <- w2[, -which.min(colSums(w2))]
  w <- cbind(1, w1, w2)   #データを結合
  
  
  #素性ベクトルの回帰係数を生成
  beta <- rep(0, ncol(x0))
  for(j in 1:ncol(x0)){
    beta[j] <- runif(1, -0.8, 1.2)
  }
  betat <- beta
  
  
  ##階層モデルのパラメータを生成
  ##ユーザーベースの階層モデルのパラメータ
  tau_u <- tau_ut <- 0.5   #標準偏差
  
  #回帰係数を設定
  alpha_u <- rep(0, ncol(u))
  for(j in 1:ncol(u)){
    if(j==1){
      alpha_u[j] <- runif(1, -1.3, -0.5)
    } else {
      alpha_u[j] <- runif(1, -0.4, 0.6)
    }
  }
  alpha_ut <- alpha_u
  
  #回帰モデルからユーザー個別の回帰パラメータを生成
  theta_u <- as.numeric(u %*% alpha_u + rnorm(hh, 0, tau_u))
  
  
  ##アイテムベースの階層モデルのパラメータ
  tau_v <- tau_vt <- 0.7   #標準偏差
  
  #回帰係数を設定
  alpha_v <- rep(0, ncol(v))
  for(j in 1:ncol(v)){
    if(j==1){
      alpha_v[j] <- runif(1, -1.2, -0.3)
    } else {
      alpha_v[j] <- runif(1, -0.6, 0.7)
    }
  }
  alpha_vt <- alpha_v
  
  #回帰モデルからユーザー個別の回帰パラメータを生成
  theta_v <- as.numeric(v %*% alpha_v + rnorm(item, 0, tau_v))
  
  
  ##コンテキストベースの階層モデルのパラメータ
  tau_w <- tau_wt <- 0.4   #標準偏差
  
  #回帰係数を設定
  kw <- 2
  alpha_w <- matrix(0, nrow=ncol(w), ncol=kw)
  for(j in 1:ncol(w)){
    if(j==1){
      alpha_w[j, ] <- runif(kw, 0.2, 0.5)
    } else {
      alpha_w[j, ] <- runif(kw, 0.2, 0.7)
    }
  }
  alpha_wt <- alpha_w
  
  #回帰モデルからユーザー個別の回帰パラメータを生成
  theta_w <- w %*% alpha_w + mvrnorm(context, rep(0, kw), tau_w^2 * diag(kw))
  theta_w1 <- theta_w[, 1]
  theta_w2 <- theta_w[, 2]
  
  ##コンテキスト依存のユーザーおよびアイテムバイアスのパラメータを生成
  #コンテキスト依存の変量効果
  tau_uw <- 0.5; tau_vw <- 0.5
  alpha_uw <- rnorm(unique(uw_id0), 0, tau_uw)
  alpha_vw <- rnorm(unique(vw_id0), 0, tau_vw)
  
  #コンテキスト依存バイアス平滑化パラメータ
  theta_uw <- alpha_uw[uw_id0] + theta_u[user_id0] * theta_w1[context_id0] 
  theta_vw <- alpha_vw[vw_id0] + theta_v[item_id0] * theta_w2[context_id0] 
  
  
  ##プロビットモデルから応答変数を生成
  #潜在効用の生成
  mu <- x0 %*% beta + theta_uw + theta_vw
  U <- rnorm(N0, mu, 1)
  
  #応答変数を生成
  y0 <- as.numeric(U > 0)
  if(mean(y0) > 0.25 & mean(y0) < 0.4) break   #break条件
}

##




##欠損ベクトルを生成
#欠損有無のベータ分布のパラメータを設定
beta1 <- rbeta(hh, 8.5, 10.0)   #ユーザー購買確率
beta2 <- rbeta(item, 6.5, 8.0)   #アイテム購買確率
beta3 <- rbeta(context, 6.0, 7.0)   #コンテキスト購買確率

#インデックスを設定
index_item <- rep(1:item, rep(context, item))
index_context <- rep(1:context, item)
n <- length(index_item)

#ベルヌーイ分布から欠損を生成
Z_list <- list()
for(i in 1:hh){
  if(i%%100==0){
    print(i)
  }
  prob <- beta1[i] * beta2[index_item] * beta3[index_context]
  Z_list[[i]] <- rbinom(n, 1, prob)
}
z_vec <- unlist(Z_list)
mean(z_vec)
sum(z_vec)








