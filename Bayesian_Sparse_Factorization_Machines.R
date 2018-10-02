#####Bayesian Sparse Factorization Machines#####
library(MASS)
library(Matrix)
library(matrixStats)
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
item <- 3000   #アイテム数
tag <- 150   #タグ数
pt <- rtpois(item, rgamma(item, 30.0, 0.225), a=1, b=Inf)   #購買接触数
N <- sum(pt)
n <- rtpois(item, 0.8, a=0, b=5)   #タグ数
k <- 10   #基底数
vec_k <- rep(1, k)

#IDを設定
item_id <- rep(1:item, pt)
pt_id <- as.numeric(unlist(tapply(1:N, item_id, rank)))
ID <- data.frame(no=1:N, id=item_id, t=pt_id)   #データの結合
item_list <- list()
for(j in 1:item){
  item_list[[j]] <- which(item_id==j)
}

##階層モデルの説明変数を設定
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

##タグを生成
#パラメータの設定
topic <- 25
omega <- extraDistr::rdirichlet(topic, rep(0.5, tag))
z <- as.numeric(rmnom(item, 1,  extraDistr::rdirichlet(item, rep(1.0, topic))) %*% 1:topic)

#多項分布からタグを生成
max_n <- max(n)
tag_list <- list()
tag_id <- matrix(0, nrow=N, ncol=max_n)
tag_data <- matrix(0, nrow=N, ncol=tag); storage.mode(tag_data) <- "integer"

for(j in 1:item){
  repeat { 
    x <- as.numeric(rmnom(1, n[j], omega[z[j], ]))
    if(max(x)==1){
      tag_list[[j]] <- (x * 1:tag)[x > 0]
      break
    }
  }
  tag_id[item_list[[j]], 1:n[j]] <- matrix(tag_list[[j]], nrow=length(item_list[[j]]), ncol=n[j], byrow=T)
  tag_data[item_list[[j]], ] <- matrix(x, nrow=length(item_list[[j]]), ncol=tag, byrow=T)
}
tag_id0 <- tag_id
tag_id0[tag_id0==0] <- tag+1

#組み合わせを作成
index_combine <- t(combn(c(1, 2, 3, 4, 5), m=2))
combine_list <- list()
combine_n <- rep(0, max(index_combine[, 1]))
for(j in 1:max(index_combine[, 1])){
  combine_list[[j]] <- index_combine[which(index_combine[, 1]==j), 2]
  combine_n[j] <- length(combine_list[[j]])
}

#交差項のインデックスを作成
index_n <- which(rowSums(tag_id > 0) >= 2)
tag_n <- length(index_n)
tag_dt <- sparseMatrix(rep(1:tag_n, max_n), 1:(tag_n*max_n), x=rep(1, tag_n*max_n), dims=c(tag_n, tag_n*max_n))


####応答変数を生成####
rp <- 0
repeat { 
  rp <- rp + 1
  print(rp)
  
  ##モデルのパラメータ
  beta <- betat <- 5.5
  sigma <- sigmat <- 0.75
  
  ##アイテムの階層モデルのパラメータ
  #標準偏差を設定
  tau_v <- tau_vt <- runif(1, 0.3, 0.75)
  
  #回帰係数を設定
  alpha_v <- rep(0, ncol(v))
  for(j in 1:ncol(v)){
    alpha_v[j] <- runif(1, -1.25, 1.25)
  }
  alpha_vt <- alpha_v
  
  #回帰モデルからアイテム個別の変量効果を生成
  theta_v <- theta_vt <- as.numeric(v %*% alpha_v) + rnorm(item, 0, tau_v)
  theta_vec1 <- theta_v[item_id]
  
  ##タグのパラメータを生成
  #タグ個別の変量効果を生成
  tau_r <- tau_rt <- runif(1, 0.25, 0.5)
  theta_r <- theta_rt <- rnorm(tag, 0, tau_r)
  theta_vec2 <- as.numeric(matrix(c(theta_r, 0)[tag_id0], nrow=N, ncol=max_n) %*% rep(1, max_n))
  
  ##交互作用のパラメータを生成
  #特徴ベクトルのパラメータを生成
  tau_g <- tau_gt <- runif(k, 0.01, 0.4) * diag(k)
  theta_g <- theta_gt <- mvrnorm(tag, rep(0, k), tau_g)
  theta_g0 <- rbind(theta_g, 0)

  #交互作用のベクトルを生成
  WH <- rep(0, N)
  for(j in 1:length(combine_n)){
    W <- theta_g0[tag_id0[index_n, j], ]
    H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% theta_g0[tag_id0[index_n, combine_list[[j]]], ])
    WH[index_n] <- WH[index_n] + as.numeric((H * W) %*% vec_k)
  }
  
  ##正規分布から評価ベクトルを生成
  mu <- beta + theta_vec1 + theta_vec2 + WH   #期待値を設定 
  y0 <- rnorm(N, mu, sigma)   #応答変数を生成

  #break条件
  if(min(y0) > -3.0 & max(y0) < 14.0 & mean(y0) > 4.5 & mean(y0) < 6.0){
    break
  }
}

#応答変数を1～10に変換する
y <- round(y0)
y[y > 10] <- 10; y[y < 1] <- 1
hist(y0, breaks=25, col="grey", main="真のスコア分布", xlab="スコア")
hist(y, breaks=25, col="grey", main="切断されたスコア分布", xlab="スコア")


####マルコフ連鎖モンテカルロ法でBayesian SFMを推定####
##アルゴリズムの設定
R <- 2000
burnin <- 500
keep <- 2
disp <- 10
iter <- 0

##事前分布の設定
#回帰パラメータの事前分布
alpha1 <- 0
alpha2 <- rep(0, ncol(v))
alpha3 <- rep(0, k)

#分散の事前分布
s0 <- 0.1
v0 <- 0.1
tau1 <- 100 * diag(ncol(v))
inv_tau1 <- solve(tau1)
tau2 <- 100 * diag(k)
inv_tau2 <- solve(tau2)

##真値の設定
#モデルパラメータ
beta <- betat
sigma <- sigmat

#階層モデルのパラメータ
alpha_v <- alpha_vt
tau_v <- tau_vt
tau_r <- tau_rt
tau_g <- tau_gt

#変量効果のパラメータ
theta_v <- theta_vt
theta_r <- theta_rt
theta_g <- theta_gt[-(tag+1), ]
theta_g0 <- rbind(theta_g, 0)

#評価ベクトルの期待値
theta_vec1 <- theta_v[item_id]
theta_vec2 <- as.numeric(matrix(c(theta_r, 0)[tag_id0], nrow=N, ncol=max_n) %*% rep(1, max_n))
WH <- rep(0, N)
for(j in 1:length(combine_n)){
  W <- theta_g0[tag_id0[index_n, j], ]
  H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% theta_g0[tag_id0[index_n, combine_list[[j]]], ])
  WH[index_n] <- WH[index_n] + as.numeric((H * W) %*% vec_k)
}
mu <- beta + theta_vec1 + theta_vec2 + WH   #期待値を設定 


##初期値の設定
#モデルパラメータ
beta <- mean(y)
sigma <- 1.0

#階層モデルのパラメータ
alpha_v <- rep(0, ncol(v))
tau_v <- 0.1
tau_r <- 0.1
tau_g <- 0.01 * diag(k)

#変量効果のパラメータ
theta_v <- as.numeric(v %*% alpha_v) + rnorm(item, 0, tau_v)
theta_r <- rnorm(tag, 0, tau_r)
theta_g <- mvrnorm(tag, rep(0, k), tau_g)
theta_g0 <- rbind(theta_g, 0)

#評価ベクトルの期待値
theta_vec1 <- theta_v[item_id]
theta_vec2 <- as.numeric(matrix(c(theta_r, 0)[tag_id0], nrow=N, ncol=max_n) %*% rep(1, max_n))
WH <- rep(0, N)
for(j in 1:length(combine_n)){
  W <- theta_g0[tag_id0[index_n, j], ]
  H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% theta_g0[tag_id0[index_n, combine_list[[j]]], ])
  WH[index_n] <- WH[index_n] + as.numeric((H * W) %*% vec_k)
}
mu <- beta + theta_vec1 + theta_vec2 + WH   #期待値を設定 


##パラメータの格納用配列
m <- 0
BETA <- matrix(0, nrow=R/keep, ncol=k)
SIGMA <- rep(0, R/keep)
THETA <- matrix(0, nrow=k-1, ncol=s)

##インデックスとデータの定数を設定
#インデックスを設定
index_list11 <- index_list12 <- list()
index_list21 <- index_list22 <- index_list23 <- list()

for(j in 1:(k-1)){
  #データを抽出
  index <- (allocation_index11==j) + (j_data11==j)
  
  #推定パラメータのインデックス
  j_index <- v_index[rowSums(v_index * cbind(v_index[, 1]==j, v_index[, 2]==j)) > 0, ]
  index_list11[[j]] <- j_index[j_index!=j]
  index_list12[[j]] <- allocation_index21[index==1]
  
  #固定パラメータのインデックス
  index1 <- as.numeric(t(allocation_index11 * (1-index)))
  index2 <- as.numeric(t(allocation_index21 * (1-index)))
  index3 <- as.numeric(t(j_data11 * (1-index)))
  index_list21[[j]] <- index1[index1 > 0]
  index_list22[[j]] <- index2[index2 > 0]
  index_list23[[j]] <- index3[index3 > 0]
}

#データの設定
uu <- t(u) %*% u
inv_uu <- solve(uu + inv_tau1)
v_array <- array(0, dim=c(N, k-2, k-1))
for(j in 1:(k-1)){
  v_array[, , j] <- v[, index_list12[[j]]]
}
v_mu <- rep(0, N)
for(j in 1:(k-1)){
  v_mu <- v_mu + v[, allocation_index21[j, ], drop=FALSE] %*% (theta[j_data11[j, ], ] * theta[allocation_index11[j, ], ]) %*% vec
}

##対数尤度の基準値
LLst <- sum(dnorm(y, mean(y), sd(y), log=TRUE))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##モデルパラメータをサンプリング
  #平均
  er <- as.numeric(y - v_mu)
  
  #多変量正規分布のパラメータを設定
  beta_mu <- inv_uu %*% t(u) %*% er   #多変量正規分布の平均ベクトル
  cov <- sigma^2 * inv_uu   #多変量正規分布の分散
  
  #多変量正規分布から回帰係数をサンプリング
  beta <- mvrnorm(1, beta_mu, cov)
  u_mu <- u %*% beta
  
  
  ##交互作用項の特徴ベクトルをサンプリング
  for(j in 1:(k-1)){
    #応答変数を設定
    er <- as.numeric(y - u_mu - v[, index_list22[[j]]] %*% (theta[index_list21[[j]], ] * theta[index_list23[[j]], ]) %*% vec)
    
    #特徴ベクトルのパラメータを設定
    x <- v_array[, , j] %*% theta[index_list11[[j]], ] 
    inv_xxv <- solve(t(x) %*% x + inv_tau2)
    theta_mu <- inv_xxv %*% t(x) %*% er   #多変量正規分布の平均べクトル
    cov <- sigma^2 * inv_xxv   #多変量正規分布の分散
    
    #多変量正規分布から回帰係数をサンプリング
    theta[j, ] <- mvrnorm(1, theta_mu, cov)
  }
  
  #交互作用の平均ベクトルを更新
  v_mu <- rep(0, N)
  for(j in 1:(k-1)){
    v_mu <- v_mu + v[, allocation_index21[j, ], drop=FALSE] %*% (theta[j_data11[j, ], ] * theta[allocation_index11[j, ], ]) %*% vec
  }
  
  ##モデルの標準偏差をサンプリング
  mu <- u_mu + v_mu
  er <- as.numeric(y - mu)
  
  #逆ガンマ分布のパラメータ
  s1 <- as.numeric(t(er) %*% er) + s0
  v1 <- N + v0
  
  #逆ガンマ分布から標準偏差をサンプリング
  sigma <- sqrt(1/rgamma(1, v1/2, s1/2))
  
  
  ##パラメータの格納とサンプリング結果の表示
  #パラメータを格納
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
    SIGMA[mkeep] <- sigma
    
    if(rp >= burnin){
      m <- m + 1
      THETA <- THETA + theta
    }
  }
  
  if(rp%%disp==0){
    #対数尤度を算出
    LL <- sum(dnorm(y, mu, sigma, log=TRUE))
    
    #サンプリング結果を表示
    print(rp)
    print(c(LL, LLst))
    print(round(rbind(beta, betat), 3))
    print(c(sigma, sigmat))
  }
}

####推定結果の確認と適合度####
##サンプリング結果の可視化
matplot(BETA, type="l", xlab="サンプリング回数", main="betaのサンプリング結果のプロット")
plot(1:length(SIGMA), SIGMA, type="l", xlab="サンプリング回数", main="sigmaのサンプリング結果のプロット")

##パラメータの事後平均
#バーンイン期間
RS1 <- burnin / keep
RS2 <- R / keep

#事後平均を計算  
beta <- colMeans(BETA[RS1:RS2, ])
theta <- THETA / m
sigma <- mean(SIGMA[RS1:RS2])

#パラメータの真値と比較
round(rbind(beta=beta, betat), 3)
round(rbind(theta=theta, thetat), 3)
round(c(sigma, sigmat), 3)


