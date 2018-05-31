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


##欠損ベクトルとIDを設定
#IDを仮設定
item_id0 <- rep(1:item, rep(context, item))
context_id0 <- rep(1:context, item)
n <- length(item_id0)

#要素ごとの出現確率
beta1 <- rbeta(hh, 3.0, 36.0)
beta2 <- rbeta(item, 2.5, 16.0)
beta3 <- rbeta(context, 0.5, 6.0)
beta_vec2 <- beta2[item_id0]
beta_vec3 <- beta3[context_id0]

#ユーザーごとにIDを作成
user_id_list <- item_id_list <- context_id_list <- list()
for(i in 1:hh){
  if(i%%100==0){
    print(i)
  }
  #欠損ベクトルを生成
  prob <- beta1[i] * beta_vec2 * beta_vec3
  deficit <- rbinom(n, 1, prob)
  index_z <- which(deficit==1)
  
  #IDを設定
  user_id_list[[i]] <- rep(i, n)[index_z]
  item_id_list[[i]] <- item_id0[index_z]
  context_id_list[[i]] <- context_id0[index_z]
}
#リストを変換
user_id <- unlist(user_id_list)
item_id <- unlist(item_id_list)
context_id <- unlist(context_id_list)
N <- length(user_id)


#ユーザー×コンテキストのID
uw_index <- paste(user_id, context_id, sep="-")
uw_id <- left_join(data.frame(id=uw_index, no_vec=1:length(uw_index)),
                   data.frame(id=unique(uw_index), no=1:length(unique(uw_index))), by="id")$no

#アイテム×コンテキストのID
vw_index <- paste(item_id, context_id, sep="-")
vw_id <- left_join(data.frame(id=vw_index, no_vec=1:length(vw_index)),
                   data.frame(id=unique(vw_index), no=1:length(unique(vw_index))), by="id")$no

##応答変数が妥当になるまでパラメータの生成を繰り返す
for(rp in 1:1000){
  print(rp)

  ##素性ベクトルを生成
  k1 <- 2; k2 <- 3; k3 <- 4
  x1 <- matrix(runif(N*k1, 0, 1), nrow=N, ncol=k1)
  x2 <- matrix(0, nrow=N, ncol=k2)
  for(j in 1:k2){
    pr <- runif(1, 0.25, 0.55)
    x2[, j] <- rbinom(N, 1, pr)
  }
  x3 <- rmnom(N, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
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
  theta_ut <- theta_u <- as.numeric(u %*% alpha_u + rnorm(hh, 0, tau_u))
  
  
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
  
  #回帰モデルからアイテム個別の回帰パラメータを生成
  theta_vt <- theta_v <- as.numeric(v %*% alpha_v + rnorm(item, 0, tau_v))
  
  
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
  
  #回帰モデルからコンテキスト個別の回帰パラメータを生成
  theta_wt <- theta_w <- w %*% alpha_w + mvrnorm(context, rep(0, kw), tau_w^2 * diag(kw))
  theta_wt1 <- theta_w1 <- theta_w[, 1]
  theta_wt2 <- theta_w2 <- theta_w[, 2]
  
  ##コンテキスト依存のユーザーおよびアイテムバイアスのパラメータを生成
  #コンテキスト依存の変量効果
  tau_uwt <- tau_uw <- 0.4; tau_vwt <- tau_vw <- 0.4
  alpha_uwt <- alpha_uw <- rnorm(unique(uw_id), 0, tau_uw)
  alpha_vwt <- alpha_vw <- rnorm(unique(vw_id), 0, tau_vw)
  
  #コンテキスト依存バイアス平滑化パラメータ
  theta_uwt <- theta_uw <- alpha_uw[uw_id] + theta_u[user_id] * theta_w1[context_id] 
  theta_vwt <- theta_vw <- alpha_vw[vw_id] + theta_v[item_id] * theta_w2[context_id] 
  
  
  ##プロビットモデルから応答変数を生成
  #潜在効用の生成
  mu <- x0 %*% beta + theta_uw + theta_vw
  U <- rnorm(N, mu, 1)
  
  #応答変数を生成
  y <- as.numeric(U > 0)
  if(mean(y) > 0.25 & mean(y) < 0.4) break   #break条件
}

#####モンテカルロEMアルゴリズムでBSHPモデルを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
tol <- 1
iter <- 1
dl <- 100
L <- 500   #モンテカルロサンプリング数

##インデックスを作成
#コンテキスト依存ユーザーインデックス
n_uw <- rep(0, length(unique(uw_id)))
uw_index <- list()
context_u <- rep(0, length(n_uw))
context_w1 <- rep(0, length(n_uw))

for(i in 1:hh){
  if(i%%100==0){
    print(i)
  }
  id <- uw_id[user_index[[i]]]
  min_id <- min(id); max_id <- max(id)
  for(j in min_id:max_id){
    uw_index[[j]] <- user_index[[i]][id==j]
    context_u[j] <- user_id[uw_index[[j]]][1]
    context_w1[j] <- context_id[uw_index[[j]]][1]
    n_uw[j] <- length(uw_index[[j]])
  }
}
N_uw <- length(n_uw)

#コンテキスト依存アイテムインデックス
n_vw <- rep(0, length(unique(vw_id)))
vw_index <- list()
context_v <- rep(0, length(n_vw))
context_w2 <- rep(0, length(n_vw))

for(i in 1:item){
  if(i%%100==0){
    print(i)
  }
  id <- vw_id[item_index[[i]]]
  unique_id <- unique(id)
  for(j in 1:length(unique_id)){
    index <- unique_id[j]
    vw_index[[index]] <- item_index[[i]][id==index]
    context_v[index] <- item_id[vw_index[[index]]][1]
    context_w2[index] <- context_id[vw_index[[index]]][1]
    n_vw[index] <- length(vw_index[[index]])
  }
}
N_vw <- length(n_vw)


##階層モデルのインデックス
#ユーザーインデックス
user_index <- list()
user_n <- rep(0, hh)
for(i in 1:hh){
  user_index[[i]] <- which(context_u==i)
  user_n[i] <- length(user_index[[i]])
}
#アイテムインデックス
item_index <- list()
item_n <- list()
for(j in 1:item){
  item_index[[j]] <- which(context_v==j)
  item_n[j] <- length(item_index[[j]])
}
#コンテキストインデックス
context_index1 <- context_index2 <- list()
context_n1 <- context_n2 <- rep(0, context)
for(j in 1:context){
  context_index1[[j]] <- which(context_w1==j)
  context_index2[[j]] <- which(context_w2==j)
  context_n1[j] <- length(context_index1[[j]])
  context_n2[j] <- length(context_index2[[j]])
}

##初期値の設定

##パラメータの真値
#素性ベクトルの回帰係数
sigma <- 1.0
beta <- betat  

#変量効果のパラメータ
theta_u <- theta_ut   #ユーザーの変量効果
theta_v <- theta_vt   #アイテムの変量効果
theta_w <- theta_wt   #コンテキストの変量効果
alpha_uw <- alpha_uwt   #コンテキスト依存のユーザーバイアス
alpha_vw <- alpha_vwt   #コンテキスト依存のユーザーバイアス
theta_uw <- theta_uwt   #コンテキスト依存のユーザーの変量効果
theta_vw <- theta_vwt   #コンテキスト依存のアイテムの変量効果

#階層モデルのパラメータを生成
alpha_u <- alpha_ut   #ユーザーの階層モデルの回帰係数
tau_v <- tau_vt   #アイテムの階層モデルの標準偏差
alpha_v <- alpha_vt   #アイテムの階層モデルの回帰係数
tau_w <- tau_wt   #コンテキストの階層モデルの標準偏差
alpha_w <- alpha_wt   #コンテキストの階層モデルの標準偏差
tau_uw <- tau_uwt   #コンテキスト依存のユーザーバイアスの標準偏差
tau_vw <- tau_vwt   #コンテキスト依存のユーザーバイアスの標準偏差

#モデルの変量効果
theta_uw <- rep(0, N_uw)
theta_vw <- rep(0, N_vw)
for(i in 1:N_uw){
  theta_uw[i] <- theta_uwt[uw_index[[i]]][1]
}
for(i in 1:N_vw){
  theta_vw[i] <- theta_vwt[vw_index[[i]]][1]
}
theta_vwt0 <- theta_vw; theta_uwt0 <- theta_uw


##切断領域を定義
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)


####モンテカルロEMアルゴリズムをパラメータを推定####
##切断正規分布から潜在効用を生成
theta_uw_vec <- theta_uw[uw_id]
theta_vw_vec <- theta_vw[vw_id]
beta_mu <- as.numeric(x0 %*% beta)
mu <- beta_mu + theta_uw_vec + theta_vw_vec   #潜在効用の期待値
U <- rtnorm(mu, sigma, a, b)   #潜在効用を生成

###モンテカルロEステップで潜在変数をサンプリング
##コンテキスト依存のユーザー変量効果を推定
#データの設定
uw_er <- U - beta_mu - theta_vw_vec   #誤差を設定

#階層モデルのパラメータ
theta_u_vec <- theta_u[context_u] 
theta_w1_vec <- theta_w1[context_w1] 
theta_vec <- theta_u_vec * theta_w1_vec 

#事後分布のパラメータを設定
uw_mu <- rep(0, N_uw)
for(i in 1:N_uw){
  uw_mu[i] <- mean(uw_er[uw_index[[i]]])
}
weights <- tau_uw^2 / (sigma/n_uw + tau_uw^2)    #重み係数
mu_par <- weights*uw_mu + (1-weights)*theta_vec   #事後分布の平均

#正規分布より事後分布をサンプリング
theta_uw <- rnorm(N_uw, mu_par, weights/n_uw)
theta_uw_vec <- theta_uw[uw_id]


##コンテキスト依存のアイテム変量効果を推定
#データの設定
vw_er <- U - beta_mu - theta_uw_vec   #誤差を設定

#階層モデルのパラメータ
theta_v_vec <- theta_v[context_v] 
theta_w2_vec <- theta_w2[context_w2] 
theta_vec <- theta_v_vec * theta_w2_vec 

#事後分布のパラメータを設定
vw_mu <- rep(0, N_vw)
for(i in 1:N_vw){
  vw_mu[i] <- mean(vw_er[vw_index[[i]]])
}
weights <- tau_vw^2 / (sigma/n_vw + tau_vw^2)   #重み係数
mu_par <- weights*vw_mu + (1-weights)*theta_vec   #事後分布の平均

cbind(1/(n_vw/sigma + 1/tau_vw^2), weights/n_vw)

#正規分布より事後分布をサンプリング
theta_vw <- rnorm(N_vw, mu_par, weights/n_vw)
theta_vw_vec <- theta_vw[vw_id]


##ユーザー変量効果を推定
#ユーザーごとに事後分布のパラメータを設定
inv_tau_u <- 1/tau_u^2
mu_par <- rep(0, hh)
sigma_par <- rep(0, hh)

for(i in 1:hh){
  index <- context_w1[user_index[[i]]]
  X <- theta_w1[index]
  Xy <- tau_uw^2 * (t(X) %*% theta_uw[user_index[[i]]])
  XXV <- tau_uw^2 * (t(X) %*% X) + inv_tau_u
  sigma_par[i] <- inv_XXV <- 1/XXV
  mu_par[i] <- inv_XXV %*% (Xy + inv_tau_u %*% theta_u[i])
}
#正規分布から事後分布をサンプリング
rnorm(hh, mu_par, sqrt(sigma_par))

max(sqrt(sigma_par))
t(X) %*% X
X
theta_u[i]

tau_uw^-2

theta_vw[user_index[[i]]]


sd(theta_vw)



