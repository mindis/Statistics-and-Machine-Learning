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
s1 <- 10; vec1 <- rep(1, s1)   #行列分解の基底数
s2 <- 5; vec2 <- rep(1, s2)   #交互作用の基底数
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
z_matrix <- matrix(0, nrow=k-1, ncol=k-1)
z_matrix[upper.tri(z_matrix)] <- 1
z_list <- list()
for(j in 1:(k-1)){
  index <- which(z_matrix[j, ]==1)
  if(length(index) > 0){
    z_list[[j]] <- cbind(j1=j, j2=which(z_matrix[j, ]==1))
  }
}
z_index <- do.call(rbind, z_list)

#交差項の入力変数を作成
z <- x[, z_index[, 1]+1] * x[, z_index[, 2]+1]

#説明変数の割当インデックス
allocation_index11 <- allocation_index12 <- matrix(0, nrow=k-1, ncol=k-2)
allocation_index21 <- allocation_index22 <- matrix(0, nrow=k-1, ncol=k-2)
for(j in 1:(k-1)){
  index <- which(rowSums(z_index==j)==1)
  allocation_index11[j, ] <- allocation_index12[j, ] <- 
    rowSums(matrix(as.numeric(z_index[index, ]!=j), nrow=length(index)) * z_index[index, ])
  allocation_index21[j, ] <- allocation_index22[j, ] <-  index
}
allocation_index11[lower.tri(allocation_index11)] <- 0
allocation_index21[lower.tri(allocation_index21)] <- 0
vec <- rep(1, s)
j_data12 <- matrix(1:(k-1), nrow=k-1, ncol=k-2) 
j_data11 <- j_data12 * (allocation_index11 > 0)


####応答変数を生成####
rp <- 0
repeat { 
  print(rp <- rp + 1)
  
  ##パラメータと応答変数を生成
  #モデルの標準偏差
  sigma <- sigmat <- 1.0
  
  #階層モデルの分散パラメータ
  Cov_x <- Cov_xt <- runif(ncol(x), 0.05, 0.25) * diag(ncol(x))
  Cov_u <- Cov_ut <- runif(s1, 0.05, 0.20) * diag(s1)
  Cov_v <- Cov_vt <- runif(s1, 0.05, 0.20) * diag(s1)
  Cov_z <- Cov_zt <- runif(s2, 0.05, 0.15) * diag(s2)
  
  #階層モデルの回帰係数を設定
  alpha_x <- matrix(0, nrow=ncol(u), ncol=ncol(x))
  alpha_u <- matrix(0, nrow=ncol(u), ncol=s1)
  alpha_v <- matrix(0, nrow=ncol(v), ncol=s1)
  alpha_z <- array(0, dim=c(ncol(u), s2, k-1))
  
  for(j in 1:ncol(u)){
    if(j==1){
      alpha_x[j, ] <- runif(ncol(x), -0.4, 0.3)
      alpha_u[j, ] <- runif(s1, -0.4, 0.2)
      alpha_z[j, , ] <- matrix(rnorm(s2*(k-1), 0, 0.2), nrow=s2, ncol=k-1)
    } else {
      alpha_x[j, ] <- runif(ncol(x), -0.4, 0.3)
      alpha_u[j, ] <- runif(s1, -0.35, 0.2)
      alpha_z[j, , ] <- matrix(rnorm(s2*(k-1), 0, 0.2), nrow=s2, ncol=k-1)
    }
  }
  for(j in 1:ncol(v)){
    if(j==1){
      alpha_v[j, ] <- runif(s2, -0.5, 0.4)
    } else {
      alpha_v[j, ] <- runif(s2, -0.4, 0.4)
    }
  }
  alpha_xt <- alpha_x; alpha_ut <- alpha_u; alpha_vt <- alpha_v; alpha_zt <- alpha_z   #真値を格納
  
  #多変量回帰モデルからユーザー個別の回帰パラメータを生成
  theta_x <- theta_xt <- u %*% alpha_x + mvrnorm(hh, rep(0, ncol(x)), Cov_x)   #変量効果のパラメータ
  theta_u <- theta_ut <- u %*% alpha_u + mvrnorm(hh, rep(0, s1), Cov_u)   #ユーザーの行列分解のパラメータ
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, s1), Cov_v)   #アイテムの行列分解のパラメータ
  theta_z <- array(0, c(hh, s2, k-1))
  for(j in 1:(k-1)){
    theta_z[, , j] <- u %*% alpha_z[, , j] + mvrnorm(hh, rep(0, s2), Cov_z)   #交互作用のパラメータ
  }
  theta_xt <- theta_x; theta_ut <- theta_u; theta_vt <- theta_v; theta_zt <- theta_z
  
  ##正規分布から効用と購買ベクトルを生成
  #変量効果のパラメータ
  x_mu <- as.numeric((x * theta_x[user_id, ]) %*% rep(1, ncol(x)))
  
  #行列分解のパラメータ
  uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)

  #交互作用のパラメータ
  z_vec <- rep(0, N)
  j_data_vec11 <- as.numeric(j_data11)[as.numeric(j_data11) > 0]; m <- length(j_data_vec11)
  allocation_vec11 <- as.numeric(allocation_index11)[as.numeric(allocation_index11) > 0]
  allocation_vec21 <- as.numeric(allocation_index21)[as.numeric(allocation_index21) > 0]
  for(j in 1:length(j_data_vec11)){
    z_vec <- z_vec + z[, allocation_vec21[j]] * (theta_z[user_id, , j_data_vec11[j]] * theta_z[user_id, , allocation_vec11[j]]) %*% vec2
  }
  z_vec <- as.numeric(z_vec)
  
  #潜在効用と応答変数を生成
  mu <- x_mu + z_vec + uv   #期待値
  U <- mu + rnorm(N, 0, sigma)   #潜在効用を生成
  
  #購買ベクトルに変換
  y <- ifelse(U > 0, 1, 0)
  if(mean(y) > 0.25 & mean(y) < 0.4) break   #break条件
}

#生成した応答変数を確認
mean(y)   #購買確率
prob <- pnorm(U, 0, sigma)   #応答確率
mean(prob[y==1]); mean(prob[y==0])   #購買有無別の応答確率
hist(U, col="grey", main="潜在効用の分布", xlab="潜在効用", breaks=25)


####マルコフ連鎖モンテカルロ法で階層ベイズFactorization Machinesを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  par <- qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma)
  return(par)
}

##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
R <- 2000
keep <- 2  
iter <- 0
burnin <- 500/keep
disp <- 10

##インデックスを設定
user_index <- item_index <- list()
ui_id <- iu_id <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
  ui_id[[i]] <- item_id[user_index[[i]]]
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
  iu_id[[j]] <- user_id[item_index[[j]]]
}

##事前分布の設定
#変量効果の階層モデルの事前分布
Deltabar1 <- matrix(0, nrow=ncol(u), ncol=k)   #階層モデルの回帰係数の事前分布の平均
ADelta1 <- 0.01 * diag(1, ncol(u))   #階層モデルの回帰係数の事前分布の分散
nu1 <- k + 1   #逆ウィシャート分布の自由度
V1 <- nu1 * diag(rep(1, k)) #逆ウィシャート分布のパラメータ

#ユーザーの行列分解の階層モデルの事前分布
Deltabar2 <- matrix(0, nrow=ncol(u), ncol=s1)   #階層モデルの回帰係数の事前分布の平均
ADelta2 <- 0.01 * diag(1, ncol(u))   #階層モデルの回帰係数の事前分布の分散
nu2 <- s1 + 1   #逆ウィシャート分布の自由度
V2 <- nu2 * diag(rep(1, s1)) #逆ウィシャート分布のパラメータ

#アイテムの行列分解の階層モデルの事前分布
Deltabar3 <- matrix(0, nrow=ncol(v), ncol=s1)   #階層モデルの回帰係数の事前分布の平均
ADelta3 <- 0.01 * diag(1, ncol(v))   #階層モデルの回帰係数の事前分布の分散
nu3 <- s1 + 1   #逆ウィシャート分布の自由度
V3 <- nu3 * diag(rep(1, s1)) #逆ウィシャート分布のパラメータ

#交互作用項の階層モデルの事前分布
Deltabar4 <- matrix(0, nrow=ncol(u), ncol=s2)   #階層モデルの回帰係数の事前分布の平均
ADelta4 <- 0.01 * diag(1, ncol(u))   #階層モデルの回帰係数の事前分布の分散
nu4 <- s2 + 1   #逆ウィシャート分布の自由度
V4 <- nu4 * diag(rep(1, s2)) #逆ウィシャート分布のパラメータ

#変量効果の階層モデルの事前分布
Deltabar1 <- matrix(0, nrow=ncol(u), ncol=k)   #階層モデルの回帰係数の事前分布の平均
ADelta1 <- 0.01 * diag(1, ncol(u))   #階層モデルの回帰係数の事前分布の分散
nu1 <- k + 1   #逆ウィシャート分布の自由度
V1 <- nu1 * diag(rep(1, k)) #逆ウィシャート分布のパラメータ


##パラメータの真値
#階層モデルのパラメータの真値
Cov_x <- Cov_xt; Cov_u <- Cov_ut; Cov_v <- Cov_vt; Cov_z <- array(Cov_zt, dim=c(s2, s2, k-1))
inv_Cov_x <- solve(Cov_x); inv_Cov_u <- solve(Cov_u); inv_Cov_v <- solve(Cov_v)
inv_Cov_z <- array(0, dim=c(s2, s2, k-1))
for(j in 1:(k-1)){
  inv_Cov_z[, , j] <- solve(Cov_z[, , j])
}
alpha_x <- alpha_xt; alpha_u <- alpha_ut; alpha_v <- alpha_vt; alpha_z <- alpha_zt
x_mu <- u %*% alpha_x; u_mu <- u %*% alpha_u; v_mu <- v %*% alpha_v
z_mu <- array(0, c(hh, s2, k-1))
for(j in 1:(k-1)){
  z_mu[, , j] <- u %*% alpha_z[, , j]
}

#モデルパラメータの真値
sigma <- sigmat
theta_x <- theta_xt; theta_u <- theta_ut; theta_v <- theta_vt; theta_z <- theta_zt
user_mu <- as.numeric((x * theta_x[user_id, ]) %*% rep(1, ncol(x)))

#行列分解と交互作用項のパラメータ
uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)
z_vec <- rep(0, N)
for(j in 1:length(j_data_vec11)){
  z_vec <- z_vec + z[, allocation_vec21[j]] * (theta_z[user_id, , j_data_vec11[j]] * theta_z[user_id, , allocation_vec11[j]]) %*% vec2
}
z_vec <- as.numeric(z_vec)


##パラメータの初期値
#階層モデルの初期値
Cov_x <- diag(0.01, ncol(x)); Cov_u <- Cov_v <- diag(0.01, s1); Cov_z <- array(diag(0.01, s2), dim=c(s2, s2, k-1))
alpha_x <- matrix(0, nrow=ncol(u), ncol=ncol(x)); x_mu <- u %*% alpha_x
alpha_u <- matrix(0, nrow=ncol(u), ncol=s1); u_mu <- u %*% alpha_u
alpha_v <- matrix(0, nrow=ncol(v), ncol=s1); v_mu <- v %*% alpha_v
alpha_z <- array(0, dim=c(ncol(u), s2, k-1)); z_mu <- array(0, dim=c(hh, s2, k-1))

#モデルパラメータの初期値
sigma <- 1
theta_x <- mvrnorm(hh, as.numeric(solve(t(x) %*% x) %*% t(x) %*% y), Cov_x) 
theta_u <- mvrnorm(hh, rep(0, s1), Cov_u)
theta_v <- mvrnorm(item, rep(0, s1), Cov_v)
theta_z <- array(0, dim=c(hh, s2, k-1))
for(j in 1:(k-1)){
  theta_z[, , j] <- mvrnorm(hh, rep(0, s2), Cov_z[, , j])
}

#行列分解と交互作用項のパラメータ
uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)
z_vec <- rep(0, N)
for(j in 1:m){
  z_vec <- z_vec + z[, allocation_vec21[j]] * (theta_z[user_id, , j_data_vec11[j]] * theta_z[user_id, , allocation_vec11[j]]) %*% vec2
}
z_vec <- as.numeric(z_vec)


##インデックスとデータの定数を設定
#インデックスを設定
index_list11 <- index_list12 <- list()
index_list21 <- index_list22 <- index_list23 <- list()

for(j in 1:(k-1)){
  #データを抽出
  index <- (allocation_index11==j) + (j_data11==j)
  
  #推定パラメータのインデックス
  j_index <- z_index[rowSums(z_index * cbind(z_index[, 1]==j, z_index[, 2]==j)) > 0, ]
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

#データの定数を設定
xx_list <- list()
for(i in 1:hh){
  xx_list[[i]] <- t(x[user_index[[i]], ]) %*% x[user_index[[i]], ]
}
z_list <- list()
for(i in 1:hh){
  z_array <- array(0, dim=c(length(user_index[[i]]), k-2, k-1))
  for(j in 1:(k-1)){
    z_array[, , j] <- z[user_index[[i]], index_list12[[j]]]
  }
  z_list[[i]] <- z_array
}

##切断領域を定義
index_y1 <- which(y==1)
index_y0 <- which(y==0)
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)

##対数尤度の基準値
prob <- mean(y)
LLst <- sum(y*log(prob)) + sum((1-y)*log(1-prob))   #対数尤度


####ギブスサンプリングでパラメータをサンプリング####
  for(rp in 1:R){
    
  ##切断正規分布より潜在効用を生成
  mu <- user_mu + uv + z_vec   #潜在効用の期待値
  U <- extraDistr::rtnorm(N, mu, sigma, a, b)   #潜在効用を生成
  
  ##ユーザーの回帰ベクトルをサンプリング
  #モデルの応答変数
  y_er <- U - uv - z_vec
  
  for(i in 1:hh){
    #回帰ベクトルの事後分布のパラメータ
    XX <- xx_list[[i]]
    Xy <- t(x[user_index[[i]], ]) %*% y_er[user_index[[i]]]
    inv_XXV <- solve(XX + inv_Cov_x)
    mu <- inv_XXV %*% (Xy + inv_Cov_x %*% x_mu[i, ])   #事後分布の平均
    
    #多変量正規分布から回帰ベクトルをサンプリング
    theta_x[i, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  user_mu <- as.numeric((x * theta_x[user_id, ]) %*% rep(1, ncol(x)))
  
  
  ##ユーザーの特徴行列をサンプリング
  #モデルの応答変数
  u_er <- U - user_mu - z_vec
  
  for(i in 1:hh){
    #特徴ベクトルの事後分布のパラメータ
    X <- theta_v[ui_id[[i]], ]
    Xy <- t(X) %*% u_er[user_index[[i]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_u)
    mu <- inv_XXV %*% (Xy + inv_Cov_u %*% u_mu[i, ])
    
    #多変量正規分布から特徴ベクトルをサンプリング
    theta_u[i, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  
  ##アイテムの特徴行列をサンプリング
  for(j in 1:item){
    
    #特徴ベクトルの事後分布のパラメータ
    X <- theta_u[iu_id[[j]], ]
    Xy <- t(X) %*% u_er[item_index[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_v)
    mu <- inv_XXV %*% (Xy + inv_Cov_v %*% v_mu[j, ])
    
    #多変量正規分布から特徴ベクトルをサンプリング
    theta_v[j, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  #行列分解のパラメータを更新
  uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)
  
  
  ##交互作用項の特徴ベクトルをサンプリング
  #モデルの応答変数
  z_er <- U - user_mu - uv
  
  for(i in 1:hh){
    #データを抽出
    zi <- z[user_index[[i]], ]
    theta_zi <- t(theta_z[i, , ])
    
      for(j in 1:(k-1)){
      #応答変数を設定
      er <- z_er[user_index[[i]]] - as.numeric(zi[, index_list22[[j]]] %*%
                                                 (theta_zi[index_list21[[j]], ] * theta_zi[index_list23[[j]], ]) %*% vec2)
      
      #交互作用項の事後分布のパラメータ
      X <- z_list[[i]][, , j] %*% theta_zi[index_list11[[j]], ]
      Xy <- t(X) %*% er
      inv_XXV <- solve(t(X) %*% X + inv_Cov_z[, , j])
      mu <- as.numeric(inv_XXV %*% (Xy + inv_Cov_z[, , j] %*% z_mu[i, , j]))   #事後分布の平均
      
      #多変量正規分布から交互作用項をサンプリング
      theta_z[i, , j] <- mvrnorm(1, mu, sigma^2*inv_XXV)
      theta_zi[j, ] <- theta_z[i, , j]
    }
  }
  
  #交互作用項のパラメータを更新
  z_vec <- rep(0, N)
  for(j in 1:m){
    z_vec <- z_vec + z[, allocation_vec21[j]] * (theta_z[user_id, , j_data_vec11[j]] * theta_z[user_id, , allocation_vec11[j]]) %*% vec2
  }
  z_vec <- as.numeric(z_vec)
  
  
  ##ユーザーの回帰ベクトルの階層モデルのパラメータをサンプリング
  #多変量回帰モデルからパラメータをサンプリング
  out <- rmultireg(theta_x, u, Deltabar1, ADelta1, nu1, V1)
  alpha_x <- out$B; x_mu <- u %*% alpha_x   
  Cov_x <- out$Sigma; inv_Cov_x <- solve(Cov_x)
  
  ##ユーザー特徴行列の階層モデルのパラメータをサンプリング
  #多変量回帰モデルからパラメータをサンプリング
  out <- rmultireg(theta_u, u, Deltabar2, ADelta2, nu2, V2)
  alpha_u <- out$B; u_mu <- u %*% alpha_u   
  Cov_u <- out$Sigma; inv_Cov_u <- solve(Cov_u)
  
  ##アイテムの特徴行列の階層モデルのパラメータをサンプリング
  #多変量回帰モデルからパラメータをサンプリング
  out <- rmultireg(theta_v, v, Deltabar3, ADelta3, nu3, V3)
  alpha_v <- out$B; v_mu <- v %*% alpha_v   
  Cov_v <- out$Sigma; inv_Cov_v <- solve(Cov_v)
  
  ##交互作用項の階層モデルのパラメータをサンプリング
  #多変量回帰モデルパラメータをサンプリング
  for(j in 1:(k-1)){
    out <- rmultireg(theta_z[, , j], u, Deltabar4, ADelta4, nu4, V4)
    alpha_z[, , j] <- out$B; z_mu[, , j] <- u %*% alpha_z[, , j]
    Cov_z[, , j] <- out$Sigma; inv_Cov_z[, , j] <- solve(Cov_z[, , j])
  }
  
  ##対数尤度を計算
  mu <- user_mu + uv + z_vec   #潜在効用の期待値
  prob <- pnorm(mu, 0, sigma)   #購買確率
  prob[prob==1] <- 0.9999999; prob[prob==0] <- 0.0000001
  LL <- sum(y*log(prob) + (1-y)*log(1-prob))   #対数尤度
  print(c(LL, LLst))
}


