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
  alpha_z <- array(0, dim=c(ncol(u), s2, ncol(z)))
  
  for(j in 1:ncol(u)){
    if(j==1){
      alpha_x[j, ] <- runif(ncol(x), -0.4, 0.3)
      alpha_u[j, ] <- runif(s1, -0.4, 0.2)
      alpha_z[j, , ] <- matrix(rnorm(s2*ncol(z), 0, 0.2), nrow=s2, ncol=ncol(z))
    } else {
      alpha_x[j, ] <- runif(ncol(x), -0.4, 0.3)
      alpha_u[j, ] <- runif(s1, -0.35, 0.2)
      alpha_z[j, , ] <- matrix(rnorm(s2*ncol(z), 0, 0.2), nrow=s2, ncol=ncol(z))
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
  theta_z <- array(0, c(hh, s2, ncol(z)))
  for(j in 1:ncol(z)){
    theta_z[, , j] <- u %*% alpha_z[, , j] + mvrnorm(hh, rep(0, s2), Cov_z)   #交互作用のパラメータ
  }
  theta_xt <- theta_x; theta_ut <- theta_u; theta_vt <- theta_v; theta_zt <- theta_z
  
  ##正規分布から効用と購買ベクトルを生成
  #変量効果のパラメータ
  x_mu <- as.numeric((x * theta_x[user_id, ]) %*% rep(1, ncol(x)))
  
  #行列分解のパラメータ
  uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)
  
  #交互作用のパラメータ
  z_mu <- rep(0, N)
  j_data_vec11 <- as.numeric(j_data11)[as.numeric(j_data11) > 0]
  allocation_vec11 <- as.numeric(allocation_index11)[as.numeric(allocation_index11) > 0]
  allocation_vec21 <- as.numeric(allocation_index21)[as.numeric(allocation_index21) > 0]
  for(j in 1:length(j_data_vec11)){
    z_mu <- z_mu + z[, allocation_vec21[j]] * (theta_z[user_id, , j_data_vec11[j]] * theta_z[user_id, , allocation_vec11[j]]) %*% vec2
  }
  z_mu <- as.numeric(z_mu)
  
  #潜在効用と応答変数を生成
  mu <- x_mu + z_mu + uv   #期待値
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
Cov_x <- Cov_xt; Cov_u <- Cov_ut; Cov_v <- Cov_vt; Cov_z <- Cov_zt
alpha_x <- alpha_xt; alpha_u <- alpha_ut; alpha_v <- alpha_vt; alpha_z <- alpha_zt
x_mu <- u %*% alpha_x; u_mu <- u %*% alpha_u; v_mu <- v %*% alpha_v
z_mu <- array(0, c(hh, s2, ncol(z)))
for(j in 1:ncol(z)){
  z_mu[, , j] <- u %*% alpha_z[, , j]
}

#モデルパラメータの真値
sigma <- sigmat
theta_x <- theta_xt; theta_u <- theta_ut; theta_v <- theta_vt; theta_z <- theta_zt
user_mu <- as.numeric((x * theta_x[user_id, ]) %*% rep(1, ncol(x)))

#行列分解と交互作用項のパラメータ
uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)
z_mu <- rep(0, N)
for(j in 1:length(j_data_vec11)){
  z_mu <- z_mu + z[, allocation_vec21[j]] * (theta_z[user_id, , j_data_vec11[j]] * theta_z[user_id, , allocation_vec11[j]]) %*% vec2
}
z_mu <- as.numeric(z_mu)


##パラメータの初期値


