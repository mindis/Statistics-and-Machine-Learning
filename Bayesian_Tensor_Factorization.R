#####Bayesian Tensor Factorization#####
library(MASS)
library(matrixStats)
library(FAdist)
library(NMF)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

####データの発生####
##データの設定
hh <- 5000   #ユーザー数
item <- 1000   #アイテム数
time <- 12   #観測機関数
N0 <- hh*item*time
k <- 10   #基底数
vec <- rep(1, k)

#IDを設定
user_id0 <- rep(1:hh, item*time)
item_id0 <- rep(rep(1:item, rep(hh, item)), time)
time_id0 <- rep(1:time, rep(hh*item, time))
id <- cbind(user_id0, item_id0, time_id0)

##テンソル分解の定義に基づきデータを生成
for(rp in 1:1000){
  print(rp)
  #多変量正規分布のパラメータを設定
  mu01 <- rep(1.2, k); mu02 <- rep(0.9, k); mu03 <- rep(0.5, k)
  tau01 <- tau02 <- tau03 <- diag(0.2, k)
  sigma <- sigmat <- 0.5

  #多変量正規分布から特徴行列を生成
  W0 <- WT0 <- mvrnorm(hh, mu01, tau01)
  H0 <- HT0 <- mvrnorm(item, mu02, tau02)
  C0 <- CT0 <- mvrnorm(time, mu03, tau03)
  
  #多変量正規分布から評点スコアを生成
  WHC0 <- array(0, dim=c(hh, item, time))
  for(j in 1:k){
    WHC0 <- WHC0 + W0[, j] %o% t(H0)[j, ] %o% t(C0)[j, ]
  }
  #whc0 <- as.numeric((W0[user_id0, ] * t(H0)[item_id0, ] * t(C0)[time_id0, ]) %*% vec)   #こちらでもok
  whc0 <- as.numeric(WHC0)   #テンソルをベクトルに変換
  y_vec0 <- whc0 + rnorm(hh*item*time, 0, sigma)   #誤差を生成
  
  #収束条件
  if(mean(y_vec0) < 5.5 & mean(y_vec0) > 4.5 & sd(y_vec0) > 1.5 & sd(y_vec0) < 2.0) break
}

#応答変数を1～10に変換する
y0 <- round(y_vec0)
y0[y0 > 10] <- 10; y0[y0 < 1] <- 1


##欠損ベクトルを生成
#欠損確率を生成
user_prob <- rbeta(hh, 10, 50)
item_prob <- rbeta(item, 15, 55)
time_prob <- rbeta(time, 60, 140)
prob <- user_prob[user_id0]*item_id0[item_id0]*time_prob[time_id0]

#ベルヌーイ分布から欠損ベクトルを生成
z_vec <- rbinom(N0, 1, prob)
N <- sum(z_vec)
y <- y0[z_vec==1]; y_vec <- y_vec0[z_vec==1]; whc <- whc0[z_vec==1]
hist(y_vec, breaks=25, col="grey", main="潜在的なスコア分布", xlab="スコア")
hist(y, breaks=25, col="grey", main="整数でのスコア分布", xlab="スコア")

#欠損ベクトルからidを再構成
user_id <- user_id0[z_vec==1]
item_id <- item_id0[z_vec==1]
time_id <- time_id0[z_vec==1]


####マルコフ連鎖モンテカルロ法でテンソル分解を推定####
##アルゴリズムの設定
R <- 2000
keep <- 2
disp <- 10
iter <- 0

##事前分布の設定
theta <- rep(0, k)
tau <- 100 * diag(k)
inv_tau <- solve(tau)
s0 <- 1.0
v0 <- 1.0

##真値の設定
W <- WT0
H <- HT0
C <- CT0
sigma <- sigmat

##初期値の設定
sigma <- sd(y)
W <- mvrnorm(hh, rep(0.8, k), 0.2 * diag(k))
H <- mvrnorm(item, rep(0.8, k), 0.2 * diag(k))
C <- mvrnorm(time, rep(0.8, k), 0.2 * diag(k))

##パラメータの格納用配列
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(item, k, R/keep))
C_array <- array(0, dim=c(time, k, R/keep))
SIGMA <- rep(0, R/keep)

##事前分布の設定
##インデックスを設定
user_index <- item_index <- time_index <- list()
ui_id <- ut_id <- iu_id <- it_id <- tu_id <- ti_id <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
  ui_id[[i]] <- item_id[user_index[[i]]]
  ut_id[[i]] <- time_id[user_index[[i]]]
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
  iu_id[[j]] <- user_id[item_index[[j]]]
  it_id[[j]] <- time_id[item_index[[j]]]
}
for(j in 1:time){
  time_index[[j]] <- which(time_id==j)
  tu_id[[j]] <- user_id[time_index[[j]]]
  ti_id[[j]] <- item_id[time_index[[j]]]
}
const1 <- hh / 1.5  #正規化定数
const2 <- item / 1.5

##対数尤度の基準値
LLst <- sum(dnorm(y, mean(y), sd(y), log=TRUE))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##ユーザー特徴行列をサンプリング
  for(i in 1:hh){
    #特徴ベクトルのパラメータ
    X <- H[ui_id[[i]], ] * C[ut_id[[i]], ]
    inv_XXV <- solve(t(X) %*% X +inv_tau)
    beta_mu <- inv_XXV %*% t(X) %*% y[user_index[[i]]]   #多変量正規分布の平均ベクトル
    cov <- sigma^2 * inv_XXV
    
    #多変量正規分布からパラメータをサンプリング
    w <- mvrnorm(1, beta_mu, cov)
    W[i, ] <- w
  }
  W <- W / matrix(colSums(W), nrow=hh, ncol=k, byrow=T) * hh
  
  ##アイテム特徴行列をサンプリング
  for(j in 1:item){
    #特徴ベクトルのパラメータ
    X <- W[iu_id[[j]], ] * C[it_id[[j]], ]
    inv_XXV <- solve(t(X) %*% X +inv_tau)
    beta_mu <- inv_XXV %*% t(X) %*% y[item_index[[j]]]   #多変量正規分布の平均ベクトル
    cov <- sigma^2 * inv_XXV
    
    #多変量正規分布からパラメータをサンプリング
    h <- mvrnorm(1, beta_mu, cov)
    H[j, ] <- h
  }
  H <- H / matrix(colSums(H), nrow=item, ncol=k, byrow=T) * item
  
  ##時間の特徴行列をサンプリング
  for(j in 1:time){
    #特徴ベクトルのパラメータ
    X <- W[tu_id[[j]], ] * H[ti_id[[j]], ]
    inv_XXV <- solve(t(X) %*% X +inv_tau)
    beta_mu <- inv_XXV %*% t(X) %*% y[time_index[[j]]]   #多変量正規分布の平均ベクトル
    cov <- sigma^2 * inv_XXV
    
    #多変量正規分布からパラメータをサンプリング
    c <- mvrnorm(1, beta_mu, cov)
    C[j, ] <- c
  }
  
  
  ##モデルの標準偏差をサンプリング
  #モデルの誤差を推定
  mu <- (W[user_id, ] * H[item_id, ] * C[time_id, ]) %*% vec   #モデルの平均
  er <- y - mu   #モデルの誤差
  
  #逆ガンマ分布のパラメータ
  s <- as.numeric(t(er) %*% er) + s0
  v <- N + v0
  
  #逆ガンマ分布から標準偏差をサンプリング
  sigma <- sqrt(1/rgamma(1, v/2, s/2))
  
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    W_array[, , mkeep] <- W
    H_array[, , mkeep] <- H
    C_array[, , mkeep] <- C
    SIGMA[mkeep] <- sigma
  }
  
  #対数尤度を計算
  if(rp%%disp==0){
    LLi <- dnorm(y, mu, sigma, log=TRUE)
    LL <- sum(LLi)
    
    #サンプリング結果の表示
    print(rp)
    print(c(LL, LLst))
    print(round(c(sigma, sigmat), 3))
  }
}

matplot(t(W_array[1, , ]), type="l")
matplot(t(H_array[1, , ]), type="l")
matplot(t(C_array[1, , ]), type="l")
C_array
W
H
C
C_array
round(W_array, 3)


sum((y - rowSums(W[user_id, ] * H[item_id, ] * C[time_id, ]))^2)


