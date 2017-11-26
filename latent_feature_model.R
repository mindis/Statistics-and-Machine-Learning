#####ベイジアン潜在因子モデル#####
library(MASS)
library(lda)
library(RMeCab)
library(bayesm)
library(extraDistr)
library(matrixStats)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(654978)

####データの発生####
#データを設定
hh <- 2000   #サンプル数
colums <- 100   #変数数
k <- 10   #潜在変数数

#バイナリ行列を発生
Z0 <- matrix(0, nrow=hh, ncol=k)
for(j in 1:k){
  p <- runif(1, 0.25, 0.5)
  Z0[, j] <- rbinom(hh, 1, p)
}
r0 <- colMeans(Z0)   #混合率

#因子行列を発生
X0 <- matrix(rnorm(hh*k, 0, 1), nrow=k, ncol=colums)

#応答変数を発生
Cov0 <- 0.5
y <- Z0 %*% X0 + matrix(rnorm(hh*colums, 0, Cov0), nrow=hh, ncol=colums)


####マルコフ連鎖モンテカルロ法で潜在特徴モデルを推定####
##アルゴリズムの設定
R <- 10000
keep <- 4
iter <- 0

##事前分布の設定
#正規分布の分散の事前分布のパラメータ
Deltabar <- matrix(0, nrow=k, ncol=colums)   #階層モデルの回帰係数の事前分布の分散
ADelta <- 0.01 * diag(rep(1, k))   #階層モデルの回帰係数の事前分布の分散
nu <- colums   #逆ウィシャート分布の自由度
V <- nu * diag(rep(1, colums)) #逆ウィシャート分布のパラメータ
s0 <- 0.01
v0 <- 0.01

##初期値の設定
#初期値の候補を生成
iter <- 5000
Z_array <- array(0, dim=c(hh, k, iter))
storage.mode <- "integer"
sq_error <- rep(0, iter)

for(i in 1:iter){
  print(i)
  for(j in 1:k){
    p <- runif(1, 0.25, 0.6)
    Z_array[, j, i] <- rbinom(hh, 1, p)
  }
  X <- solve(t(Z_array[, , i]) %*% Z_array[, , i]) %*% t(Z_array[, , i]) %*% y
  sq_error[i] <- sum((y - Z_array[, , i] %*% X)^2)
}
min(sq_error)


#ベストな初期値を選択
bt <- which.min(sq_error)
Z <- Z_array[, , bt]
X <- solve(t(Z) %*% Z) %*% t(Z) %*% y
Cov <- 0.1
r <- colMeans(Z)
rm(Z_array)

##サンプリング結果の格納用配列
Zi <- array(0, dim=c(hh, k, R/keep))
FACTOR <- array(0, dim=c(k, colums, R/keep))
SIGMA <- rep(0, R/keep)
storage.mode(Zi) <- "integer"
gc(); gc()


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  z1_rate <- matrix(0, nrow=hh, ncol=k)
  
  ##潜在変数Zをサンプリング
  for(j in 1:k){
    
    #パターン別に対数尤度を計算
    z0 <- z1 <- Z
    z0[, j] <- 0; z1[, j] <- 1
    LLi0 <- rowSums(dnorm(y, z0 %*% X, Cov, log=TRUE))
    LLi1 <- rowSums(dnorm(y, z1 %*% X, Cov, log=TRUE))
    
    #logsumexpの尤度を計算
    LLi <- cbind(LLi0, LLi1)
    
    LLi_max <- matrix(apply(LLi, 1, max), nrow=hh, 2)
    r_matrix <- matrix(c(1-r[j], r[j]), nrow=hh, ncol=2, byrow=T)  

    #割当確率のパラメータを設定
    expl <- r_matrix * exp(LLi - LLi_max)
    expl_log <- log(expl)
    expl_max <- log(max(expl))
    z1_rate[, j] <- exp(expl_log[, 2] - (log(rowSums(exp(expl_log - expl_max))) + expl_max))   #セグメント割当確率
    
    #ベルヌーイ分布より潜在変数を生成
    Z[, j] <- as.integer(z1_rate[, j] > runif(hh))
  }
 
  #混合率を更新
  r <- colMeans(Z)
  
  ##多変量回帰モデルで因子行列Xを更新
  out <- rmultireg(Y=y, X=Z, Bbar=Deltabar, A=ADelta, nu=nu, V=V)
  X <- out$B
  
  ##逆ガンマ分布から標準偏差を更新
  er <- as.numeric(y - Z %*% X)
  s <- s0 + t(er) %*% er
  v <- v0 + hh*colums
  Cov <- sqrt(1/(rgamma(1, v/2, s/2)))   #逆ガンマ分布からsigma^2をサンプリング

  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    Zi[, , mkeep] <- Z
    FACTOR[, , mkeep] <- X
    SIGMA[mkeep] <- Cov
      
    #サンプリング結果の表示
    print(rp)
    print(round(rbind(r, r0), 3))
    print(round(c(Cov, Cov0), 3))
    print(cbind(Z[1:20, ], Z0[1:20, ]))
  }
}

####サンプリング結果の確認と適合度の確認####
burnin <- 250   #バーンイン期間(1000サンプルまで)
RS <- R/keep 

##サンプリング結果の可視化
matplot(t(FACTOR[, 1, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(FACTOR[, 2, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(FACTOR[, 3, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
plot(1:RS, SIGMA, type="l", xlab="サンプリング回数")

##事後要約値
cbind(apply(Zi[, , burnin:RS], c(1, 2), mean), Z0)
round(cbind(t(apply(FACTOR[, , burnin:RS], c(1, 2), mean)), t(X0)), 2)
round(c(mean(SIGMA[burnin:RS]), Cov0), 3)


pnorm(100)
pnorm(-100)
