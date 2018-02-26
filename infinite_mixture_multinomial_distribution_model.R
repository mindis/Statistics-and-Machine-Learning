#####無限混合多項分布モデル#####
library(MASS)
library(mclust)
library(reshape2)
library(bayesm)
detach("package:bayesm", unload=TRUE)
library(ExtDist)
library(extraDistr)
library(matrixStats)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(1853)

####データの発生####
#データの設定
n <- 1000   #セグメントあたりのサンプル数
seg <- 6   #セグメント数
N <- n*seg   #総サンプル数
k <- 12   #変数数
w <- rpois(N, rgamma(N, 60, 0.6))   #サンプルあたりの頻度

#セグメントの設定
seg_id <- rep(1:seg, rep(n, seg))


####応答変数の発生####
##セグメントごとにパラメータの設定
alpha <- rep(0.15, k)
p <- extraDistr::rdirichlet(seg, alpha)

##多項分布よりデータを発生
Data <- matrix(0, nrow=N, ncol=k)
for(i in 1:seg){
  index <- which(seg_id==i)
  Data[index, ] <- rmnom(length(index), w[index], p[i, ])
}
colnames(Data) <- 1:k
storage.mode(Data) <- "integer"

####マルコフ連鎖モンテカルロ法で無限次元混合多項分布モデルを推定####
##アルゴリズムの設定
R <- 10000
keep <- 2
disp <- 20
sbeta <- 1.5
iter <- 0

##事前分布の設定
tau <- 1   #ディクレリ分布の事前分布
alpha <- 1   #CRPの事前分布

##初期値の設定
seg0 <- 2   #初期セグメントは2つ
r <- c(0.5, 0.5)   #混合率の初期値
par_mean <- rep(1/ncol(Data), ncol(Data))   #CRP用のパラメータ

#初期セグメントを設定
z <- matrix(0, nrow=N, ncol=seg0)
out <- kmeans(Data, seg0)   #kmeans法

#セグメント割当
z0 <- out$cluster
for(i in 1:seg0){z[z0==i, i] <- 1}   

#セグメントごとのパラメータ
oldpar0 <- extraDistr::rdirichlet(2, colSums(Data)/sum(Data)+1)
oldpar <- (oldpar0+0.001) / matrix(rowSums(oldpar0+0.001), nrow=2, ncol=k, byrow=T)

#パラメータの格納用配列
max_seg <- 20
Z <- matrix(0, nrow=N, ncol=max_seg)
P <- array(0, dim=c(max_seg, k, R/keep))
storage.mode(Z) <- "integer"

#データの設定
const <- lfactorial(rowSums(Data)) - rowSums(lfactorial(Data))


####MCMCでパラメータをサンプリング####
for(rp in 1:R){
  
  ##多項分布の混合尤度を計算
  #パラメータごとに対数尤度を計算
  LLind0 <- Data %*% t(log(oldpar))
  
  #新しい潜在変数の尤度の計算と尤度の結合
  LL_new <- Data %*% log(par_mean)   #新しい潜在変数の対数尤度
  LLi0 <- cbind(LLind0, LL_new)
  LLi <- exp(LLi0 - rowMaxs(LLi0))   #尤度に変換
  
  ##CRPの計算
  gamma0 <- cbind(matrix(colSums(z), nrow=N, ncol=ncol(z), byrow=T) - z, alpha)
  gamma1 <- LLi * gamma0/(N-1-alpha)
  
  ##多項分布より潜在変数をサンプリング
  z_rate <- gamma1 / rowSums(gamma1)   #潜在変数zの割当確率
  z <- rmnom(N, 1, z_rate)
  z <- z[, colSums(z) > 0]

  ##多項分布のパラメータを更新
  dir_par <- t(t(Data) %*% z) + tau   #ディクレリ分布のパラメータ
  oldpar <- extraDistr::rdirichlet(ncol(z), dir_par)   #多項分布からパラメータを生成


  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    if(rp >= R/2){Z[, 1:ncol(z)] <- Z[, 1:ncol(z)] + z}   #繰り返し数が最大反復数の半分を超えたらパラメータを格納
    P[1:nrow(oldpar), , mkeep] <- oldpar
    
    if(rp%%disp==0){
      print(rp)
      print(colSums(z))
      #print(round(rbind(oldpar, p), 3))
    }
  }
}

####サンプリング結果の可視化と要約####
#バーンイン期間
burnin1 <- R/(keep+2)   
burnin2 <- 1000

##サンプリング結果をプロット
matplot(t(P[1, , 1:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(P[2, , 1:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(P[3, , 1:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(P[4, , 1:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(P[5, , 1:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(P[6, , 1:(R/keep)]), type="l", ylab="パラメータ")

##サンプリング結果の事後平均
mcmc_seg <- sum(colSums(Z) > 10000)   #推定されたセグメント数

#潜在変数zの推定量
round(Z_mu <- (Z/rowSums(Z))[, colSums(Z) > 0], 3)   #潜在変数の割当確率
colnames(Z_mu) <- 1:ncol(Z_mu)
round(colMeans(Z_mu), 3)   #混合率

#多項分布のパラメータの推定量
p_mu <- matrix(0, nrow=mcmc_seg, ncol=k)
for(i in 1:mcmc_seg){
  p_mu[i, ] <- colMeans(t(P[i, , burnin1:(R/keep)]))
}
round(rbind(p_mu, p), 3)   #真のパラメータと比較



