#####Mean Shift Model#####
options(warn=2)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(SMC)
library(RcppSMC)
library(SMC)
library(KFAS)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(93441)

####データの発生####
##データの設定
d <- 1000   #文書数
v <- 300   #語彙数
a <- rpois(d, rgamma(d, 15, 0.7))   #区切り文章数
f1 <- sum(a)   #総文章数
w <- rtpois(f1, 8, 5, Inf)   #区切り文章ごとの単語数
f2 <- sum(w)   #総語彙数

##IDの設定
w_id1 <- rep(1:d, a)


##パラメータの設定
#事前分布のパラメータ
alpha01 <- rep(0.2, v)   #ディリクレ分布のパラメータ
alpha11 <- 20   #ベータ分布のパラメータ
beta11 <- 100

#パラメータを生成
gamma <- rbeta(d, alpha11, beta11)   #切換え確率のパラメータ

#モデルに基づきデータを生成
WX <- matrix(0, nrow=f1, ncol=v)
storage.mode(WX) <- "integer"
Z1_list <- list()
Z2_list <- list()
theta1_list <- list()
theta2_list <- list()

for(i in 1:d){
  
  #切換え変数を生成
  z <- rbinom(a[i], 1, gamma[i])
  z[1] <- 1   #1文目はz=1
  z_no <- cumsum(z)   #切換え履歴に変換
  
  #ディリクレ分布からパラメータを生成
  theta <- extraDistr::rdirichlet(sum(z), alpha01)
  theta_all <- theta[z_no, ]
  
  #多項分布から単語を生成
  words <- rmnom(a[i], w[w_id1==i], theta_all)
  
  #データを格納
  WX[w_id1==i, ] <- words
  Z1_list[[i]] <- z
  Z2_list[[i]] <- z_no
  theta1_list[[i]] <- theta
  theta2_list[[i]] <- theta_all
}

#リストを変換
Z <- unlist(Z1_list)
Zn <- unlist(Z2_list)
theta <- do.call(rbind, theta1_list)
theta_all <- do.call(rbind, theta2_list)
rm(theta1_list); rm(theta2_list)


####粒子フィルタでMultinomial Mean Shift modelを推定####


##粒子フィルタの設定
s <- 3000   #粒子数
alpha1 <- 0.1
beta1 <- rep(1/v, v)

#パラメータを生成
par <- WX[1, ]
beta0 <- (par + alpha1) / sum(par + alpha1)

#尤度を更新
Li0 <- as.numeric(WX[2, ] %*% t(log(beta0)))
Li1 <- as.numeric(WX[2, ] %*% t(log(beta1)))
Li <- cbind(Li0, Li1)

#パラメータをリサンプリング
resample_par <- exp(Li0 + Li1 - max(Li0 + Li1))
resample_par / sum(resample_par)



#潜在変数zの割当確率を更新
z_par <- exp(Li - rowMaxs(Li))
z_rate <- z_par / rowSums(z_par)
z_rate
rbinom(s, z_rate[, 2])


WX[2, ] %*% log((WX[1, ] + alpha1) / sum(WX[1, ] + alpha1))
WX[2, ] %*% log(rep(alpha1, v) / (v * alpha1))


