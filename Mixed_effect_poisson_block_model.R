#####mixed effect poisson block model#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(flexmix)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

set.seed(506832)

####データの発生####
##データの設定
d <- 150   #アイテム数
k <- 7   #潜在変数数
N <- d*(d-1)/2   #総サンプル数
vec <- rep(1, k)

##IDを設定
id1 <- id2 <- c()
for(i in 1:(d-1)){
  id1 <- c(id1, rep(i, length((i+1):d)))
  id2 <- c(id2, (i+1):d)
}

##潜在変数の生成
#ディリクレ分布からパラメータを生成
alpha0 <- rep(0.1, k)
theta <- thetat <- extraDistr::rdirichlet(d, alpha0)
Z1 <- rmnom(N, 1, theta[id1, ])
Z2 <- rmnom(N, 1, theta[id2, ])
z1_vec <- as.numeric(Z1 %*% 1:k); z2_vec <- as.numeric(Z2 %*% 1:k)


##応答変数の生成
#パラメータを生成
cov <- covt <- 0.75
beta <- betat <- 0.8
alphat <- alpha <- rnorm(d, beta, cov)   #変量効果のパラメータ
phi0 <- matrix(rnorm(k*k, 0, 0.85), nrow=k, ncol=k)   #潜在変数のパラメータ
phi0[upper.tri(phi0)] <- 0
phi <- phi0 + t(phi0)
diag(phi) <- diag(phi0)
phit <- phi

#ポアソン分布の平均構造
mu <- alpha[id1] + alpha[id2] + (phi[z1_vec, ] * Z2) %*% vec
lambda <- exp(mu)


#ポアソン分布から応答変数を生成
y <- rpois(N, lambda)
sum(y); mean(y)
hist(y, xlab="頻度", main="アイテム間の出現頻度", col="grey", breaks=25)


####マルコフ連鎖モンテカルロ法でmixed effect poisson block modelを推定####
##変量効果ポアソン回帰モデルの対数尤度
loglike1 <- function(alpha, theta, y, y_factorial, z1, Z2, vec, id1, id2){
  
  #尤度を定義する
  lambda <- exp(alpha[id1] + alpha[id2] + (phi[z1, ] * Z2) %*% vec)   #平均構造
  LLi <- as.numeric(y*log(lambda)-lambda - y_factorial)   #対数尤度
  LL <- sum(LLi)   #対数尤度の和
  
  #結果を返す
  LL_value <- list(LLi=LLi, LL=LL)
  return(LL_value)
}

loglike2 <- function(alpha, phi, y, y_factorial, index, id1, id2){
  #尤度を定義する
  lambda <- exp(alpha[id1[index]] + alpha[id2[index]] + phi)
  LL <- sum(y[index]*log(lambda)-lambda - y_factorial[index])
  return(LL)
}


##アルゴリズムの設定
R <- 10000
keep <- 4  
iter <- 0
burnin <- 1000/keep
disp <- 100

##事前分布の設定
#回帰係数の事前分布
alpha01 <- 0
beta01 <- 0
tau01 <- 0.01

#変量効果の事前分布
alpha02 <- 0
s02 <- 0.01
v02 <- 0.01

#ディリクレ分布の事前分布
alpha03 <- 0.1


##パラメータの真値
theta <- thetat
phi <- phit
alpha <- alphat
beta <- betat
cov <- covt
Zi1 <- Z1; Zi2 <- Z2
z_vec1 <- as.numeric(Zi1 %*% 1:k)
z_vec2 <- as.numeric(Zi2 %*% 1:k)


##初期値の設定
#変量効果の初期値
cov <- 0.5   #階層モデルの標準偏差の初期値
beta <- 0.8   #階層モデルの平均の初期値
mu <- rep(0, d)
for(i in 1:d){
  mu[i] <- mean(y[which(id1==i | id2==i)])
}
alpha <- rnorm(d, 0, 0.5)
rank_mu <- ceiling(rank(mu))
alpha <- sort(x, decreasing=TRUE)[rank_mu]   #変量効果の初期値

#潜在変数のパラメータ
phi0 <- matrix(rnorm(k*k, 0, 0.5), nrow=k, ncol=k)   #潜在変数のパラメータ
phi0[upper.tri(phi0)] <- 0
phi <- phi0 + t(phi0)

#潜在変数の初期値
theta <- extraDistr::rdirichlet(d, rep(0.2, k))
Zi1 <- rmnom(N, 1, theta[id1, ])
Zi2 <- rmnom(N, 1, theta[id2, ])
z1_vec <- as.numeric(Zi1 %*% 1:k); z2_vec <- as.numeric(Zi2 %*% 1:k)


##パラメータの格納用配列
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, k, R/keep))
ALPHA <- matrix(0, nrow=R/keep, ncol=d)
BETA <- rep(0, R/keep)
COV <- rep(0, R/keep)
SEG2 <- SEG1 <- matrix(0, nrow=N, ncol=k)

##定数を計算
y_factorial <- lfactorial(y)
Y_factorial <- matrix(y_factorial, nrow=N, ncol=k*k)
upper_tri <- matrix(as.logical(upper.tri(phi) + diag(1, k)), nrow=k, ncol=k)
vec <- rep(1, k)

##インデックスを設定
#入力データのインデックス
item_list <- list()
seg_list1 <- seg_list2 <- list()
index_list1 <- index_list2 <- list()
for(i in 1:d){
  item_list[[i]] <- which(id1==i | id2==i)
  seg_list1[[i]] <- as.numeric(id1[item_list[[i]]]!=i)
  seg_list2[[i]] <- as.numeric(id2[item_list[[i]]]!=i)
  index_list1[[i]] <- which(item_list[[i]] * seg_list1[[i]] != 0)
  index_list2[[i]] <- which(item_list[[i]] * seg_list2[[i]] != 0)
}
item_vec <- rep(1, d-1)

#パラメータのインデックス
par <- c()
for(i in 1:k){
  for(j in i:k){
    par <- c(par, i, j)
  }
}
index_par <- matrix(par, nrow=k*(k-1)/2+k, ncol=2, byrow=T)


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
    
  ##メトロポリスヘイスティング法で変量効果をサンプリング
  #新しいパラメータをサンプリング
  alphad <- alpha
  alphan <- alphad + rnorm(d, 0, 0.1)
  
  #事前分布の誤差
  er_new <- alphan - beta 
  er_old <- alphad - beta
  
  #対数尤度と対数事前分布を設定
  lognew0 <- loglike1(alphan, phi, y, y_factorial, z1_vec, Zi2, vec, id1, id2)$LLi
  logold0 <- loglike1(alphad, phi, y, y_factorial, z1_vec, Zi2, vec, id1, id2)$LLi
  logpnew1 <- -0.5 * (er_new^2 / cov)
  logpold1 <- -0.5 * (er_old^2 / cov)
  
  #アイテムごとに対数尤度の和を取る
  lognew1 <- logold1 <- rep(0, d)
  for(i in 1:d){
    lognew1[i] <- sum(lognew0[item_list[[i]]])
    logold1[i] <- sum(logold0[item_list[[i]]])
  }
  
  #MHサンプリング
  rand <- runif(d)   #一様分布から乱数を発生
  LLind_diff <- exp(lognew1 + logpnew1 - logold1 - logpold1)   #採択率を計算
  LLind_diff <- ifelse(LLind_diff==Inf, 1, ifelse(LLind_diff==-Inf, 0, LLind_diff))
  alpha2 <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- (alpha2 >= rand)*1 + (alpha2 < rand)*0
  alpha <- flag*alphan + (1-flag)*alphad   #alphaがrandを上回っていたら採択
  
  
  ##階層モデルのパラメータをサンプリング
  #正規分布から平均パラメータをサンプリング
  #beta_mu <- d/(d + tau01) * mean(alpha)
  #beta <- rnorm(1, beta_mu, cov/(d + tau01))
  
  #逆ガンマ分布から標準偏差をサンプリング
  s <- s02 + sum((alpha - mean(alpha))^2)
  v <- v02 + d
  cov <- sqrt(1/rgamma(1, v/2, s/2))
  
  
  ##ギブスサンプリングで潜在変数をサンプリング
  #多項分布から潜在変数をサンプリング
  for(i in 1:d){
    
    #インデックスを抽出
    index <- item_list[[i]]
    z1_allocation <- seg_list1[[i]]
    z2_allocation <- seg_list2[[i]]
    
    #ポアソン分布の平均構造
    Zi_pairs <- Zi1[index, ] * z1_allocation + Zi2[index, ] * z2_allocation   #対となる潜在変数
    Zi_mu <- phi[as.numeric(Zi_pairs %*% 1:k), ]
    lambda <- exp(matrix(alpha[id1[index]] + alpha[id2[index]], nrow=length(index), ncol=k) + Zi_mu)
    
    #対数尤度と潜在変数の割当確率
    LLi <- y[index]*log(lambda)-lambda - y_factorial[index]   #対数尤度
    z_par <- exp(LLi - rowMaxs(LLi)) * matrix(theta[i, ], nrow=length(index), ncol=k, byrow=T)
    z_rate <- z_par / rowSums(z_par)
    
    #多項分布より潜在変数をサンプリング
    index1 <- index_list1[[i]]; index2 <- index_list2[[i]]
    if(length(index1) > 0){
      Zi1[index1, ] <- rmnom(length(index1), 1, z_rate[index1, ])
    } 
    if(length(index2) > 0){
      Zi2[index2, ] <- rmnom(length(index2), 1, z_rate[index2, ])
    }
  }
  #潜在変数行列を変換
  z1_vec <- as.numeric(Zi1 %*% 1:k)
  z2_vec <- as.numeric(Zi2 %*% 1:k)
  Zi1_T <- t(Zi1); Zi2_T <- t(Zi2)
  
  
  ##MH法で潜在変数のパラメータをサンプリング
  for(i in 1:nrow(index_par)){
    
    #インデックスを設定
    index_phi <- index_par[i, ]
    index_z <- which(Zi1[, index_phi[1]]*Zi2[, index_phi[2]]==1)
    
    #新しいパラメータをサンプリング
    phid <- phi[index_phi[1], index_phi[2]]
    phin <- phid + rnorm(1, 0, 0.1)
    
    #対数尤度と対数事前分布を設定
    lognew2 <- loglike2(alpha, phin, y, y_factorial, index_z, id1, id2)
    logold2 <- loglike2(alpha, phid, y, y_factorial, index_z, id1, id2)
    logpnew2 <- -0.5 * (phin^2 / (1/tau01))
    logpold2 <- -0.5 * (phid^2 / (1/tau01))
    
    #MHサンプリング
    rand <- runif(1)   #一様分布から乱数を発生
    LLind_diff <- exp(lognew2 + logpnew2 - logold2 - logpold2)   #採択率を計算
    if(LLind_diff==Inf) LLind_diff <- 1; if(LLind_diff==-Inf) LLind_diff <- 0   #Infの場合の処理
    alpha2 <- (LLind_diff >= 1)*1 + (LLind_diff < 1)*LLind_diff
    
    #alphaの値に基づき新しいbetaを採択するかどうかを決定
    flag <- (alpha2 >= rand)*1 + (alpha2 < rand)*0
    phi[index_phi[1], index_phi[2]] <- flag*phin + (1-flag)*phid   #alphaがrandを上回っていたら採択
  }
  phi[lower.tri(phi)] <- phi[upper.tri(phi)]   #パラメータを対象行列に変更
  phi <- phi - mean(phi)   #パラメータを中央化
  
  
  ##ディリクレ分布からトピック分布をサンプリング
  wsum0 <- matrix(0, nrow=d, ncol=k)
  for(i in 1:d){
    index1 <- item_list[[i]] * seg_list1[[i]]
    index2 <- item_list[[i]] * seg_list2[[i]]
    wsum0[i, ] <- cbind(Zi1_T[, index2], Zi2_T[, index1]) %*% item_vec
  }
  wsum <- wsum0 + alpha03   #ディリクレ分布のパラメータ
  theta <- extraDistr::rdirichlet(d, wsum)   #ディリクレ分布からトピックを生成
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    ALPHA[mkeep, ] <- alpha
    BETA[mkeep] <- beta
    COV[mkeep] <- cov
  }
  
  #トピック割当はバーンイン期間を超えたら格納する
  if(rp%%keep==0 & rp >= burnin){
    SEG1 <- SEG1 + Zi1
    SEG2 <- SEG2 + Zi2
  }
  
  if(rp%%disp==0){
    #サンプリング結果の表示
    print(rp)
    print(sum(lognew1))
    print(round(cbind(phi, phit), 3))
  }
}
matplot(ALPHA[, 1:10], type="l")
matplot(t(PHI[1, , ]), type="l")
plot(1:(R/keep), COV, type="l")
