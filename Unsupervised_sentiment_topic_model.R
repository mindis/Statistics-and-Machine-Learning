#####Unsupervised Sentiment topic model#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(364023)

####データの発生####
##データの設定
k1 <- 15   #トピック数
k2 <- 3   #センチメント数
d <- 3000   #文書数
r <- 5   #評価数
v <- 1000   #語彙数
index_v1 <- 1:700   #一般語の語彙割当
index_v2 <- 701:v   #極性語の語彙割当
v1 <- length(index_v1)
v2 <- length(index_v2)
spl <- split(index_v2, 1:3)
w <- rpois(d, rgamma(d, 75, 0.55))   #単語数
f <- sum(w)   #総単語数

##IDの設定
t_id <- c()
d_id <- rep(1:d, w)
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}

##パラメータの生成
#ディリクレ分布とベータ分布の事前分布
alpha01 <- rep(0.1, k1) 
alpha02 <- c(0.1, 0.25, 0.2, 0.25, 0.1)
alpha11 <- c(rep(0.1, v1), rep(0.001, v2))
beta01 <- 100.0; beta02 <- 40.0

#極性分布の事前分布
alpha12 <- matrix(0, nrow=r, ncol=v)
alpha12[1, ] <- c(rep(0.001, v1), rep(0.75, length(spl[[1]])), rep(0.05, length(spl[[2]])), rep(0.0025, length(spl[[3]])))
alpha12[2, ] <- c(rep(0.001, v1), rep(0.6, length(spl[[1]])), rep(0.15, length(spl[[2]])), rep(0.025, length(spl[[3]])))
alpha12[3, ] <- c(rep(0.001, v1), rep(0.2, length(spl[[1]])), rep(1.0, length(spl[[2]])), rep(0.2, length(spl[[3]])))
alpha12[4, ] <- c(rep(0.001, v1), rep(0.025, length(spl[[1]])), rep(0.15, length(spl[[2]])), rep(0.6, length(spl[[3]])))
alpha12[5, ] <- c(rep(0.001, v1), rep(0.0025, length(spl[[1]])), rep(0.05, length(spl[[2]])), rep(0.75, length(spl[[3]])))


for(rp in 1:1000){
  print(rp)
  
  #事前分布からパラメータを生成
  theta <- extraDistr::rdirichlet(d, alpha01)
  ohm <- extraDistr::rdirichlet(d, alpha02)
  phi <- extraDistr::rdirichlet(k1, alpha11)
  omega <- extraDistr::rdirichlet(r, alpha12)
  lambda <- rbeta(d, beta01, beta02)
  
  
  ##モデルに基づきデータを生成
  word_list <- list()
  Z1_list <- Z2_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  y <- matrix(0, nrow=d, ncol=r)
  y_vec <- rep(0, d)
  
  for(i in 1:d){
    
    #一般語か極性語かどうかを生成
    z1 <- rbinom(w[i], 1, lambda[i])
    index <- which(z1==1)
    
    #評価値を生成
    y[i, ] <- rmnom(1, 1, ohm)
    y_vec[i] <- as.numeric(y[i, ] %*% 1:r)
    
    #一般語のトピックを生成
    z2 <- matrix(0, nrow=w[i], ncol=k1)
    z2[index, ] <- rmnom(sum(z1), 1, theta[i, ])
    z2_vec <- as.numeric(z2 %*% 1:k1)
    
    #単語を生成
    word <- matrix(0, nrow=w[i], ncol=v)
    word[index, ] <- rmnom(sum(z1), 1, phi[z2_vec[index], ])
    word[-index, ] <- rmnom(sum(1-z1), 1, omega[y_vec[i], ])
    
    #データを格納
    WX[i, ] <- colSums(word)
    wd_list[[i]] <- as.numeric(word %*% 1:v)
    Z1_list[[i]] <- z1
    Z2_list[[i]] <- z2
  }
  if(min(colSums(WX)) > 0) break
}

##リストを変換
wd <- unlist(wd_list)
Z1 <- unlist(Z1_list)
Z2 <- do.call(rbind, Z2_list)


####マルコフ連鎖モンテカルロ法でUnsupervised Sentiment topic modelを推定####
##トピックモデルの尤度と負担率を定義する関数
burden_fr <- function(theta, phi, wd, w, k){
  #負担係数を計算
  Bur <- theta[w, ] * t(phi)[wd, ]   #尤度
  Br <- Bur / rowSums(Bur)   #負担率
  r <- colSums(Br) / sum(Br)   #混合率
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##アルゴリズムの設定
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##インデックスの設定
doc_list <- list()
doc_vec <- list()
wd_list <- list()
wd_vec <- list()

for(i in 1:d){
  doc_list[[i]] <- which(d_id==i)
  doc_vec[[i]] <- rep(1, length(doc_list[[i]]))
}
for(j in 1:v){
  wd_list[[j]] <- which(wd==j)
  wd_vec[[j]] <- rep(1, length(wd_list[[j]]))
}

##事前分布の設定
#トピックモデルの事前分布
alpha01 <- 0.1
alpha02 <- 0.01

#ベータ分布の事前分布
beta01 <- 1
beta02 <- 1

##パラメータの真値
beta1 <- betat1; beta2 <- betat2
theta1 <- thetat1; theta2 <- thetat2
phi1 <- phit1; phi2 <- phit2
omega1 <- omegat1; omega2 <- omegat2

##初期値の設定
beta1 <- beta2 <- rep(0.5, d)
theta1 <- extraDistr::rdirichlet(d, rep(1.0, k1))
theta2 <- extraDistr::rdirichlet(d, rep(1.0, k2))
phi1 <- extraDistr::rdirichlet(k1, c(rep(1.5, v11), rep(0.1, v12)))
phi2 <- extraDistr::rdirichlet(k2, c(rep(0.1, v11), rep(1.5, v12)))
omega1 <- extraDistr::rdirichlet(k2, c(rep(1.5, v21), rep(0.1, v22)))
omega2 <- as.numeric(extraDistr::rdirichlet(1, c(rep(0.1, v21), rep(1.5, v22))))


##パラメータの格納用配列
BETA1 <- matrix(0, nrow=R/keep, ncol=d)
BETA2 <- matrix(0, nrow=R/keep, ncol=d)
THETA1 <- array(0, dim=c(d, k1, R/keep))
THETA2 <- array(0, dim=c(d, k2, R/keep))
PHI1 <- array(0, dim=c(k1, v1, R/keep))
PHI2 <- array(0, dim=c(k2, v1, R/keep))
OMEGA1 <- array(0, dim=c(k2, v2, R/keep))
OMEGA2 <- matrix(0, nrow=R/keep, ncol=v2)
FLAG1 <- rep(0, f1)
FLAG2 <- rep(0, f2)
SEG11 <- matrix(0, nrow=f1, ncol=k1)
SEG12 <- matrix(0, nrow=f1, ncol=k2)
SEG21 <- matrix(0, nrow=f2, ncol=k2)

##基準の対数尤度を計算
#ユニグラムモデルの対数尤度
LLst <- sum(WX1 %*% log(colSums(WX1)/f1)) + sum(WX2 %*% log(colSums(WX2)/f2))

#真値での対数尤度
phit11 <- (phit1+10^-100) / rowSums(phit1+10^-100); phit12 <- (phit2+10^-100) / rowSums(phit2+10^-100)
omegat11 <- (omegat1+10^-100) / rowSums(omegat1+10^-100); omegat12 <- (omegat2+10^-100) / sum(omegat2+10^-100)
LLbest1 <- sum(log((1-flag1)*rowSums(thetat1[d_id1, ]*t(phit11)[wd1, ]) + flag1*rowSums(thetat2[d_id1, ]*t(phit12)[wd1, ])))
LLbest2 <- sum(log(flag2*rowSums(thetat2[d_id2, ]*t(omegat11)[wd2, ]) + (1-flag2)*omegat12[wd2]))      
LLbest <- LLbest1 + LLbest2


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##文書1のトピックの生成過程をサンプリング
  #期待尤度を定義
  phi1_T <- t(phi1); phi2_T <- t(phi2)
  Lho11 <- theta1[d_id1, ] * phi1_T[wd1, ]   
  Lho12 <- theta2[d_id1, ] * phi2_T[wd1, ] 
  topic_par11 <- Lho11 %*% topic_vec1   #独自成分の期待尤度
  topic_par12 <- Lho12 %*% topic_vec2   #共通部分の期待尤度
  
  #ベルヌーイ分布から生成過程をサンプリング
  beta_par1 <- (1-beta1)[d_id1] * topic_par11
  beta_par2 <- beta1[d_id1] * topic_par12
  flag_rate1<- beta_par2 / (beta_par1+beta_par2)   #潜在変数の割当確率
  flag_vec1 <- rbinom(f1, 1, flag_rate1)
  
  
  ##文書2のトピックの生成過程をサンプリング
  #期待尤度を定義
  omega1_T <- t(omega1)
  Lho21 <- theta2[d_id2, ] * omega1_T[wd2, ]
  topic_par21 <- Lho21 %*% topic_vec2   #共通部分の期待尤度
  topic_par22 <- omega2[wd2]   #独自成分の期待尤度
  
  #ベルヌーイ分布から生成過程をサンプリング
  beta_par1 <- beta2[d_id2] * topic_par21
  beta_par2 <- (1-beta2)[d_id2] * topic_par22
  flag_rate2 <- beta_par1 / (beta_par1+beta_par2)   #潜在変数の割当確率
  flag_vec2 <- rbinom(f2, 1, flag_rate2)
  
  
  ##混合率をサンプリング
  #ベータ分布のパラメータ
  freq1 <- tapply(flag_vec1, d_id1, sum)
  freq2 <- tapply(flag_vec2, d_id2, sum)
  
  #ベータ分布からパラメータをサンプリング
  beta1 <- rbeta(d, freq1+beta01, w1-freq1+beta02)
  beta2 <- rbeta(d, freq2+beta01, w2-freq2+beta02)
  
  
  ##文書1のトピックをサンプリング
  #データの設定
  Zi11 <- matrix(0, nrow=f1, ncol=k1)
  Zi12 <- matrix(0, nrow=f1, ncol=k2)
  index1 <- which(flag_vec1==0)
  
  #多項分布から独自トピックをサンプリング
  z_rate11 <- (Lho11 / as.numeric(topic_par11))[index1, ]   #トピックの割当確率
  Zi11[index1, ] <- rmnom(length(index1), 1, z_rate11)   #トピックのサンプリング
  zi11_vec <- as.numeric(Zi11 %*% 1:k1)
  Zi11_T <- t(Zi11)
  
  #多項分布から共通トピックをサンプリング
  z_rate12 <- (Lho12 / as.numeric(topic_par12))[-index1, ]   #トピックの割当確率
  Zi12[-index1, ] <- rmnom(f1-length(index1), 1, z_rate12)   #トピックのサンプリング
  zi12_vec <- as.numeric(Zi12 %*% 1:k2)
  Zi12_T <- t(Zi12)
  
  
  ##文書2のトピックをサンプリング
  #データの設定
  Zi21 <- matrix(0, nrow=f2, ncol=k2)
  index2 <- which(flag_vec2==1)
  
  #多項分布から独自トピックをサンプリング
  z_rate21 <- (Lho21 / as.numeric(topic_par21))[index2, ]   #トピックの割当確率
  Zi21[index2, ] <- rmnom(length(index2), 1, z_rate21)   #トピックのサンプリング
  zi21_vec <- as.numeric(Zi21 %*% 1:k2)
  Zi21_T <- t(Zi21)
  
  
  ##トピック分布のパラメータをサンプリング
  #文書1の独自トピック分布をサンプリング
  wsum0 <- matrix(0, nrow=d, ncol=k1)
  for(i in 1:d){
    wsum0[i, ] <- Zi11_T[, doc1_list[[i]]] %*% doc1_vec[[i]]
  }
  wsum <- wsum0 + alpha01   #ディリクレ分布のパラメータ
  theta1 <- extraDistr::rdirichlet(d, wsum)   #ディリクレ分布からパラメータをサンプリング
  
  #共通トピック分布をサンプリング
  wsum0 <- matrix(0, nrow=d, ncol=k2)
  for(i in 1:d){
    wsum0[i, ] <- Zi12_T[, doc1_list[[i]]] %*% doc1_vec[[i]] + Zi21_T[, doc2_list[[i]]] %*% doc2_vec[[i]]
  }
  wsum <- wsum0 + alpha01   #ディリクレ分布のパラメータ
  theta2 <- extraDistr::rdirichlet(d, wsum)   #ディリクレ分布からパラメータをサンプリング
  
  
  ##単語分布のパラメータをサンプリング
  #文書1の単語分布をサンプリング
  vsum01 <- matrix(0, nrow=k1, ncol=v1)
  vsum02 <- matrix(0, nrow=k2, ncol=v1)
  for(j in 1:v1){
    vsum01[, j] <- Zi11_T[, wd1_list[[j]], drop=FALSE] %*% wd1_vec[[j]]
    vsum02[, j] <- Zi12_T[, wd1_list[[j]], drop=FALSE] %*% wd1_vec[[j]]
  }
  vsum1 <- vsum01 + alpha02; vsum2 <- vsum02 + alpha02   #ディリクレ分布のパラメータ
  phi1 <- extraDistr::rdirichlet(k1, vsum1)   #ディリクレ分布から独自トピックの単語分布をサンプリング
  phi2 <- extraDistr::rdirichlet(k2, vsum2)   #ディリクレ分布から共通トピックの単語分布をサンプリング
  
  #文書2の単語分布をサンプリング
  vsum01 <- matrix(0, nrow=k2, ncol=v2)
  vsum02 <- rep(0, v2)
  for(j in 1:v2){
    vsum01[, j] <- Zi21_T[, wd2_list[[j]], drop=FALSE] %*% wd2_vec[[j]]
    vsum02[j] <- (1-flag2[wd2_list[[j]], drop=FALSE]) %*% wd2_vec[[j]]
  }
  vsum1 <- vsum01 + alpha02; vsum2 <- vsum02 + alpha02   #ディリクレ分布のパラメータ
  omega1 <- extraDistr::rdirichlet(k2, vsum1)   #ディリクレ分布から共通トピックの単語分布をサンプリング
  omega2 <- as.numeric(extraDistr::rdirichlet(1, vsum2))   #ディリクレ分布から一般語の単語分布をサンプリング
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    BETA1[mkeep, ] <- beta1
    BETA2[mkeep, ] <- beta2
    THETA1[, , mkeep] <- theta1
    THETA2[, , mkeep] <- theta2
    PHI1[, , mkeep] <- phi1
    PHI2[, , mkeep] <- phi2
    OMEGA1[, , mkeep] <- omega1
    OMEGA2[mkeep, ] <- omega2
  } 
  
  #トピック割当はバーンイン期間を超えたら格納する
  if(rp%%keep==0 & rp >= burnin){
    FLAG1 <- FLAG1 + flag1; FLAG2 <- FLAG2 + flag2
    SEG11 <- SEG11 + Zi11; SEG12 <- SEG12 + Zi12
    SEG21 <- SEG21 + Zi21
  }
  
  if(rp%%disp==0){
    #対数尤度を計算
    LL1 <- sum(log((1-flag_rate1)*topic_par11 + flag_rate1*topic_par12))
    LL2 <- sum(log(flag_rate2*topic_par21 + (1-flag_rate2)*topic_par22))
    
    #サンプリング結果を確認
    print(rp)
    print(c(LL1+LL2, LL1, LL2, LLbest, LLbest1, LLbest2, LLst))
    print(round(c(mean(beta1), mean(beta2), mean(betat1), mean(betat2)), 3))
    print(round(cbind(phi1[, c(1:5, (v11+1):(v11+5))], phit1[, c(1:5, (v11+1):(v11+5))]), 3))
    print(round(rbind(omega2[c(1:10, (v21+1):(v21+10))], omegat2[c(1:10, (v21+1):(v21+10))]), 3))
  }
}

####サンプリング結果の要約と可視化####
burnin <- 2000/keep
RS <- R/keep

##サンプリング結果の可視化
#トピック分布の可視化
matplot(t(THETA1[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA1[10, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA1[100, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA1[1000, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA2[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA2[10, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA2[100, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA2[1000, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")

#単語分布の可視化
matplot(t(PHI1[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI1[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI2[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI2[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(OMEGA1[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(OMEGA1[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(OMEGA2, type="l", xlab="サンプリング回数", ylab="パラメータ")


##サンプリング結果の事後分布
#トピック生成過程の事後平均
round(cbind(colMeans(BETA1[burnin:RS, ]), betat1), 3)
round(cbind(colMeans(BETA2[burnin:RS, ]), betat2), 3)

#トピック分布の事後平均
round(cbind(apply(THETA1[, , burnin:RS], c(1, 2), mean), thetat1), 3)
round(cbind(apply(THETA2[, , burnin:RS], c(1, 2), mean), thetat2), 3)

#単語分布の事後平均
round(cbind(t(apply(PHI1[, , burnin:RS], c(1, 2), mean)), t(phit1)), 3)
round(cbind(t(apply(PHI2[, , burnin:RS], c(1, 2), mean)), t(phit2)), 3)
round(cbind(t(apply(OMEGA1[, , burnin:RS], c(1, 2), mean)), t(omegat1)), 3)
round(cbind(colMeans(OMEGA2[burnin:RS, ]), omegat2), 3)


##潜在変数のサンプリング結果の事後分布
seg11_rate <- SEG11 / rowSums(SEG11); seg12_rate <- SEG12 / rowSums(SEG12)
seg21_rate <- SEG21 / rowSums(SEG21)
seg11_rate[is.nan(seg11_rate)] <- 0; seg12_rate[is.nan(seg12_rate)] <- 0
seg21_rate[is.nan(seg21_rate)] <- 0

#トピック割当結果を比較
round(cbind(rowSums(SEG11), seg11_rate, Z11), 3)
round(cbind(rowSums(SEG12), seg12_rate, Z12), 3)
round(cbind(rowSums(SEG21), seg21_rate, Z21), 3)


