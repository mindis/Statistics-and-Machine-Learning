#####Four Levels Pachinko Allocation Model#####
options(warn=2)
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

#set.seed(93441)

####データの発生####
##データの設定
L <- 2
k1 <- 4   #上位トピック数
k2 <- 15   #下位トピック数
d <- 2500   #文書数
v <- 1000   #語彙数
w <- rpois(d, rgamma(d, 60, 0.4))   #文書あたりの単語数
f <- sum(w)   #総単語数

#文書IDの設定
d_id <- rep(1:d, w)
a_id <- c()
for(i in 1:d){a_id <- c(a_id, 1:w[i])}


##パラメータの設定
#ディレクリ分布のパラメータを設定
alpha1 <- rep(0.2, k1)
alpha2 <- rep(0.2, k2)
beta1 <- rep(0.15, k2)
beta2 <- rep(0.085, v)

##モデルに基づきデータを生成
for(rp in 1:1000){
  print(rp)
  
  #ディレクリ分布からパラメータを生成
  theta1 <- thetat1 <- extraDistr::rdirichlet(d, alpha1)
  theta2 <- thetat2 <- array(0, dim=c(d, k2, k1))
  for(j in 1:k1){
    theta2[, , j] <- thetat2[, , j] <- extraDistr::rdirichlet(d, alpha2)
  }
  gamma <- gammat <- extraDistr::rdirichlet(k1, beta1)
  phi <- phit <- extraDistr::rdirichlet(k2, beta2)
  
  ##文書ごとにトピックと単語を生成
  Z1_list <- list()
  Z2_list <- list()
  word_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  
  for(i in 1:d){
    #上位トピックを生成
    z1 <- rmnom(w[i], 1, theta1[i, ])
    z1_vec <- as.numeric(z1 %*% 1:k1)
    
    #下位トピックを生成
    z2 <- rmnom(w[i], 1, t(theta2[i, , z1_vec]))
    z2_vec <- as.numeric(z2 %*% 1:k2)
    
    #単語を生成
    word <- rmnom(w[i], 1, phi[z2_vec, ])
    
    #データを格納
    Z1_list[[i]] <- z1
    Z2_list[[i]] <- z2
    word_list[[i]] <- word
    WX[i, ] <- colSums(word)
  }
  if(min(colSums(WX)) > 0){
    break
  }
}

#データを変換
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)
storage.mode(WX) <- "integer"
sparse_data <- as(do.call(rbind, word_list), "CsparseMatrix")
wd <- as.numeric(sparse_data %*% 1:v)
rm(word_list); rm(Z1_list); rm(Z2_list)
gc(); gc()


####マルコフ連鎖モンテカルロ法でPAMを推定####
##単語ごとに尤度と負担率を計算する関数
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
burnin <- 1000/keep
disp <- 10

##インデックスの設定
doc_list <- list()
doc_vec <- list()
for(i in 1:d){
  doc_list[[i]] <- which(d_id==i)
  doc_vec[[i]] <- rep(1, length(doc_list[[i]]))
}
wd_list <- list()
wd_vec <- list()
for(j in 1:v){
  wd_list[[j]] <- which(wd==j)
  wd_vec[[j]] <- rep(1, length(wd_list[[j]]))
}

##事前分布の設定
alpha1 <- alpha2 <- 1
beta1 <- 0.1

##パラメータの真値
theta1 <- thetat1
theta2 <- thetat2
phi <- phit
Zi1 <- Z1 
Zi2 <- Z2

##初期値の設定
theta1 <- extraDistr::rdirichlet(d, rep(2.0, k1))
theta2 <- array(0, dim=c(d, k2, k1))
for(j in 1:k1){
  theta2[, , j] <- extraDistr::rdirichlet(d, rep(2.0, k2))
}
phi <- extraDistr::rdirichlet(k2, rep(2.0, v))
Zi1 <- rmnom(f, 1, theta1[d_id, ])
z1_vec <- as.numeric(Z1 %*% 1:k1)

##パラメータの格納用配列
THETA1 <- array(0, dim=c(d, k1, R/keep))
THETA2 <- array(0, dim=c(d, k2, k1, R/keep))
PHI <- array(0, dim=c(k2, v, R/keep))
SEG1 <- matrix(0, nrow=f, ncol=k1)
SEG2 <- matrix(0, nrow=f, ncol=k2)


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){

  ##下位トピックをサンプリング
  #上位トピック割当に基づき下位トピック分布を設定
  theta_s <- matrix(0, nrow=f, ncol=k2)
  for(j in 1:k1){
    theta_s <- theta_s + theta2[d_id, , j] * Zi1[, j]
  }
  
  #トピック分布の尤度と負担率を推定
  word_par <- theta_s * t(phi)[wd, ]
  word_rate <- word_par / rowSums(word_par)
  
  #多項分布よりトピックをサンプリング
  Zi2 <- rmnom(f, 1, word_rate)
  Zi2_T <- t(Zi2)
  z2_vec <- as.numeric(Zi2 %*% 1:k2)
  
  
  ##上位トピックをサンプリング
  #下位トピック割当に基づき上位のトピック出現分布を設定
  theta_s <- matrix(0, nrow=f, ncol=k1)
  for(j in 1:k2){
    theta_s <- theta_s + theta2[d_id, j, ] * Zi2[, j] 
  }
  
  #トピック分布の尤度と負担率を推定
  topic_par <- theta1[d_id, ] * theta_s
  topic_rate <- topic_par / rowSums(topic_par)
  
  #多項分布よりトピックをサンプリング
  Zi1 <- rmnom(f, 1, topic_rate)
  Zi1_T <- t(Zi1)
  z1_vec <- as.numeric(Zi1 %*% 1:k1)
  
  
  ##パラメータをサンプリング
  ##上位トピックのパラメータをサンプリング
  #ディクレリ分布のパラメータ
  wsum01 <- matrix(0, nrow=d, ncol=k1)
  for(i in 1:d){
    wsum01[i, ] <- Zi1_T[, doc_list[[i]]] %*% doc_vec[[i]]
  }
  wsum1 <- wsum01 + alpha1 
  theta1 <- extraDistr::rdirichlet(d, wsum1)   #ディクレリ分布からトピック分布をサンプリング
  
  ##下位トピックのパラメータをサンプリング
  #ディクレリ分布のパラメータ
  for(j in 1:k1){
    wsum02 <- matrix(0, nrow=d, ncol=k2)
    for(i in 1:d){
      Zi2_W <- Zi2_T[, doc_list[[i]]] * matrix(Zi1_T[j, doc_list[[i]]], nrow=k2, ncol=length(doc_list[[i]]), byrow=T)
      wsum02[i, ] <- Zi2_W %*% doc_vec[[i]]
    }
    wsum2 <- wsum02 + alpha2
    theta2[, , j] <- extraDistr::rdirichlet(d, wsum2)   #ディクレリ分布からトピック分布をサンプリング
  }
  
  ##単語分布のパラメータをサンプリング
  #ディクレリ分布のパラメータ
  vsum0 <- matrix(0, nrow=k2, ncol=v)
  for(j in 1:v){
    vsum0[, j] <- Zi2_T[, wd_list[[j]], drop=FALSE] %*% wd_vec[[j]] 
  }
  vsum <- vsum0 + beta1
  phi <- extraDistr::rdirichlet(k2, vsum)   #ディクレリ分布から単語分布をサンプリング
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA1[, , mkeep] <- theta1
    THETA2[, , , mkeep] <- theta2
    PHI[, , mkeep] <- phi
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(rp%%keep==0 & rp >= burnin){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
    }
    
    if(rp%%disp==0){
      #サンプリング結果を確認
      print(rp)
      print(c(sum(log(rowSums(word_par))), sum(log(rowSums(topic_par)))))
      print(round(cbind(theta1[1:10, ], thetat1[1:10, ]), 3))
      print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
    }
  }
}

####サンプリング結果の可視化と要約####
burnin <- 1000/keep
RS <- R/keep

##サンプリング結果の可視化
#上位トピックのトピック分布のサンプリング結果
matplot(t(THETA1[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="上位トピックのサンプリング結果")
matplot(t(THETA1[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="上位トピックのサンプリング結果")
matplot(t(THETA1[10, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="上位トピックのサンプリング結果")
matplot(t(THETA1[15, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="上位トピックのサンプリング結果")
matplot(t(THETA1[20, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="上位トピックのサンプリング結果")

#下位トピックのトピック分布のサンプリング結果
matplot(t(THETA2[1, , 1, ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="下位トピックのサンプリング結果")
matplot(t(THETA2[1, , 2, ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="下位トピックのサンプリング結果")
matplot(t(THETA2[1, , 3, ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="下位トピックのサンプリング結果")
matplot(t(THETA2[1, , 4, ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="下位トピックのサンプリング結果")
matplot(t(THETA2[10, , 1, ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="下位トピックのサンプリング結果")
matplot(t(THETA2[10, , 2, ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="下位トピックのサンプリング結果")
matplot(t(THETA2[10, , 3, ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="下位トピックのサンプリング結果")
matplot(t(THETA2[10, , 4, ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="下位トピックのサンプリング結果")

#単語分布のサンプリング結果
matplot(t(PHI[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="単語分布のサンプリング結果")
matplot(t(PHI[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="単語分布のサンプリング結果")
matplot(t(PHI[10, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="単語分布のサンプリング結果")
matplot(t(PHI[15, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="単語分布のサンプリング結果")

##事後分布の要約統計量
#トピック分布の事後平均
topic_mu1 <- apply(THETA1[, , burnin:RS], c(1, 2), mean)
topic_mu2 <- array(0, dim=c(d, k2, k1))
for(j in 1:k1){
  topic_mu2[, , j] <- apply(THETA2[, , j, burnin:RS], c(1, 2), mean)
}
round(topic_mu1, 3)
round(topic_mu2[, , 1], 3)

#単語分布の事後平均
round(phi_mu <- t(apply(PHI[, , burnin:RS], c(1, 2), mean)), 3)
round(cbind(phi_mu, t(phit)), 3)

#トピック割当の事後分布
topic_rate1 <- SEG1 / rowSums(SEG1) 
topic_allocation1 <- apply(topic_rate1, 1, which.max)
round(data.frame(真値=Z1 %*% 1:k1, 推定=topic_allocation1, z=topic_rate1), 3)

topic_rate2 <- SEG2 / rowSums(SEG2)
topic_allocation2 <- apply(topic_rate2, 1, which.max)
round(data.frame(真値=Z2 %*% 1:k2, 推定=topic_allocation2, z=topic_rate2), 3)
