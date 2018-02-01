#####Multinomial Topic Regression Model####
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

#set.seed(5723)

####データの発生####
dataset <- 2   #データセット数
k <- 15   #トピック数
hh <- 2000   #ユーザー数
page <- 400   #ページ数
item <- 1000   #アイテム数
w1 <- rpois(hh, rgamma(hh, 35, 0.3))   #1人あたいのページ閲覧数
w2 <- rpois(hh, rgamma(hh, 32, 0.3))   #1人あたりのアイテム購買数
f1 <- sum(w1)   #総ページ閲覧数
f2 <- sum(w2)   #総アイテム購買数


#IDの設定
u_id1 <- rep(1:hh, w2)
u_id2 <- rep(1:hh, w1)
t_id1 <- c()
t_id2 <- c()
for(i in 1:hh){
  t_id1 <- c(t_id1, 1:w1[i])
  t_id2 <- c(t_id2, 1:w2[i])
}

##デザイン行列の設定
Data <- cbind(rep(1, hh*dataset), rep(0:(dataset-1), rep(hh, dataset)))
index_logit <- rep(1:hh, dataset)
index_data <- matrix(1:(hh*dataset), nrow=hh, ncol=dataset)

##パラメータの設定
for(rp in 1:100){
  print(rp)
  #トピック分布のパラメータ
  alpha01 <- alphat01 <- cbind(mvrnorm(hh, rep(0, k-1), diag(3.5, k-1)), 0)
  alpha02 <- alphat02 <- cbind(mvrnorm(hh, rep(0, k-1), diag(0.75, k-1)), 0)
  alpha11 <- rep(0.15, page)
  alpha12 <- rep(0.1, item)

  #パラメータを生成
  logit <- Data[, 1]*alpha01[index_logit, ] + Data[, 2]*alpha02[index_logit, ]   #ロジットの計算
  theta1 <- thetat1 <- exp(logit[index_data[, 1], ]) / rowSums(exp(logit[index_data[, 1], ]))
  theta2 <- thetat2 <- exp(logit[index_data[, 2], ]) / rowSums(exp(logit[index_data[, 2], ]))
  phi0 <- t(extraDistr::rdirichlet(page, rep(0.01, k))) * 
                (matrix(extraDistr::rdirichlet(1, rep(2.0, page)), nrow=k, ncol=page, byrow=T))
  phi <- phit <- phi0 / rowSums(phi0)
  gamma0 <- t(extraDistr::rdirichlet(item, rep(0.01, k))) * 
                (matrix(extraDistr::rdirichlet(1, rep(2.0, item)), nrow=k, ncol=item, byrow=T))
  gamma <- gammat <- gamma0 / rowSums(gamma0)
  
  ##モデルにもとづき単語を生成する
  WX1 <- matrix(0, nrow=hh, ncol=page)
  WX2 <- matrix(0, nrow=hh, ncol=item)
  Z1_list <- list()
  Z2_list <- list()
  wd1_list <- list()
  wd2_list <- list()
  
  for(i in 1:hh){
    #多項分布からトピックを生成
    z1 <- rmnom(w1[i], 1, theta1[i, ])
    z1_vec <- as.numeric(z1 %*% 1:k)
    z2 <- rmnom(w2[i], 1, theta2[i, ])
    z2_vec <- as.numeric(z2 %*% 1:k)
    
    #トピック分布から履歴を生成
    data1 <- rmnom(w1[i], 1, phi[z1_vec, ])
    data2 <- rmnom(w2[i], 1, gamma[z2_vec, ])
    
    #データを格納
    WX1[i, ] <- colSums(data1)
    WX2[i, ] <- colSums(data2)
    Z1_list[[i]] <- z1
    Z2_list[[i]] <- z2
    wd1_list[[i]] <- as.numeric(data1 %*% 1:page)
    wd2_list[[i]] <- as.numeric(data2 %*% 1:item)
  }
  if(min(colSums(WX1)) > 0 & min(colSums(WX2)) > 0) break
}

#リストを変換
wd1 <- unlist(wd1_list)
wd2 <- unlist(wd2_list)
z1 <- unlist(Z1_list)
z2 <- unlist(Z2_list)
sparse_data1 <- as(WX1, "CsparseMatrix")
sparse_data2 <- as(WX2, "CsparseMatrix")


##インデックスを作成
user_list1 <- list()
user_list2 <- list()
page_list <- list()
item_list <- list()
for(i in 1:hh){
  user_list1[[i]] <- which(u_id1==i)
  user_list2[[i]] <- which(u_id2==i)
}
for(i in 1:page){page_list[[i]] <- which(wd1==i)}
for(i in 1:item){item_list[[i]] <- which(wd2==i)}


####マルコフ連鎖モンテカルロ法でMultinomial Topic Regression Modelを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #負担係数の格納用
  for(j in 1:k){
    #負担係数を計算
    Bi <- rep(theta[, j], w) * phi[j, wd]   #尤度
    Bur[, j] <- Bi   
  }
  Br <- Bur / rowSums(Bur)   #負担率の計算
  bval <- list(Br=Br, Bur=Bur)
  return(bval)
}

##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10


##初期値を設定
alpha1 <- cbind(mvrnorm(hh, rep(0, k-1), diag(1.0, k-1)), 0)
alpha2 <- cbind(mvrnorm(hh, rep(0, k-1), diag(0.25, k-1)), 0)
logit <- Data[, 1]*alpha1[index_logit, ] + Data[, 2]*alpha2[index_logit, ]   #ロジットの計算
theta1 <- exp(logit[index_data[, 1], ]) / rowSums(exp(logit[index_data[, 1], ]))
theta2 <- exp(logit[index_data[, 2], ]) / rowSums(exp(logit[index_data[, 2], ]))
phi0 <- extraDistr::rdirichlet(k, colSums(WX1)/sum(WX1)*10) + 0.0001
phi <- phi0 / rowSums(phi0)
gamma0 <- extraDistr::rdirichlet(k, colSums(WX2)/sum(WX2)*10) + 0.0001
gamma <- gamma0 / rowSums(gamma0)

##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- 0.5
alpha02 <- 0.5

#変量効果の事前分布
nu <- k-1   #逆ウィシャート分布の自由度
V <- nu * diag(rep(1, (dataset*k)-dataset))


##パラメータの格納用配列
ALPHA1 <- array(0, dim=c(hh, k-1, R/keep))
ALPHA2 <- array(0, dim=c(hh, k-1, R/keep))
SIGMA <- matrix(0, nrow=R/keep, ncol=(dataset*k)-dataset) 
PHI <- array(0, dim=c(k, page, R/keep))
GAMMA <- array(0, dim=c(k, item, R/keep))
SEG1 <- matrix(0, nrow=f1, ncol=k)
SEG2 <- matrix(0, nrow=f2, ncol=k)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"

#対数尤度の基準値
LLst <- sum(WX1 %*% log(colSums(WX1)/f1) + WX2 %*% log(colSums(WX2)/f2))

##パラメータの真値
alpha1 <- alpha01
alpha2 <- alpha02
theta1 <- thetat1
theta2 <- thetat2
phi <- phit
gamma <- gammat


####ギブスサンプリングでHTMモデルのパラメータをサンプリング####
for(rp in 1:R){
  
  ##データセットごとにトピックを生成
  #トピック尤度を計算
  page_par <- burden_fr(theta1, phi, wd1, w1, k)
  item_par <- burden_fr(theta2, gamma, wd2, w2, k)
  page_rate <- page_par$Br
  item_rate <- item_par$Br
  
  #トピックの割当確率からトピックを生成
  Zi1 <- rmnom(f1, 1, page_rate)
  Zi2 <- rmnom(f2, 1, item_rate)
  z1_vec <- as.numeric(Zi1 %*% 1:k)
  z2_vec <- as.numeric(Zi2 %*% 1:k)
  
  
  ##HMMで文章単位のセグメントを生成
  #文章単位でのトピック頻度行列を作成
  HMM_data <- matrix(0, nrow=a, ncol=k2)
  for(i in 1:a){
    HMM_data[i, ] <- vec_list[[i]] %*% Zi2[s_list[[i]], , drop=FALSE]
  }
  
  #潜在変数ごとに尤度を推定
  theta_log <- log(t(theta3))
  LLi0 <- HMM_data %*% theta_log   #対数尤度
  LLi_max <- rowMaxs(LLi0)
  LLi <- exp(LLi0 - LLi_max)   #尤度に変換
  
  #セグメント割当確率の推定とセグメントの生成
  z_rate1 <- matrix(0, nrow=a, ncol=k1)
  Zi1 <- matrix(0, nrow=a, ncol=k1)
  z1_vec <- rep(0, a)
  rf02 <- matrix(0, nrow=k1, ncol=k1) 
  
  for(j in 1:max_time){
    if(j==1){
      #セグメントの割当確率
      LLs <- matrix(theta1, nrow=length(index_t11), ncol=k1, byrow=T) * LLi[index_t11, ]   #重み付き尤度
      z_rate1[index_t11, ] <- LLs / rowSums(LLs)   #割当確率
      
      #多項分布よりセグメントを生成
      Zi1[index_t11, ] <- rmnom(length(index_t11), 1, z_rate1[index_t11, ])
      z1_vec[index_t11] <- as.numeric(Zi1[index_t11, ] %*% 1:k1)
      
      #混合率のパラメータを更新
      rf01 <- colSums(Zi1[index_t11, ])
      
    } else {
      
      #セグメントの割当確率
      index <- index_t22[[j]]
      LLs <- theta2[z1_vec[index_t21[[j]]], , drop=FALSE] * LLi[index, , drop=FALSE]   #重み付き尤度
      z_rate1[index, ] <- LLs / rowSums(LLs)   #割当確率
      
      #多項分布よりセグメントを生成
      Zi1[index, ] <- rmnom(length(index), 1, z_rate1[index, ])
      z1_vec[index] <- as.numeric(Zi1[index, ] %*% 1:k1)
      
      #混合率のパラメータを更新
      rf02 <- rf02 + t(Zi1[index_t21[[j]], , drop=FALSE]) %*% Zi1[index, , drop=FALSE]   #マルコフ推移
    }
  }
  
  ##パラメータをサンプリング
  #ディクレリ分布からHMMの混合率をサンプリング
  rf11 <- colSums(Zi1[index_t11, ]) + beta01
  rf12 <- rf02 + alpha01
  theta1 <- extraDistr::rdirichlet(1, rf11)
  theta2 <- extraDistr::rdirichlet(k1, rf12)
  
  #トピック分布のパラメータをサンプリング
  wf0 <- matrix(0, nrow=k1, ncol=k2)
  for(j in 1:k1){
    wf0[j, ] <- colSums(HMM_data * Zi1[, j])
  }
  wf <- wf0 + beta02
  theta3 <- extraDistr::rdirichlet(k1, wf)
  
  #単語分布phiをサンプリング
  vf0 <- matrix(0, nrow=k2, ncol=v)
  for(j in 1:v){
    vf0[, j] <- colSums(Zi2[word_list[[j]], , drop=FALSE])
  }
  vf <- vf0 + alpha02
  phi <- extraDistr::rdirichlet(k2, vf)
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA1[mkeep, ] <- theta1
    THETA2[, , mkeep] <- theta2
    THETA3[, , mkeep] <- theta3
    PHI[, , mkeep] <- phi
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(mkeep >= burnin & rp%%keep==0){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
    }
    
    #サンプリング結果を確認
    if(rp%%disp==0){
      print(rp)
      print(c(sum(log(rowSums(word_par))), LLst))
      print(round(rbind(theta1, thetat1), 3))
      print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
    }
  }
}


####サンプリング結果の可視化と要約####
burnin <- 2000/keep   #バーンイン期間
RS <- R/keep

##サンプリング結果の可視化
#HMMの初期分布のサンプリング結果
matplot(THETA1, type="l", xlab="サンプリング数", ylab="パラメータ")

#HMMのパラメータのサンプリング結果
matplot(t(THETA2[1, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA2[5, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA2[15, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")

#文書のトピック分布のサンプリング結果
matplot(t(THETA3[1, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA3[5, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA3[15, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")

#単語の出現確率のサンプリング結果
matplot(t(PHI[, 1, ]), type="l", ylab="パラメータ", main="トピック1の単語の出現率のサンプリング結果")
matplot(t(PHI[, 100, ]), type="l", ylab="パラメータ", main="トピック2の単語の出現率のサンプリング結果")
matplot(t(PHI[, 200, ]), type="l", ylab="パラメータ", main="トピック2の単語の出現率のサンプリング結果")
matplot(t(PHI[, 300, ]), type="l", ylab="パラメータ", main="トピック2の単語の出現率のサンプリング結果")
matplot(t(PHI[, 400, ]), type="l", ylab="パラメータ", main="トピック3の単語の出現率のサンプリング結果")
matplot(t(PHI[, 500, ]), type="l", ylab="パラメータ", main="トピック4の単語の出現率のサンプリング結果")


##サンプリング結果の要約推定量
#トピック分布の事後推定量
topic_mu <- apply(THETA[, , burnin:(R/keep)], c(1, 2), mean)   #トピック分布の事後平均
round(cbind(topic_mu, thetat), 3)
round(topic_sd <- apply(THETA[, , burnin:(R/keep)], c(1, 2), sd), 3)   #トピック分布の事後標準偏差

#単語出現確率の事後推定量
word_mu <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #単語の出現率の事後平均
word <- round(t(rbind(word_mu, phit)), 3)
colnames(word) <- 1:ncol(word)
word

##トピックの事後分布の要約
round(cbind(z1, seg1_mu <- SEG1 / length(burnin:RS)), 3)
round(cbind(z2, seg2_mu <- SEG2 / rowSums(SEG2)), 3)







