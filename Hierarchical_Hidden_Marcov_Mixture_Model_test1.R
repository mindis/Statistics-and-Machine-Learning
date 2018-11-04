#####Hierarchical Hidden Marcov Mixture Model#####
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
library(data.table)
library(ggplot2)

'%!in%' <- function(a,b) ! a %in% b

#set.seed(5723)

####データの発生####
##データの設定
k1 <- 4   #文脈数
k2 <- 7   #状態数
k3 <- 10   #トピック数
d <- 3000   #文書数
s <- rep(300, k1)   #文脈ごとの語彙数
v <- sum(s)   #総語彙数
v1 <- c(1, (cumsum(s)+1)[-k1])
v2 <- cumsum(s)
w <- rpois(d, rgamma(d, 70, 0.5))   #文書ごとの単語数
f <- sum(w)   #総語彙数

##IDの設定
d_id <- rep(1:d, w)
t_id <- as.numeric(unlist(tapply(1:f, d_id, rank)))
index_end <- as.numeric(tapply(1:f, d_id, max))

##パラメータの事前分布の設定
#マルコフ推移確率のディレクリ分布の事前分布を設定
alpha11 <- rep(5.0, k1)
alpha12 <- matrix(1.5, nrow=k1, ncol=k1)
diag(alpha12) <- 10^-300
alpha21 <- rep(1.5, k2)
alpha22 <- matrix(0.2, nrow=k2, ncol=k2)
diag(alpha22) <- 1.0

#切り替え確率のベータ分布の事前分布を設定
s0 <- 7.5
v0 <- 50.0

#トピック分布のディリクレ分布の事前分布を設定
alpha3 <- array(0.15, dim=c(k2, k3, k1))

#トピックごとの単語分布のディリクレ分布の事前分布を設定
gamma <- matrix(0.0001, nrow=k1, ncol=v)
for(j in 1:k1){
  gamma[j, v1[j]:v2[j]] <- 0.1
}

##すべての単語が生成されるまでループ
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  
  ##事前分布からパラメータを生成
  #ディリクレ分布からマルコフ推移確率を生成
  theta11 <- thetat11 <- as.numeric(extraDistr::rdirichlet(1, alpha11))
  theta12 <- thetat12 <- extraDistr::rdirichlet(k1, alpha12)
  theta21 <- thetat21 <- extraDistr::rdirichlet(k1, alpha21)
  theta22 <- array(0, dim=c(k2, k2, k1))
  for(j in 1:k1){
    theta22[, , j] <- extraDistr::rdirichlet(k2, alpha22)
  }
  thetat22 <- theta22
  
  #ベータ分布から切り替え変数を生成
  beta <- betat <- rbeta(d, s0, v0)
  
  #ディリクレ分布からトピック分布を生成
  theta3 <- array(0, dim=c(k2, k3, k1))
  for(j in 1:k1){
    theta3[, , j] <- extraDistr::rdirichlet(k2, alpha3[, , j])
    
  }
  thetat3 <- theta3
  
  #ディリクレ分布から単語分布を生成
  phi <- array(0, dim=c(k3, v, k1))
  for(j in 1:k1){
    repeat {
      phi[, , j] <- extraDistr::rdirichlet(k3, gamma[j, ])
      if(min(colMaxs(phi[, v1[j]:v2[j], j])) > (k3^2*2)/f) break
    }
  }
  phit <- phi

  ##HHMMMモデルの仮定からデータを生成
  #データの格納用
  z1_list <- list()
  z21_list <- list()
  z22_list <- list()
  z3_list <- list()
  wd_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  
  for(i in 1:d){
    
    #切換え変数を生成
    z1_vec <- rbinom(w[i], 1, beta[i])
    z1_vec[1] <- 1
    
    ##多項分布よりマルコフ状態推移を生成
    z21_vec <- rep(0, w[i])
    z22_vec <- rep(0, w[i])
    z3_vec <- rep(0, w[i])
    word_vec <- rep(0, w[i])
    words <- matrix(0, nrow=w[i], ncol=v)
    
    for(j in 1:w[i]){
      if(j==1){
        #上位階層の状態推移を生成
        z21 <- rmnom(1, 1, theta11)
        z21_vec[j] <- as.numeric(z21 %*% 1:k1)
        
        #下位階層の状態推移を生成
        z22 <- rmnom(1, 1, theta21[z21_vec[j], ])
        z22_vec[j] <- as.numeric(z22 %*% 1:k2)
        
      } else {
        
        if(z1_vec[j]==1){
          #上位階層の状態推移を生成
          z21 <- rmnom(1, 1, theta12[z21_vec[j-1], ])
          z21_vec[j] <- as.numeric(z21 %*% 1:k1)
          
          #下位階層の状態推移を生成
          z22 <- rmnom(1, 1, theta21[z21_vec[j], ])
          z22_vec[j] <- as.numeric(z22 %*% 1:k2)
          
        } else {
          
          #上位階層の状態推移を生成
          z21_vec[j] <- z21_vec[j-1]
          
          #下位階層の状態推移を生成
          z22 <- rmnom(1, 1, theta22[z22_vec[j-1], , z21_vec[j]])
          z22_vec[j] <- as.numeric(z22 %*% 1:k2)
        }
      }
      ##トピックと単語を生成
      #多項分布からトピックを生成
      z3 <- rmnom(1, 1, theta3[z22_vec[j], , z21_vec[j]])
      z3_vec[j] <- as.numeric(z3 %*% 1:k3)
      
      #トピックから単語を生成
      word <- rmnom(1, 1, phi[z3_vec[j], , z21_vec[j]])
      words[j, ] <- word
      word_vec[j] <- as.numeric(word %*% 1:v)
    }
    
    #生成したデータを格納
    z1_list[[i]] <- z1_vec
    z21_list[[i]] <- z21_vec
    z22_list[[i]] <- z22_vec
    z3_list[[i]] <- z3_vec
    wd_list[[i]] <- word_vec
    WX[i, ] <- colSums(words)
  }
  
  if(min(colSums(WX)) > 0){
    break 
  }
}

##データを変換
z1 <- unlist(z1_list)
z21 <- unlist(z21_list)
z22 <- unlist(z22_list)
z3 <- unlist(z3_list)
wd <- unlist(wd_list)
sparse_data <- sparseMatrix(1:f, wd, dims=c(f, v))
sparse_data_T <- t(sparse_data)


####マルコフ連鎖モンテカルロ法でHHMMモデルを推定####
##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##事前分布の設定
#ハイパーパラメータの事前分布
s0 <- 0.5
v0 <- 0.5
alpha011 <- 10.0
alpha012 <- 1.0
alpha02 <- 1.0
beta01 <- 1.0
beta02 <- 0.1


##MCMC用インデックスを作成
#単語順序のインデックス
max_word <- max(t_id)
index_t11 <- which(t_id==1)
index_list_t21 <- index_list_t22 <- list()
n_vec <- c(length(index_t11), rep(0, max_word-1))
for(j in 2:max_word){
  index_list_t21[[j]] <- which(t_id==j)-1
  index_list_t22[[j]] <- which(t_id==j)
  n_vec[j] <- length(index_t22[[j]])
}
index_t21 <- unlist(index_list_t21); index_t22 <- unlist(index_list_t22)
vec_k1 <- rep(1, k1); vec_k2 <- rep(1, k2); vec_k3 <- rep(1, k3)

#文書および単語のインデックス
doc_list <- doc_vec <- list()
wd_list <- wd_vec <- list()
for(i in 1:d){
  doc_list[[i]] <- which(d_id==i)
  doc_vec[[i]] <- rep(1, length(doc_list[[i]]))
}


##真値を設定
#パラメータの真値
const <- lfactorial(w) - rowSums(lfactorial(WX))   #多項分布の密度関数の対数尤度の定数
beta <- betat
theta11 <- thetat11
theta12 <- thetat12
theta21 <- thetat21
theta22 <- thetat22
theta3 <- thetat3
phi <- phit

#トピック割当の真値
Zi1 <- z1
Zi21 <- Z21 <- matrix(as.numeric(table(1:f, z21)), nrow=f, ncol=k1); z21_vec <- z21
Zi22 <- Z22 <- matrix(as.numeric(table(1:f, z22)), nrow=f, ncol=k2); z22_vec <- z22
Zi3 <- Z3 <- matrix(as.numeric(table(1:f, z3)), nrow=f, ncol=k3); z3_vec <- z3


##初期値を設定
#パラメータの初期値を生成
const <- lfactorial(w) - rowSums(lfactorial(WX))   #多項分布の密度関数の対数尤度の定数
beta <- rep(0.5, d)
theta11 <- rep(1/k1, k1)
theta12 <- matrix(1/(k1-1), nrow=k1, ncol=k1); diag(theta12) <- 0
theta21 <- extraDistr::rdirichlet(k1, rep(100.0, k2))
theta22 <- array(0, dim=c(k2, k2, k1))
theta3 <- array(0, dim=c(k2, k3, k1))
phi <- array(0, dim=c(k3, v, k1))
for(j in 1:k1){
  theta22[, , j] <- extraDistr::rdirichlet(k2, rep(100.0, k2))
  theta3[, , j] <- extraDistr::rdirichlet(k2, rep(100.0, k3))
  phi[, , j] <- extraDistr::rdirichlet(k3, rep(100.0, v))
}

#トピック割当の初期値
Zi1 <- rep(0, f); Zi21 <- matrix(0, nrow=f, ncol=k1)
Zi1[index_t11] <- 1; Zi1[-index_t11] <- rbinom(f-length(index_t11), 1, 0.5)
#Zi21[Zi1==1, ] <- rmnom(sum(Zi1), 1, theta11); Zi21 <- Zi21[Zi1==1, ][cumsum(Zi1==1), ]; z21_vec <- as.numeric(Zi21 %*% 1:k1)
Zi21 <- rmnom(f, 1, rep(100.0, k1)); z21_vec <- as.numeric(Zi21 %*% 1:k1)
Zi22 <- rmnom(f, 1, rep(100.0, k2)); z22_vec <- as.numeric(Zi22 %*% 1:k2)
Zi3 <- rmnom(f, 1, rep(100.0, k3)); z3_vec <- as.numeric(Zi3 %*% 1:k3)


##パラメータの保存用配列
BETA <- matrix(0, nrow=R/keep, ncol=d)
THETA11 <- matrix(0, nrow=R/keep, k1)
THETA12 <- array(0, dim=c(k1, k1, R/keep))
THETA21 <- array(0, dim=c(k1, k2, R/keep))
THETA22 <- array(0, dim=c(k2, k2, k1, R/keep))
PHI <- array(0, dim=c(k3, v, k1, R/keep))
SEG1 <- rep(0, f)
SEG21 <- matrix(0, nrow=f, ncol=k1)
SEG22 <- matrix(0, nrow=f, ncol=k2)
SEG3 <- matrix(0, nrow=f, ncol=k3)



#対数尤度の基準値
LLst <- sum(sparse_data %*% log(colSums(sparse_data) / sum(sparse_data)))
Z21_freq <- as.numeric((plyr::count(Z21))$freq)


####MCMCでパラメータをサンプリング####
for(rp in 1:R){
  
  ##単語ごとの尤度を設定
  #トピックごとの期待尤度を推定
  LLt <- array(0, dim=c(f, k2, k1))
  phi_wd <- array(0, dim=c(f, k3, k1)) 
  for(j in 1:k1){
    phi_wd[, , j] <- t(phi[, , j])[wd, ]
    LLt[, , j] <- phi_wd[, , j] %*% t(theta3[, , j]) 
  }
  
  ##上位階層の潜在状態を生成
  #上位階層の混合率
  beta_vec <- beta[d_id]   #切り替え確率をベクトル化
  pi_dt11 <- pi_dt12 <- matrix(1, nrow=f, ncol=k1)
  pi_dt11[index_t11, ] <- matrix(theta11, nrow=d, ncol=k1, byrow=T)   #文書の先頭の混合率
  pi_dt11[index_t22, ] <- beta_vec[index_t21] * theta12[z21_vec[index_t21], ] + (1-beta_vec[index_t21]) * Zi21[index_t21, ]
  pi_dt12[index_t21, ] <- beta_vec[index_t22] * theta12[z21_vec[index_t22], ] + (1-beta_vec[index_t22]) * Zi21[index_t22, ]
  
  
  #対数尤度を上位階層の期待尤度に変換
  par <- matrix(0, nrow=f, ncol=k1)
  for(j in 1:k1){
    par[index_t11, j] <- LLt[index_t11, , j] %*% theta21[j, ]
    par[index_t22, j] <- (LLt[index_t22, , j] %*% t(theta22[, , j]) * Zi22[index_t21, ]) %*% vec_k2
  }
  
  
  #多項分布から上位階層の状態を生成
  Lho1 <- pi_dt11 * pi_dt12 * par
  par_rate1 <- Lho1 / as.numeric(Lho1 %*% vec_k1)   #潜在変数の割当確率
  Zi21 <- rmnom(f, 1, par_rate1)   #多項分布から状態を生成
  z21_vec <- as.numeric(Zi21 %*% 1:k1)
  
  ##ベータ分布から混合率を生成
  #パラメータの格納用配列
  Zi1 <- rep(0, f); Zi1[index_t11] <- 1
  tau1 <- tau2 <- rep(0, d)
  
  for(i in 1:d){
    #文脈の切り替え数をカウント
    Zi1[doc_list[[i]][-1]] <- 1 - as.numeric(z21_vec[doc_list[[i]][-w[i]]]==z21_vec[doc_list[[i]][-1]])
    S <- t(Zi21[doc_list[[i]][-w[i]], ]) %*% Zi21[doc_list[[i]][-1], ]
    diag_S <- sum(diag(S))
    
    #ベータ分布のパラメータ
    tau1[i] <- (w[i]-1) - diag_S + s0
    tau2[i] <- diag_S + v0
  }
  #ベータ分布からパラメータを生成
  beta <- rbeta(d, tau1, tau2)
  index_s0 <- which(Zi1==0); index_s1 <- which(Zi1==1)

  
  ##下位階層の潜在状態を生成
  #下位階層の混合率
  pi_dt21 <- pi_dt22 <- matrix(0, nrow=f, ncol=k2)
  for(j in 1:k1){
    pi_dt21[index_t22, ] <- pi_dt21[index_t22, ] + theta22[z22_vec[index_t21], , j] * Zi21[index_t21, j]   #1単語前の混合率
    pi_dt22[index_t21, ] <- pi_dt22[index_t21, ] + t(theta22[, , j])[z22_vec[index_t22], ] * Zi21[index_t22, j]   #1単語後の混合率
  }
  pi_dt21[index_s1, ] <- theta21[z21_vec[index_s1], ]   #文書と文脈の先頭を設定
  pi_dt22[sort(c(index_s1[index_s1 %!in% index_t11]-1, index_end)), ] <- 1   #文書の末尾を設定
  
  #多項分布から上位階層の状態を生成
  Zi22 <- par_rate2 <- matrix(0, nrow=f, ncol=k2)
  index_z <- list()
  for(j in 1:k1){
    index_z[[j]] <- which(z21_vec==j)
    if(length(index_z[[j]])==0){
      next
    }
    Lho2 <- pi_dt21[index_z[[j]], ] * pi_dt22[index_z[[j]], ] * LLt[index_z[[j]], , j]
    par_rate2[index_z[[j]], ] <- Lho2 / as.numeric(Lho2 %*% vec_k2)   #潜在変数の割当確率
    Zi22[index_z[[j]], ] <- rmnom(length(index_z[[j]]), 1, par_rate2[index_z[[j]], ])   #多項分布から状態を生成
  }
  z22_vec <- as.numeric(Zi22 %*% 1:k2)
  
  
  ##状態に基づきトピックを生成
  #トピック尤度を設定
  Lho_topic <- matrix(0, nrow=f, ncol=k3)
  for(j in 1:k1){
    Lho_topic[index_z[[j]], ] <- t(phi[, , j])[wd[index_z[[j]]], ] * theta3[z22_vec[index_z[[j]]], , j]
  }
  
  #多項分布からトピックを生成
  topic_rate <- Lho_topic / as.numeric(Lho_topic %*% vec_k3)   #トピックの割当確率
  Zi3 <- rmnom(f, 1, topic_rate)
  z3_vec <- as.numeric(Zi3 %*% 1:k3)

  
  ##ディリクレ分布から上位階層の混合率を生成
  #上位階層の混合率のパラメータを設定
  first_vector1 <- colSums2(Zi21[index_t11, ]) + alpha011
  transition_matrix1 <- t(Zi21[index_t21, ]) %*% Zi21[index_t22, ] + alpha012
  diag(transition_matrix1) <- 0 
  
  #ディリクレ分布から混合率を生成
  theta11 <- as.numeric(extraDistr::rdirichlet(1, first_vector1))
  for(j in 1:k1){
    theta12[j, -j] <- as.numeric(extraDistr::rdirichlet(1, transition_matrix1[j, -j]))
  }
  
  
  ##ディリクレ分布から下位階層の混合率を生成
  #下位階層の混合率のパラメータを設定
  index_trans <- index_s0[Zi1[index_s0-1]==0]
  first_vector2 <- t(Zi21[index_s1, ]) %*% Zi22[index_s1, ] + alpha02
  transition_matrix2 <- theta22 <- array(0, dim=c(k2, k2, k1))
  for(j in 1:k1){
    transition_matrix2[, , j] <- t(Zi22[index_trans-1, ]) %*% (Zi21[index_trans, j] * Zi22[index_trans, ]) + alpha02
    theta22[, , j] <- extraDistr::rdirichlet(k2, transition_matrix2[, , j])   
  }

  
  ##トピック分布と単語分布を生成
  for(j in 1:k1){
    #ディリクレ分布からトピック分布を生成
    wsum <- t(Zi21[, j] * Zi22) %*% Zi3 + beta01
    theta3[, , j] <- extraDistr::rdirichlet(k2, wsum)

    #ディリクレ分布から単語分布を生成
    vsum <- as.matrix(t(Zi21[, j] * Zi3) %*% sparse_data + beta02)
    phi[, , j] <- extraDistr::rdirichlet(k3, vsum)
  }

  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
    THETA11[mkeep, ] <- theta11
    THETA12[, , mkeep] <- theta12
    THETA21[, , mkeep] <- theta21
    THETA22[, , , mkeep] <- theta22
    PHI[, , , mkeep] <- phi
  }  

  #トピック割当はバーンイン期間を超えたら格納する
  if(rp%%keep==0 & rp >= burnin){
    SEG1 <- SEG1 + Zi1
    SEG21 <- SEG21 + Zi21
    SEG22 <- SEG22 + Zi22
    SEG3 <- SEG3 + Zi3
  }
  
  if(rp%%disp==0){
    #対数尤度を計算
    LL <- sum(log(rowSums(Lho_topic)))
    
    #サンプリング結果を確認
    print(rp)
    print(c(LL, LLst))
    print(rbind(Z21_freq=colSums(Zi21), Z21_freq))
  }
}
