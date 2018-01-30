#####Hierarchical Hidden Marcov Model#####
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

#set.seed(5723)

####データの発生####
##データの設定
k1 <- 3   #混合数
k2 <- 10   #HMMの状態数
hh <- 3000   #ユーザー数
item <- 1000   #総アイテム数
pt <- rpois(hh, rgamma(hh, 13, 0.5))   #ユーザーごとの購買機会数
pt[pt < 5] <- ceiling(runif(sum(pt < 5), 5, 10))
hhpt <- sum(pt)   #総レコード数
w <- rpois(hhpt, rgamma(hhpt, 10, 0.75))   #期間ごとの購買数
w[w < 5] <- ceiling(runif(sum(w < 5), 5, 10))
f <- sum(w)   #総購買数

#IDの設定
u_id <- rep(1:hh, pt)
t_id <- c()
for(i in 1:hh){t_id <- c(t_id, 1:pt[i])}

#インデックスを作成
id_list <- list()
for(i in 1:hh){id_list[[i]] <- which(u_id==i)}


##パラメータの設定
#ディレクリ分布のパラメータを設定
alpha01 <- rep(1.5, k1)
alpha02 <- matrix(1.5, nrow=k1, ncol=k1)
alpha03 <- rep(1.5, k2)
alpha04 <- matrix(0.2, nrow=k2, ncol=k2)
diag(alpha04) <- 2.0

#ディレクリ分布からパラメータを生成
beta1 <- betat1 <- rbeta(hh, 2.5, 7.0)
theta1 <- thetat1 <- extraDistr::rdirichlet(1, alpha01)
theta0 <- extraDistr::rdirichlet(k1, alpha02)
diag(theta0) <- 0

theta2 <- thetat2 <- theta0/rowSums(theta0)
theta3 <- thetat3 <- matrix(0, nrow=k1, ncol=k2)
theta4 <- thetat4 <- array(0, dim=c(k2, k2, k1))
phi <- phit <- array(0, dim=c(k2, item, k1))
for(i in 1:1000){
  for(j in 1:k1){
    theta3[j, ] <- thetat3[j, ] <- as.numeric(extraDistr::rdirichlet(1, alpha03))
    theta4[, , j] <- thetat4[, , j] <- extraDistr::rdirichlet(k2, alpha04)
    phi0 <- t(extraDistr::rdirichlet(item, rep(0.015, k2))) * 
                    (matrix(extraDistr::rdirichlet(1, rep(2.0, item)), nrow=k2, ncol=item, byrow=T))
    phi[, , j] <- phit[, , j] <- phi0 / rowSums(phi0)
    if(sum(phi==0)==0) break
  }
}


##HHMMモデルに基づきデータを生成
z1_list <- list()
z2_list <- list()
z3_list <- list()
Data <- matrix(0, nrow=hhpt, ncol=item)

for(i in 1:hh){
  #切換え変数を生成
  z1_vec <- rbinom(pt[i], 1, beta1[i])
  z1_vec[1] <- 1
  
  #多項分布よりマルコフ状態推移を生成
  z2_vec <- rep(0, pt[i])
  z3_vec <- rep(0, pt[i])
  data <- matrix(0, nrow=pt[i], ncol=item)
  
  for(j in 1:pt[i]){
    if(j==1){
      #上位階層の状態推移を生成
      z2 <- rmnom(1, 1, theta1)
      z2_vec[j] <- as.numeric(z2 %*% 1:k1)    
      
      #下位の階層の状態推移を生成
      z3 <- rmnom(1, 1, theta3[z2_vec[j], ])
      z3_vec[j] <- as.numeric(z3 %*% 1:k2)
      
    } else {
      if(z1_vec[j]==1){
        #上位階層の状態推移を生成
        z2 <- rmnom(1, 1, theta2[z2_vec[j-1], ])
        z2_vec[j] <- as.numeric(z2 %*% 1:k1)
        
        #下位の階層の状態推移を生成
        z3 <- rmnom(1, 1, theta3[z2_vec[j], ])
        z3_vec[j] <- as.numeric(z3 %*% 1:k2)
        
      } else {
        #上位階層の状態推移を生成
        z2_vec[j] <- z2_vec[j-1]
        
        #下位の階層の状態推移を生成
        z3 <- rmnom(1, 1, theta4[z3_vec[j-1], , z2_vec[j]])
        z3_vec[j] <- as.numeric(z3 %*% 1:k2)
      }
    }
    #マルコフ状態推移にもとづきアイテム購買を生成
    data[j, ] <- rmnom(1, w[u_id==i & t_id==j], phi[z3_vec[j], , z2_vec[j]])
  }
  
  #生成したデータを格納
  z1_list[[i]] <- z1_vec
  z2_list[[i]] <- z2_vec
  z3_list[[i]] <- z3_vec
  Data[id_list[[i]], ] <- data
}

#リストを変換
Z1 <- unlist(z1_list)
Z2 <- unlist(z2_list)
Z3 <- unlist(z3_list)
storage.mode(Data) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")
sparse_data_T <- t(sparse_data)

####マルコフ連鎖モンテカルロ法でHHMMモデルを推定####
##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##パラメータの真値
beta1 <- betat1
theta1 <- thetat1
theta2 <- thetat2
theta3 <- thetat3
theta4 <- thetat4
phi <- phit
z1_vec <- Z1


##HHMMモデルの初期値を設定
#多項分布の密度関数の対数尤度の定数
const <- lfactorial(w) - rowSums(lfactorial(Data))   

#パラメータの初期値を設定
beta1 <- rep(0.3, hh)
theta1 <- extraDistr::rdirichlet(1, rep(2.0, k1))
theta0 <- extraDistr::rdirichlet(k1, rep(2.0, k1))
diag(theta0) <- 0
theta2 <- theta0/rowSums(theta0)

alpha <- matrix(0.5, nrow=k2, ncol=k2)
diag(alpha) <- 1.5
theta3 <- matrix(0, nrow=k1, ncol=k2)
theta4 <- array(0, dim=c(k2, k2, k1))
phi <- array(0, dim=c(k2, item, k1))
for(j in 1:k1){
  theta3[j, ] <- thetat3[j, ] <- as.numeric(extraDistr::rdirichlet(1, rep(2.0, k2)))
  theta4[, , j] <- thetat4[, , j] <- extraDistr::rdirichlet(k2, alpha)
  phi0 <- extraDistr::rdirichlet(k2, colSums(Data)/sum(Data)*item/10) + 0.0001
  phi[, , j] <- phi0/rowSums(phi0)
}


##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- 0.01
alpha02 <- 0.1
beta01 <- 0.01
beta02 <- 1
beta03 <- 0.25


##サンプリング結果の保存用配列
THETA1 <- matrix(0, nrow=R/keep, ncol=k1)
THETA2 <- array(0, dim=c(k1, k1, R/keep))
THETA3 <- array(0, dim=c(k1, k2, R/keep))
THETA4 <- array(0, dim=c(k2, k2, k1, R/keep))
PHI <- array(0, dim=c(k2, item, k1, R/keep))
SEG1 <- rep(0, hhpt)
SEG2 <- matrix(0, nrow=hhpt, ncol=k1)
SEG3 <- matrix(0, nrow=hhpt, ncol=k2)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"
storage.mode(SEG3) <- "integer"


##MCMC用インデックスを作成
max_time <- max(t_id)
index_t11 <- which(t_id==1)
index_t21 <- list()
index_t22 <- list()
for(j in 2:max_time){
  index_t21[[j]] <- which(t_id==j)-1
  index_t22[[j]] <- which(t_id==j)
}
index_k2 <- matrix(1:(k2*k1), nrow=k1, ncol=k2, byrow=T)


#対数尤度の基準値
LLst <- sum(dmnom(Data, w, colSums(Data)/sum(Data), log=TRUE))

##パラメータの真値
beta1 <- betat1
theta1 <- thetat1
theta2 <- thetat2
theta3 <- thetat3
theta4 <- thetat4
phi <- phit
z1_vec <- Z1

####MCMCでパラメータをサンプリング####
for(rp in 1:R){

  ##パターンごとに対数尤度を計算
  #配列の設定
  rf12 <- matrix(0, nrow=k1, ncol=k1)
  rf22 <- array(0, dim=c(k2, k2, k1))
  Zi1 <- rep(0, hhpt)
  Zi2 <- matrix(0, nrow=hhpt, ncol=k1)
  z2_vec <- rep(0, hhpt)
  LLi1 <- matrix(0, nrow=hhpt, ncol=k1*k2)
  LLm <- matrix(0, nrow=hhpt, ncol=k1*k2)
  
  for(j in 1:k1){
    log_theta3 <- matrix(log(theta3[j, ]), nrow=hhpt, ncol=k2, byrow=T)
    LLm[, index_k2[j, ]] <- as.matrix(const + sparse_data %*% t(log(phi[, , j])))
    LLi1[, index_k2[j, ]] <- log_theta3 + LLm[, index_k2[j, ]]
  }
  
  for(pd in 1:max_time){
    if(pd==1){
      ##1期目の購買期間の潜在状態を生成
      ##上位階層の潜在状態を生成
      #対数尤度を上位階層の期待尤度に変換
      LLt1 <- LLi1[index_t11, ]
      par0 <- exp(LLt1 - rowMaxs(LLt1))
      par1 <- matrix(0, nrow=hh, ncol=k1)
      for(j in 1:k1){
        par1[, j] <- rowSums(par0[, index_k2[j, ]])
      }
      
      #多項分布から上位階層の状態を生成
      r <- matrix(theta1, nrow=hh, ncol=k1, byrow=T)
      par_rate1 <- r * par1 / rowSums(r * par1)   #潜在変数の割当確率
      Zi2[index_t11, ] <- rmnom(hh, 1, par_rate1)   #状態を生成
      z2_vec[index_t11] <- as.numeric(Zi2[index_t11, ] %*% 1:k1)
      
      
      ##下位階層の潜在状態を生成
      #上位階層の状態に応じて尤度を計算
      LLt2 <- matrix(0, nrow=hh, ncol=k2)
      zi2 <- Zi2[index_t11, ]
      for(j in 1:k1){
       LLt2 <- LLt2 + LLt1[, index_k2[j, ]] * Zi2[index_t11, j] * zi2[, j]
      }
      
      #多項分布から下位階層の状態を生成
      Zi3 <- matrix(0, nrow=hhpt, ncol=k2)
      z3_vec <- rep(0, hhpt)
      par2 <- exp(LLt2 - rowMaxs(LLt2))
      par_rate2 <- par2 / rowSums(par2)   #潜在変数の割当確率
      Zi3[index_t11, ] <- rmnom(hh, 1, par_rate2)   #状態を生成
      z3_vec[index_t11] <- as.numeric(Zi3[index_t11, ] %*% 1:k2)
      
    } else {
      
      ##2期目以降のマルコフ推移に基づき潜在状態を生成
      ##潜在状態が切り替わるかどうかを生成
      #データの設定
      index1 <- index_t21[[pd]]
      index2 <- index_t22[[pd]]
      n <- length(index2)
      zi2_j <- Zi2[index1, , drop=FALSE]
      zi3_j <- Zi3[index1, , drop=FALSE]
      z3_j <- z3_vec[index1]
      
      #期待尤度を計算
      LLz1 <- matrix(0, nrow=n, ncol=k1*k2)
      LLi2 <- matrix(0, nrow=n, ncol=k1*k2) 
      LLt2 <- matrix(0, nrow=n, ncol=k2)
      theta4_par <- matrix(0, nrow=n, ncol=k2) 
      
      for(l in 1:k1){
        theta4_par <- theta4_par + theta4[z3_j, , l] * zi2_j[, l]
        LLz1[, index_k2[l, ]] <- LLi1[index2, index_k2[l, ]] * (1-zi2_j[, l])
        LLo <- (log(theta4[z3_vec[index1], , l]) + LLm[index2, index_k2[l, ]]) * zi2_j[, l] 
        LLt2 <- LLt2 + LLo
        LLi2[, index_k2[l, ]] <- LLz1[, index_k2[l, ]] + LLo
      }
      
      #上位階層の潜在状態切換え変数を生成
      max_z <- rowMaxs(LLi2)
      LLz0 <- exp(LLz1 - max_z)
      LLz <- matrix(0, nrow=n, ncol=k1+1)
      for(l in 1:k1){
        LLz[, l] <- rowSums(LLz0[, index_k2[l, ], drop=FALSE])
      }
      
      LLz[, ncol(LLz)] <- rowSums(exp(LLt2 - max_z))
      LLz[, -ncol(LLz)] <- LLz[, -ncol(LLz)] * (1-zi2_j)
      r <- beta1[u_id[index2]]   #混合率
      beta_par <- r * rowSums(LLz[, 1:k1, drop=FALSE])
      beta_rate <- beta_par / (beta_par + (1-r)*LLz[, ncol(LLz), drop=FALSE])   #潜在変数の割当確率
      Zi1[index2] <- rbinom(n, 1, beta_rate)   #二項分布から切換え変数を生成
      index_z1 <- which(Zi1[index2]==1)
      
      
      ##切換え変数を条件として上位階層を生成
      if(length(index_z1)==0){
        Zi2[index2, ] <- Zi2[index1, ]
      } else {
        r <- theta2[(Zi2 %*% 1:k1)[index1], ]
        Zi2[index2, ] <- Zi2[index1, ]   #1期前の潜在状態を繰り越す
        par_rate1 <- r * LLz[, -ncol(LLz), drop=FALSE] / rowSums(r * LLz[, -ncol(LLz), drop=FALSE])   #潜在状態の割当確率
        
        if(nrow(par_rate1)==1){
          Zi2[index2, ] <- rmnom(length(index_z1), 1, par_rate1[index_z1, ])   #潜在状態を生成
          z2_vec[index2] <- as.numeric(Zi2[index2, ] %*% 1:k1)
        } else {
          Zi2[index2, ][index_z1, ] <- rmnom(length(index_z1), 1, par_rate1[index_z1, ])   #潜在状態を生成
          z2_vec[index2][index_z1] <- as.numeric(Zi2[index2, ][index_z1, ] %*% 1:k1)
        }
        #上位階層のマルコフ推移行列を更新
        rf12 <- rf12 + t(Zi2[index1, , drop=FALSE]) %*% (Zi2[index2, , drop=FALSE] * Zi1[index2])
      }
  
      
      ##下位階層の潜在状態を生成
      #上位階層の状態に応じて尤度を計算
      LLo <- matrix(0, nrow=n, ncol=k2)
      for(j in 1:k1){
        LLo <- LLo + (LLi1[index2, , drop=FALSE] * Zi1[index2])[, index_k2[j, ], drop=FALSE] * Zi2[index2, j]
      }
      LLt3 <- LLo + LLt2*(1-Zi1[index2])
      
      #多項分布から下位階層の状態を生成
      par2 <- exp(LLt3 - rowMaxs(LLt3))
      par_rate2 <- par2 / rowSums(par2)   #潜在変数の割当確率
      Zi3[index2, ] <- rmnom(n, 1, par_rate2)   #状態を生成
      z3_vec[index2] <- as.numeric(Zi3[index2, ] %*% 1:k2)
      
      #マルコフ推移行列を更新
      for(j in 1:k1){
        rf22[, , j] <- rf22[, , j] + t(Zi3[index1, , drop=FALSE]) %*% 
                        (Zi3[index2, , drop=FALSE] * Zi2[index2, j] * (1-Zi1[index2]))
      }
    }
  }
  
  ##パラメータの更新
  #ベータ分布から混合率をサンプリング
  par1 <- as.numeric(tapply(Zi1[-index_t11], u_id[-index_t11], sum))
  par2 <- pt-1-par1
  beta1 <- rbeta(hh, par1+alpha01, par2+beta01)
  
  #ディクレリ分布から上位階層の混合率をサンプリング
  wf1 <- colSums(Zi2[index_t11, ]) + beta02
  wf2 <- rf12 + beta02
  theta1 <- extraDistr::rdirichlet(1, wf1)
  theta2 <- extraDistr::rdirichlet(k1, wf2)
  diag(theta2) <- 0
  theta2 <- theta2 / rowSums(theta2)
  
  #ディクレリ分布から下位階層の混合率をサンプリング
  for(j in 1:k1){
    wf3 <- colSums(Zi3[index_t11, ] * Zi2[index_t11, j]) + colSums(Zi3 * Zi1 * Zi2[, j]) + beta02
    wf4 <- rf22[, , j] + beta02
    theta3[j, ] <- extraDistr::rdirichlet(1, wf3)
    theta4[, , j] <- extraDistr::rdirichlet(k2, wf4)
  }
  
  #ディクレリ分布からアイテム購買確率のパラメータをサンプリング
  for(j in 1:k1){
    vf <- t(sparse_data_T %*% (Zi3 * Zi2[, j])) + beta03
    phi[, , j] <- extraDistr::rdirichlet(k2, vf)
  }

  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA1[mkeep, ] <- theta1
    THETA2[, , mkeep] <- theta2
    THETA3[, , mkeep] <- theta3
    THETA4[, , , mkeep] <- theta4
    PHI[, , , mkeep] <- phi
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(mkeep >= burnin & rp%%keep==0){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
      SEG3 <- SEG3 + Zi3
    }
    
    #サンプリング結果を確認
    if(rp%%disp==0){
      print(rp)
      #print(c(sum(log(rowSums(word_par))), LLst))
      print(round(rbind(theta1, thetat1), 3))
      print(round(cbind(theta2, thetat2), 3))
      print(round(rbind(theta3, thetat3), 3))
    }
  }
}

