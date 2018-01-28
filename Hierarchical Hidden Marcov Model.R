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
    phi[, , j] <- phit[, , j] <- phi0
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
alpha01 <- 0.1
beta01 <- 1


##サンプリング結果の保存用配列
THETA1 <- matrix(0, nrow=R/keep, ncol=k1)
THETA2 <- array(0, dim=c(k1, k1, R/keep))
THETA3 <- array(0, dim=c(k1, k2, R/keep))
THETA4 <- array(0, dim=c(k2, k2, k1, R/keep))
PHI <- array(0, dim=c(k2, item, R/keep))
SEG1 <- rep(0, hhpt)
SEG2 <- matrix(0, nrow=hhpt, ncol=k1)
SEG3 <- matrix(0, nrow=hhpt, ncol=k2)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"
storage.mode(SEG3) <- "integer"


##MCMC推定用配列
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


####MCMCでパラメータをサンプリング####

##1期目の購買期間の潜在状態を生成
#パターンごとに対数尤度を計算
LLi1 <- matrix(0, nrow=hhpt, ncol=k1*k2)
LLm <- matrix(0, nrow=hhpt, ncol=k1*k2)
for(j in 1:k1){
  log_theta3 <- matrix(log(theta3[j, ]), nrow=hhpt, ncol=k2, byrow=T)
  LLm[, index_k2[j, ]] <- as.matrix(const + sparse_data %*% t(log(phi[, , j])))
  LLi1[, index_k2[j, ]] <- log_theta3 + LLm[, index_k2[j, ]]
}

##上位階層の潜在状態を生成
#対数尤度を上位階層の期待尤度に変換
LLt1 <- LLi1[index_t11, ]
par0 <- exp(LLt1 - rowMaxs(LLt1))
par1 <- matrix(0, nrow=hh, ncol=k1)
for(j in 1:k1){
  par1[, j] <- rowSums(par0[, index_k2[j, ]])
}

#多項分布から上位階層の状態を生成
Zi2 <- matrix(0, nrow=hhpt, ncol=k1)
z2_vec <- rep(0, hhpt)
par_rate1 <- par1 / rowSums(par1)   #潜在変数の割当確率
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
cbind(z3_vec[index_t11], Z3[index_t11])

##2期目以降のマルコフ推移に基づき潜在状態を生成
##潜在状態が切り替わるかどうかを生成
#データの設定
j <- 2
index1 <- index_t21[[j]]
index2 <- index_t22[[j]]
n <- length(index2)
zi2_j <- Zi2[index1, ]
zi3_j <- Zi3[index1, ]
z3_j <- z3_vec[index1]

#期待尤度を計算
LLz1 <- matrix(0, nrow=n, ncol=k2*2)
LLi2 <- matrix(0, nrow=hhpt, ncol=k2) 
LLo <- matrix(0, nrow=n, ncol=k2)
theta4_par <- matrix(0, nrow=n, ncol=k2) 

for(l in 1:k1){
  LLz1 <- LLz1 + LLi1[index2, as.numeric(index_k2[-l, ])] * (1-zi2_j[, l]) 
  LLi2[index2, ] <- LLm[index2, index_k2[l, ]] * zi2_j[, l]
  LLo <- LLo + LLi2[index2, ] * zi3_j
  theta4_par <- theta4_par + theta4[z3_j, , l] * zi2_j[, l]
}
LLt2 <- rowSums(log(theta4_par) + LLo)


#上位階層の潜在状態切換えを生成
Zi1 <- rep(0, hhpt)
LLz0 <- cbind(LLz1, LLt2)
max_z <- rowMaxs(LLz0)
LLz <- exp(LLz0 - max_z)
beta_rate0 <- LLz / rowSums(LLz)   #潜在変数の割当確率
beta_rate <- rowSums(beta_rate0[, 1:(k2*2)])
Zi1[index2] <- rbinom(n, 1, beta_rate)


cbind(beta_rate, Z1[index2])


