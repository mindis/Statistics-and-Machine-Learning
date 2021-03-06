#####ガウス過程連続空間トピックモデル#####
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
gc()
#set.seed(5723)

####データの発生####
k <- 10   #トピック数
d <- 2000   #ユーザー数
v <- 500   #語彙数
w <- rpois(d, rgamma(d, 45, 0.3))   #1人あたりのページ閲覧数
f <- sum(w)   #総語彙数

#IDを設定
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}


##パラメータの設定
G0 <- GT0 <- extraDistr::rdirichlet(1, rep(2.0, v))   #単語出現率
u <- ut <- mvrnorm(d, rep(0, k), diag(0.75, k))
phi <- phit <- mvrnorm(v, rep(0, k), diag(0.75, k))


##データの生成
word_list <- list()
word_vec_list <- list()
WX <- matrix(0, nrow=d, ncol=v)

for(i in 1:d){

  #ディクレリ-多項分布から単語を生成
  alpha <- G0 * exp(u[i, ] %*% t(phi))
  words <- extraDistr::rdirmnom(w[i], 1, alpha)
  words_vec <- as.numeric(words %*% 1:v)
  
  #生成した単語を格納
  WX[i, ] <- colSums(words)
  word_list[[i]] <- words
  word_vec_list[[i]] <- words_vec
}

#リストを変換
WX_T <- t(WX)
word_vec <- unlist(word_vec_list)
word_data <- do.call(rbind, word_list)
sparse_data <- as(word_data, "CsparseMatrix")
storage.mode(WX) <- "integer"
storage.mode(WX) <- "integer"
rm(word_data); rm(word_list)


##インデックスを作成
dw_list <- list()
for(j in 1:v){
  index <- which(word_vec==j)
  dw_list[[j]] <- d_id[index]
}


####マルコフ連鎖モンテカルロ法で連続空間トピックモデルを推定####
##アルゴリズムの設定
R <- 10000
keep <- 5 
iter <- 0
burnin <- 2500/keep
disp <- 10

#データの設定
rej2 <- rep(0, v)
WX_vec <- as.numeric(WX)
WX_T_vec  <- as.numeric(t(WX))

#インデックスを設定
index_nzeros1 <- which(as.numeric(WX) > 0)
index_nzeros2 <- which(as.numeric(WX_T) > 0)
index_word <- list()
for(j in 1:v){
  index_word[[j]] <- which(WX[, j] > 0)
}

##事前分布の設定
cov <- diag(k)
inv_cov2 <- inv_cov1 <- solve(cov)
mu <- rep(0, k)
sigma <- diag(k)


##パラメータの真値
G0 <- GT0
u <- ut
phi <- phit

##パラメータの初期値
G0 <- colSums(WX) / sum(WX)
K0 <- rowSums(WX) / sum(WX)
u <- mvrnorm(d, rep(0, k), diag(0.01, k))
phi <- mvrnorm(v, rep(0, k), diag(0.01, k))


##パラメータの保存用配列
logl <- rep(0, R/keep)
U <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(v, k, R/keep))

##対数尤度の基準値
#ユニグラムモデルの対数尤度
par <- colSums(WX)/sum(WX)
LLst <- sum(WX %*% log(par))

#ベストな対数尤度
alpha <- matrix(GT0, nrow=d, ncol=v, byrow=T) * exp(ut %*% t(phit))
LLbest <- sum(lgamma(rowSums(alpha)) - lgamma(rowSums(alpha + WX))) + sum(lgamma(alpha + WX) - lgamma(alpha))



####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##ユーザートピック行列Uをサンプリング
  #新しいパラメータをサンプリング
  u_old <- u
  u_new <- u_old + mvrnorm(d, rep(0, k), diag(0.01, k))
  alphan <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u_new %*% t(phi))
  alphad <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u_old %*% t(phi))
  
  #Polya分布のパラメータを計算
  dirn_vec <- dird_vec <- rep(0, d*v)
  alphan_vec <- as.numeric(alphan)
  alphad_vec <- as.numeric(alphad)
  dirn_vec[index_nzeros1] <- lgamma(alphan_vec[index_nzeros1] + WX_vec[index_nzeros1]) - lgamma(alphan_vec[index_nzeros1])
  dird_vec[index_nzeros1] <- lgamma(alphad_vec[index_nzeros1] + WX_vec[index_nzeros1]) - lgamma(alphad_vec[index_nzeros1])
  
  #対数尤度と対数事前分布を計算
  dir_new <- lgamma(rowSums(alphan)) - lgamma(rowSums(alphan + WX))
  dir_old <- lgamma(rowSums(alphad)) - lgamma(rowSums(alphad + WX))
  lognew1 <- dir_new + rowSums(matrix(dirn_vec, nrow=d, ncol=v))
  logold1 <- dir_old + rowSums(matrix(dird_vec, nrow=d, ncol=v))
  logpnew1 <- -0.5 * rowSums(u_new %*% inv_cov1 * u_new)
  logpold1 <- -0.5 * rowSums(u_old %*% inv_cov1 * u_old)
  
  ##MHサンプリング
  rand <- runif(d)   #一様分布から乱数を発生
  LL_diff <- lognew1 + logpnew1 - logold1 - logpold1
  LL_diff[LL_diff > 100] <- 100
  LLind_diff <- exp(LL_diff)   #採択率を計算
  tau <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #tauの値に基づき新しいbetaを採択するかどうかを決定
  flag <- matrix(((tau >= rand)*1 + (tau < rand)*0), nrow=d, ncol=k)
  rej1 <- mean(flag[, 1])
  u <- flag*u_new + (1-flag)*u_old   #alphaがrandを上回っていたら採択


  ##単語分布のパラメータphiをサンプリング
  index_vec <- index_word[[j]]
  
  #単語ごとにMHサンプリングを実行
  #新しいパラメータをサンプリング
  phid <- phi
  phin <- phid + mvrnorm(v, rep(0, k), diag(0.05, k))
  alphan <- matrix(K0, nrow=v, ncol=d) * exp(phin %*% t(u))
  alphad <- matrix(K0, nrow=v, ncol=d) * exp(phid %*% t(u))
  
  #Polya分布のパラメータを計算
  dirn_vec <- dird_vec <- rep(0, d*v)
  alphan_vec <- as.numeric(alphan)
  alphad_vec <- as.numeric(alphad)
  dirn_vec[index_nzeros2] <- lgamma(alphan_vec[index_nzeros2] + WX_T_vec[index_nzeros2]) - lgamma(alphan_vec[index_nzeros2])
  dird_vec[index_nzeros2] <- lgamma(alphad_vec[index_nzeros2] + WX_T_vec[index_nzeros2]) - lgamma(alphad_vec[index_nzeros2])
  
  #対数尤度と対数事前分布を計算
  dir_new <- lgamma(rowSums(alphan)) - lgamma(rowSums(alphan + WX_T))
  dir_old <- lgamma(rowSums(alphad)) - lgamma(rowSums(alphad + WX_T))
  lognew2 <- dir_new + rowSums(matrix(dirn_vec, nrow=v, ncol=d))
  logold2 <- dir_old + rowSums(matrix(dird_vec, nrow=v, ncol=d))
  logpnew2 <- -0.5 * rowSums(phin %*% inv_cov2 * phin)
  logpold2 <- -0.5 * rowSums(phid %*% inv_cov2 * phid)
  
  
  ##MHサンプリング
  rand <- runif(v)   #一様分布から乱数を発生
  LL_diff <- lognew2 + logpnew2 - logold2 - logpold2
  LL_diff[LL_diff > 100] <- 100
  LLind_diff <- exp(LL_diff)   #採択率を計算
  tau <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #tauの値に基づき新しいbetaを採択するかどうかを決定
  flag <- matrix(((tau >= rand)*1 + (tau < rand)*0), nrow=v, ncol=k)
  rej2 <- mean(flag[, 1])
  phi <- flag*phin + (1-flag)*phid   #alphaがrandを上回っていたら採択

  #パラメータを正規化
  u <- scale(u)
  phi <- scale(phi)
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    logl[mkeep] <- sum(lognew1)
    U[, , mkeep] <- u
    PHI[, , mkeep] <- phi
    
    #サンプリング結果を確認
    if(rp%%disp==0){
      alpha <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u %*% t(phi))
      LL <- sum(lgamma(rowSums(alpha)) - lgamma(rowSums(alpha + WX))) + sum(lgamma(alpha + WX) - lgamma(alpha))
   
      print(rp)
      print(c(LL, LLbest, LLst))
      print(round(c(rej1, mean(rej2)), 3))
      print(round(cbind(u[1:5, ], ut[1:5, ]), 2))
      print(round(cbind(phi[1:5, ], phit[1:5, ]), 2))
    }
  }
}

####サンプリング結果の可視化と要約####
##サンプリング結果をプロット
matplot(t(U[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="トピック分布のサンプリング結果")
matplot(t(U[2, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="トピック分布のサンプリング結果")
matplot(t(U[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="トピック分布のサンプリング結果")
matplot(t(U[4, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="トピック分布のサンプリング結果")
matplot(t(U[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="トピック分布のサンプリング結果")

matplot(t(PHI[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="トピック分布のサンプリング結果")
matplot(t(PHI[2, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="トピック分布のサンプリング結果")
matplot(t(PHI[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="トピック分布のサンプリング結果")
matplot(t(PHI[4, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="トピック分布のサンプリング結果")
matplot(t(PHI[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="トピック分布のサンプリング結果")

