#####階層ディレクトリLDAモデル#####
options(warn=0)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####データの発生####
##データの設定
k1 <- 15   #一般語のトピック数
k2 <- 20   #ディレクトリのトピック数
dir <- 50   #ディレクトリ数
d <- 3500   #文書数
v1 <- 600    #ディレクトリ構造に関係のある語彙数
v2 <- 600   #ディレクトリ構造に関係のない語彙数
v <- v1 + v2   #総語彙数
w <- rpois(d, rgamma(d, 65, 0.5))   #文書あたりの単語数
f <- sum(w)   #総単語数

##IDの設定
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}

##ディレクトリの割当を設定
dir_freq <- rtpois(d, 1.0, 0, 5)   #文書あたりのディレクトリ数
dir_id <- rep(1:d, dir_freq)   #ディレクトリのid
dir_n <- length(dir_id)
dir_index <- list()
for(i in 1:d){
  dir_index[[i]] <- which(dir_id==i)
}

#ディレクトリの生成
dir_prob <- as.numeric(extraDistr::rdirichlet(1, rep(2.5, dir)))
dir_data <- matrix(0, nrow=dir_n, ncol=dir)

for(i in 1:d){
  repeat{
    x <- rmnom(dir_freq[i], 1, dir_prob)
    if(max(colSums(x))==1){
      index <- dir_index[[i]]
      x <- x[order(as.numeric(x %*% 1:dir)), , drop=FALSE]
      dir_data[index, ] <- x
      break
    }
  }
}
#ディレクトリをベクトルに変換
dir_vec <- as.numeric(dir_data %*% 1:dir)


##パラメータの設定
#ディリクレ分布の事前分布
alpha11 <- rep(0.15, k1)
alpha12 <- rep(0.1, k2)
alpha21 <- c(rep(0.001, length(1:v1)), rep(0.1, length(1:v2)))
alpha22 <- c(rep(0.075, length(1:v1)), rep(0.001, length(1:v2)))

##すべての単語が出現するまでデータの生成を続ける
for(rp in 1:1000){
  print(rp)
  
  #ディリクレ分布からパラメータを生成
  theta1 <- thetat1 <- extraDistr::rdirichlet(d, alpha11)
  theta2 <- thetat2 <- extraDistr::rdirichlet(dir, alpha12)
  phi1 <- phit1 <- extraDistr::rdirichlet(k1, alpha21)
  phi2 <- phit2 <- extraDistr::rdirichlet(k2, alpha22)
  
  #スイッチング変数を生成
  gamma_list <- list()
  for(i in 1:d){
    if(dir_freq[i]==1){
      par <- c(6.25, 5.0)
      gamma_list[[i]] <- rbeta(1, par[1], par[2])
    } else {
      par <- c(5.0, runif(dir_freq[i], 1.0, 4.5))
      gamma_list[[i]] <- as.numeric(extraDistr::rdirichlet(1, par))
    }
  }

  ##モデルに基づきデータを生成
  word_list <- wd_list <- Z1 <- z1_list <- z21_list <- z22_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  
  for(i in 1:d){
    #スイッチング変数を生成
    n <- dir_freq[i] + 1
    if(dir_freq[i]==1){
      z1 <- rbinom(w[i], 1, gamma_list[[i]])
      Z1[[i]] <- cbind(z1, 1-z1)
      z1_list[[i]] <- as.numeric((Z1[[i]] * matrix(c(1, dir_vec[dir_index[[i]]]), nrow=w[i], ncol=n, byrow=T)) %*% rep(1, n))
    } else {
      Z1[[i]] <- rmnom(w[i], 1, gamma_list[[i]])
      z1_list[[i]] <- as.numeric((Z1[[i]] * matrix(c(1, dir_vec[dir_index[[i]]]), nrow=w[i], ncol=n, byrow=T)) %*% rep(1, n))
    }
    
    #多項分布より一般語のトピックを生成
    z21 <- matrix(0, nrow=w[i], ncol=k1)
    index1 <- which(Z1[[i]][, 1]==1)
    z21[index1, ] <- rmnom(length(index1), 1, theta1[i, ])
    z21_vec <- as.numeric(z21 %*% 1:k1)
    
    #多項分布よりディレクトリのトピックを生成
    z22 <- matrix(0, nrow=w[i], ncol=k2)
    index2 <- which(Z1[[i]][, 1]==0)
    z22[index2, ] <- rmnom(length(index2), 1, theta2[z1_list[[i]][index2], ])
    z22_vec <- as.numeric(z22 %*% 1:k2)
    
    #トピックおよびディレクトリから単語を生成
    word <- matrix(0, nrow=w[i], ncol=v)
    word[index1, ] <- rmnom(length(index1), 1, phi1[z21_vec[index1], ])   #トピックから単語を生成
    word[index2, ] <- rmnom(length(index2), 1, phi2[z22_vec[index2], ])   #ディレクトリから単語を生成
    wd <- as.numeric(word %*% 1:v)
    storage.mode(word) <- "integer"
    
    #データを格納
    z21_list[[i]] <- z21
    z22_list[[i]] <- z22
    wd_list[[i]] <- wd
    word_list[[i]] <- word
    WX[i, ] <- colSums(word)
  }
  #全単語が出現していたらbreak
  if(min(colSums(WX) > 0)) break
}

##リストを変換
wd <- unlist(wd_list)
Z21 <- do.call(rbind, z21_list)
Z22 <- do.call(rbind, z22_list)
z21_vec <- as.numeric(Z21 %*% 1:k1)
z22_vec <- as.numeric(Z22 %*% 1:k2)
sparse_data <- sparseMatrix(i=1:f, j=wd, x=rep(1, f), dims=c(f, v))
sparse_data_T <- t(sparse_data)
rm(word_list); rm(wd_list); rm(z21_list); rm(z22_list)
gc(); gc()


#####マルコフ連鎖モンテカルロ法でDLDAを推定####
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

##事前分布の設定
alpha1 <- 0.1
alpha2 <- 0.1
beta1 <- 1
beta2 <- 1

##初期値の設定
#トピック分布の初期値
theta1 <- extraDistr::rdirichlet(d, rep(1.0, k1))
theta2 <- extraDistr::rdirichlet(d, rep(1.0, k2))
phi1 <- extraDistr::rdirichlet(k1, rep(1.0, v))
phi2 <- extraDistr::rdirichlet(k2, rep(1.0, v))

#スイッチング分布の初期値
gamma <- list()
for(i in 1:d){
  if(dir_freq[i]==1){
    gamma[[i]] <- 0.5
  } else {
    n <- dir_freq[i]+1
    gamma[[i]] <- rep(1/n, n)
  }
}

##真値の設定
theta1 <- thetat1
theta2 <- thetat2
phi1 <- phit1
phi2 <- phit2
gamma <- gamma_list

##パラメータの保存用配列
THETA1 <- array(0, dim=c(d, k1, R/keep))
THETA2 <- array(0, dim=c(dir, k2, R/keep))
PHI1 <- array(0, dim=c(k1, v, R/keep))
PHI2 <- array(0, dim=c(k2, v, R/keep))
SEG11 <- rep(0, f)
SEG12 <- matrix(0, nrow=f, ncol=dir)
SEG21 <- matrix(0, nrow=f, ncol=k1)
SEG22 <- matrix(0, nrow=f, ncol=k2)
storage.mode(SEG12) <- "integer"
storage.mode(SEG21) <- "integer"
storage.mode(SEG22) <- "integer"

##インデックスを設定
#文書と単語のインデックスを作成
doc_list <- doc_vec <- list()
wd_list <- wd_vec <- list()
for(i in 1:d){
  doc_list[[i]] <- which(d_id==i)
  doc_vec[[i]] <- rep(1, length(doc_list[[i]]))
}
for(j in 1:v){
  wd_list[[j]] <- which(wd==j)
  wd_vec[[j]] <- rep(1, length(wd_list[[j]]))
}


##データの設定
#ディレクトリの割当を設定
dir_z <- matrix(0, nrow=d, ncol=dir)
dir_list1 <- dir_list2 <- list()
for(i in 1:d){
  dir_z[i, ] <- colSums(dir_data[dir_index[[i]], , drop=FALSE])
  dir_list1[[i]] <- (dir_z[i, ] * 1:dir)[dir_z[i, ] > 0]
  dir_list2[[i]] <- matrix(dir_list1[[i]], nrow=w[i], ncol=dir_freq[i], byrow=T)
}
dir_Z <- dir_z[d_id, ]
storage.mode(dir_Z) <- "integer"

#ディレクトリ数ごとにディレクトリを作成
max_freq <- max(dir_freq)
dir_no <- dir_Z * matrix(1:dir, nrow=f, ncol=dir, byrow=T)
freq_index1 <- freq_index2 <- list()
freq_word <- rep(0, max_freq)

for(j in 1:max_freq){
  x <- as.numeric(t(dir_Z * matrix(dir_freq[d_id], nrow=f, ncol=dir)))
  freq_index1[[j]] <- which(dir_freq[d_id]==j)
  freq_index2[[j]] <- which(x[x!=0]==j)
  freq_word[j] <- length(freq_index[[j]])/j
}
x <- as.numeric(t(dir_no)); dir_v <- x[x!=0]   #ディレクトリ数に合わせたディレクトリベクトル
x <- as.numeric(t(dir_Z * matrix(1:f, nrow=f, ncol=dir))); wd_v <- wd[x[x!=0]]   #ディレクトリ数に合わせた単語ベクトル
vec1 <- rep(1, k1); vec2 <- rep(1, k2)

##対数尤度の基準値
LLst <- sum(sparse_data %*% log(colMeans(sparse_data)))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##単語ごとのスイッチング変数を生成
  #トピックとディレクトリの期待尤度
  Lho1 <- theta1[d_id, ] * t(phi1)[wd, ]
  Lho2 <- theta2[dir_v, ] * t(phi2)[wd_v, ]
  
  Li1 <- as.numeric(Lho1 %*% vec1)   #トピックの期待尤度
  Li2 <- matrix(0, nrow=f, ncol=max_freq)   #ディレクトリの期待尤度
  for(j in 1:max_freq){
    Li2[freq_index1[[j]], 1:j] <- matrix(as.numeric(Lho2[freq_index2[[j]], ] %*% vec2), nrow=freq_word[[j]], ncol=j, byrow=T)
  }
  
  #ベルヌーイ分布あるいは多項分布よりスイッチング変数を生成
  Zi11 <- Zi12 <- rep(0, f)
  Lho_list <- list()
  
  for(i in 1:d){
    if(dir_freq[i]==1){
      
      #潜在変数zの設定
      omega <- matrix(c(gamma[[i]], 1-gamma[[i]]), nrow=w[i], ncol=dir_freq[i]+1, byrow=T)
      z_par <- omega * cbind(Li1[doc_list[[i]]], Li2[doc_list[[i]], 1:dir_freq[i]])
      z_rate <- z_par[, 1] / rowSums(z_par)
      Lho_list[[i]] <- z_par
      
      #ベルヌーイ分布よりスイッチング変数を生成
      z1 <- rbinom(w[i], 1, z_rate)
      Zi11[doc_list[[i]]] <- z1   #トピックに関係のある単語
      Zi12[doc_list[[i]]] <- (1-z1) * dir_list1[[i]]   #ディレクトリに関係のある単語
      
      #ベータ分布から混合率をサンプリング
      z_freq <- t(z1) %*% doc_vec[[i]]
      gamma[[i]] <- rbeta(1, z_freq+beta1, w[i]-z_freq+beta2)
      
    } else {
      
      #潜在変数zの設定
      omega <- matrix(gamma[[i]], nrow=w[i], ncol=dir_freq[i]+1, byrow=T)
      z_par <- omega * cbind(Li1[doc_list[[i]]], Li2[doc_list[[i]], 1:dir_freq[i]])
      z_rate <- z_par / rowSums(z_par)
      Lho_list[[i]] <- z_par
      
      z1 <- rmnom(w[i], 1, z_rate)   #スイッチング変数を生成
      Zi11[doc_list[[i]]] <- z1[, 1]   #トピックに関係のある単語
      Zi12[doc_list[[i]]] <- (z1[, -1] * dir_list2[[i]]) %*% rep(1, length(dir_list1[[i]]))   #ディレクトリに関係のある単語

      #ディリクレ分布から混合率をサンプリング「
      z_freq <- as.numeric(t(z1) %*% doc_vec[[i]])
      gamma[[i]] <- as.numeric(extraDistr::rdirichlet(1, as.numeric(t(z1) %*% doc_vec[[i]]) + alpha1))
    }
  }

  #生成したスイッチング変数のインデックスを作成
  index_z11 <- which(Zi11==1)
  
  
  ##多項分布からトピックをサンプリング
  #一般語トピックをサンプリング
  Zi21 <- matrix(0, nrow=f, ncol=k1)
  z_rate <- Lho1[index_z11, ] / Li1[index_z11]   #トピックの割当確率
  Zi21[index_z11, ] <- rmnom(length(index_z11), 1, z_rate)   #トピックをサンプリング
  Zi21_T <- t(Zi21)
  
  #ディレクトリトピックをサンプリング
  Lho2
  
  Zi12
  
  
  ##トピック分布のパラメータをサンプリング
  wsum0 <- matrix(0, nrow=d, ncol=k)
  for(i in 1:d){
    wsum0[i, ] <- Zi2_T[, doc_list[[i]], drop=FALSE] %*% doc_vec[[i]]
  }
  wsum <- wsum0 + alpha1   #ディリクレ分布のパラメータ
  theta <- extraDistr::rdirichlet(d, wsum)   #パラメータをサンプリング
  
  
  ##単語分布のパラメータをサンプリング
  #トピックの単語分布をサンプリング
  vsum0 <- matrix(0, nrow=k, ncol=v)
  for(j in 1:v){
    vsum0[, j] <- Zi2_T[, wd_list[[j]], drop=FALSE] %*% wd_vec[[j]]
  }
  vsum <- vsum0 + alpha2   
  phi1 <- extraDistr::rdirichlet(k, vsum)   #パラメータをサンプリング
  
  #ディレクトリの単語分布をサンプリング
  Zi0 <- Zi12[-index_z11]
  sparse_data0 <- sparse_data[-index_z11, ]
  dsum0 <- matrix(0, nrow=dir2, ncol=v)
  for(j in 1:dir2){
    dsum0[j, ] <- colSums(sparse_data0[Zi0==j, ]) 
  }
  dsum <- dsum0 + alpha2   #ディリクレ分布のパラメータ
  phi2 <- extraDistr::rdirichlet(dir2, dsum)   #パラメータをサンプリング
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    PHI1[, , mkeep] <- phi1
    PHI2[, , mkeep] <- phi2
    THETA[, , mkeep] <- theta
  }  
  
  #トピック割当はバーンイン期間を超えたら格納する
  if(rp%%keep==0 & rp >= burnin){
    SEG0 <- matrix(0, nrow=f, ncol=dir2)
    for(i in 1:f){
      if(Zi12[i]==0) next
      SEG0[i, Zi12[i]] <- 1
    }
    SEG11 <- SEG11 + Zi11
    SEG12 <- SEG12 + SEG0
    SEG2 <- SEG2 + Zi2
  }
  
  if(rp%%disp==0){
    #対数尤度を計算
    Lho <- rep(0, d)
    for(i in 1:d){
      Lho[i] <- sum(log(rowSums(Lho_list[[i]])))
    }
    #サンプリング結果を確認
    print(rp)
    print(c(sum(Lho), LLst))
    print(mean(Zi11))
    print(round(cbind(phi2[, (v11-4):(v11+5)], phit2[, (v11-4):(v11+5)]), 3))
  }
}





