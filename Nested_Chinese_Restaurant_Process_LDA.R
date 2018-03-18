#####Nested Chinese Restaurant Process LDA#####
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
L <- 3   #階層数
k1 <- kt1 <- 1   #レベル1の階層数
k2 <- kt2 <- 4   #レベル2の階層数
k3 <- kt3 <- rtpois(k2, 3, a=1, b=Inf)   #レベル3の階層数
k <- sum(c(k1, k2, k3))   #総トピック数
d <- 2000   #文書数
v <- 1000   #語彙数
w <- rpois(d, rgamma(d, 85, 0.5))   #単語数
f <- sum(w)   #単語数

#IDを設定
d_id <- rep(1:d, w)

##データの生成
for(rp in 1:1000){
  print(rp)
  
  #ディレクリ分布のパラメータを設定
  alpha1 <- alphat1 <- c(0.1, 0.2, 0.225)
  alpha2 <- alphat2 <- rep(10.0, k2)
  alpha3 <- alphat3 <- list()
  for(j in 1:k2){
    alpha3[[j]] <- alphat3[[j]] <- rep(3.0, k3[j])
  }
  beta <- betat <- rep(0.085, v)
  
  #ディレクリ分布からパラメータを生成
  theta1 <- thetat1 <- extraDistr::rdirichlet(d, alpha1)
  theta2 <- thetat2 <- as.numeric(extraDistr::rdirichlet(1, alpha2))
  theta3 <- thetat3 <- list()
  for(j in 1:k2){
    theta3[[j]] <- thetat3[[j]] <- as.numeric(extraDistr::rdirichlet(1, alpha3[[j]]))
  }
  phi1 <- phit1 <- as.numeric(extraDistr::rdirichlet(k1, beta))
  phi2 <- phit2 <- extraDistr::rdirichlet(k2, beta)
  phi3 <- phit3 <- list()
  for(j in 1:k2){
    phi3[[j]] <- phit3[[j]] <- extraDistr::rdirichlet(k3[j], beta)
  }
  phi <- phit <- rbind(phi1, phi2, do.call(rbind, phi3))
  
  ##生成過程に基づき単語を生成
  Z1 <- matrix(0, nrow=d, ncol=L)
  Z1[, 1] <- 1
  Z12_list <- list()
  Z13_list <- list()
  Z2_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  data_list <- list()
  word_list <- list()
  
  for(i in 1:d){
    #ノードを生成
    z12 <- rmnom(1, 1, theta2) 
    Z1[i, 2] <- as.numeric(z12 %*% 1:k2)
    z13 <- rmnom(1, 1, theta3[[Z1[i, 2]]])
    Z1[i, 3] <- as.numeric(z13 %*% 1:k3[Z1[i, 2]])
    
    #トピックのレベルを生成
    z2 <- rmnom(w[i], 1, theta1[i, ])
    z2_vec <- as.numeric(z2 %*% 1:L)
    
    #レベルごとに単語を生成
    index <- list()
    words <- matrix(0, nrow=w[i], ncol=v)
    for(j in 1:L){
      index[[j]] <- which(z2_vec==j)
    }
    words[index[[1]], ] <- rmnom(length(index[[1]]), 1, phi1)
    words[index[[2]], ] <- rmnom(length(index[[2]]), 1, phi2[Z1[i, 2], ])
    words[index[[3]], ] <- rmnom(length(index[[3]]), 1, phi3[[Z1[i, 2]]][Z1[i, 3], ])  
    
    
    #データを格納
    Z12_list[[i]] <- z12
    Z2_list[[i]] <- z2
    WX[i, ] <- colSums(words)
    data_list[[i]] <- words
    word_list[[i]] <- as.numeric(words %*% 1:v)
  }
  if(min(colSums(WX)) > 0) break
}

#リストを変換
Z12 <- do.call(rbind, Z12_list)
Z2 <- do.call(rbind, Z2_list)
word_vec <- unlist(word_list)
Data <- do.call(rbind, data_list)
storage.mode(Data) <- "integer"
storage.mode(Z2) <- "integer"
storage.mode(WX) <- "integer"
sparse_wx <- as(WX, "CsparseMatrix")
sparse_data <- as(Data, "CsparseMatrix")
rm(data_list); rm(Z2_list); rm(word_list)
gc(); gc()


####マルコフ連鎖モンテカルロ法でnCRP-LDAを推定####
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
alpha1 <- 1
alpha2 <- 1
alpha3 <- 0.1
beta1 <- 1   #CRPの事前分布
beta2 <- 1

##パラメータの真値
#トピック数の真値
k1 <- kt1; k2 <- kt2; k3 <- kt3

#トピックモデルのパラメータの真値
theta1 <- thetat1; theta2 <- thetat2; theta3 <- thetat3
phi1 <- phit1; phi2 <- phit2; phi3 <- phit3
Zi2 <- Level <- Z2; Zi12 <- Z12
Zi13 <- list()
Zi13_table <- table(1:d, Z1[, 3])
for(j in 1:k2){ 
  Zi13[[j]] <- (matrix(as.numeric(Zi13_table), nrow=d, ncol=ncol(Zi13_table)) * Zi12[, j])[, 1:k3[j]]
}


#トピック数の初期値
k1 <- 1; k2 <- 2; k3 <- rep(2, k2)

#パラメータの初期値
theta1 <- extraDistr::rdirichlet(d, rep(1, L))
phi1 <- extraDistr::rdirichlet(1, rep(1.0, v))
theta2 <- as.numeric(extraDistr::rdirichlet(1, rep(1, 2)))
phi2 <- extraDistr::rdirichlet(k2, rep(1.0, v))
theta3 <- phi3 <- list()
for(j in 1:k2){
  theta3[[j]] <- extraDistr::rdirichlet(1, rep(1.0, k3[j]))
  phi3[[j]] <- extraDistr::rdirichlet(k3[j], rep(1.0, v))
}

#トピック割当の初期値
Zi2 <- rmnom(f, 1, rep(1/L, L))
Zi2 <- Z2
Zi12 <- rmnom(d, 1, rep(1/k2, k2))
Zi13 <- list()
for(j in 1:k2){
  Zi13[[j]] <- rmnom(d, 1, rep(0.5, k3[j]))
}


##インデックスを作成
doc_vec <- doc_list <- list()
for(i in 1:d){
  doc_list[[i]] <- which(d_id==i)
  doc_vec[[i]] <- rep(1, length(doc_list[[i]]))
}


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##レベル2のパスをサンプリング
  #新しい潜在変数のパラメータを結合
  #par0 <- extraDistr::rdirichlet(1, colSums((sparse_data * Zi2[, 2])[doc_list[[sample(1:d, 1)]], ]) + alpha3)
  par0 <- phit2[sample(1:kt2, 1), ]
  phi2_new <- rbind(phi2, par0)
  
  #文書ごとに対数尤度を設定
  LLho2 <- as.matrix(t((sparse_data * Zi2[, 2]) %*% t(log(phi2_new))))
  LLi2_T <- matrix(0, ncol=d, nrow=nrow(phi2_new))
  for(i in 1:d){
    LLi2_T[, i] <- LLho2[, doc_list[[i]]] %*% doc_vec[[i]]
  }
  LLi2 <- t(LLi2_T)
  Li2 <- exp(LLi2 - rowMaxs(LLi2))   #尤度に変換
  
  #CRPの設定
  gamma0 <- cbind(matrix(colSums(Zi12), nrow=d, ncol=k2, byrow=T) - Zi12, beta1)
  gamma1 <- Li2 * gamma0/(d-1-beta)
  
  
  #潜在変数zの割当確率の推定とZのサンプリング
  z12_rate <- gamma1 / rowSums(gamma1)
  Zi12 <- rmnom(d, 1, z12_rate)
  
  #新しいトピックの生成がないなら新しいトピックは消滅
  if(min(colSums(Zi12)) <= 10){
    index <- which(colSums(Zi12) <= 10)
    Zi12 <- Zi12[, -index]
    z12_vec <- as.numeric(Zi12 %*% 1:(ncol(Zi12)))
    k2 <- ncol(Zi12)   #トピック数を更新
    
  } else {
    
    #新しいトピックを生成
    z12_vec <- as.numeric(Zi12 %*% 1:(k2+1))
    k2 <- sum(colSums(Zi12) > 0)   #トピック数を更新
  }
  
  ##レベル3のパスのサンプリング
  #レベル3の新しいパラメータをサンプリング
  if(k2 > length(phi3)){
    par <- colSums(sparse_data * Zi12[d_id, 3] * Zi2[, 3]) + alpha3
    phi3[[k2]] <- as.numeric(extraDistr::rdirichlet(1, par))
    Zi13[[k2]] <- matrix(1, nrow=d, ncol=1)
    k3[k2] <- 1
  }
  
  #レベル2のノードごとに新しいノードをサンプリング
  z13_vec <- list()
  Zi13_list <- Zi13
  Zi13 <- list()
  
  for(j in 1:k2){
    #新しい潜在変数のパラメータを結合
    #par0 <- extraDistr::rdirichlet(1, colSums((sparse_data * Zi12[d_id, j] * Zi2[, 3])[doc_list[[sample(1:d, 1)]], ]) + alpha3)
    sparse_data * Zi12[d_id, j]
    
    phi3_new <- rbind(phi3[[j]], par0)
    
    #文書ごとに対数尤度を設定
    LLho3 <- as.matrix(t((sparse_data * Zi12[d_id, j] * Zi2[, 3]) %*% t(log(phi3_new))))
    LLi3_T <- matrix(0, nrow=nrow(phi3_new), ncol=d)
    for(i in 1:d){
      LLi3_T[, i] <- LLho3[, doc_list[[i]]] %*% doc_vec[[i]]
    }
    LLi3 <- t(LLi3_T)
    Li3 <- exp(LLi3 - rowMaxs(LLi3))   #尤度に変換
    
    #CRPの設定
    gamma0 <- cbind(matrix(colSums(Zi13_list[[j]]), nrow=d, ncol=k3[j], byrow=T) - Zi13_list[[j]], beta2)
    gamma1 <- Li3 * gamma0/(d-1-beta)
    
    #潜在変数zの割当確率の推定とZのサンプリング
    z13_rate <- gamma1 / rowSums(gamma1)
    Zi13[[j]] <- rmnom(d, 1, z13_rate)
    
    #新しいトピックの生成がないなら新しいトピックは消滅
    if(min(colSums(Zi13[[j]])) <= 5){
      index <- which(colSums(Zi13[[j]]) <= 5)
      Zi13[[j]] <- Zi13[[j]][, -index, drop=FALSE]
      z13_vec[[j]] <- as.numeric(Zi13[[j]] %*% 1:ncol(Zi13[[j]]))
      k3[j] <- ncol(Zi13[[j]])   #トピック数を更新
      
    } else {
      
      #新しいトピックを生成
      z13_vec[[j]] <- as.numeric(Zi13[[j]] %*% 1:(k3[j]+1))
      k3[j] <- sum(colSums(Zi13[[j]]) > 0)   #トピック数を更新
    }
  }

  ##パラメータを更新
  #レベル1の単語分布をサンプリング
  vsum11 <- colSums(sparse_data * Zi2[, 1]) + alpha3
  phi1 <- as.numeric(extraDistr::rdirichlet(1, vsum11))
  
  #レベル2の単語分布をサンプリング
  vsum12 <- as.matrix(t(Zi12[d_id, ]) %*% (sparse_data * Zi2[, 2])) + alpha3
  phi2 <- extraDistr::rdirichlet(nrow(vsum12), vsum12)
  
  #レベル3の単語分布をサンプリング
  phi3 <- list()
  for(j in 1:k2){
    vsum13 <- as.matrix(t(Zi13[[j]][d_id, ]) %*% (sparse_data * Zi2[, 3]) + alpha3)
    phi3[[j]] <- extraDistr::rdirichlet(nrow(vsum13), vsum13)
  }
  
  
  ##トピック割当をサンプリング
  #トピック割当の尤度を設定
  Zi12_d <- Zi12[d_id, ]
  par1 <- phit1[word_vec]
  par2 <- rowSums(t(phi2)[word_vec, ] * Zi12_d)
  par0 <- c()
  for(j in 1:length(phi3)){
    par0 <- cbind(par0, rowSums(t(phi3[[j]])[word_vec, , drop=FALSE] * Zi13[[j]][d_id, ]))
  }
  par3 <- rowSums(par0 * Zi12_d)
  z_par <- theta1[d_id, ] * cbind(par1, par2, par3)   #トピック尤度
  
  #多項分布からトピック割当をサンプリング
  z_rate <- z_par / rowSums(z_par)   #レベル割当確率
  Zi2 <- rmnom(f, 1, z_rate)   #多項分布からレベル割当を生成
  Zi2 <- Z2
  Zi2_T <- t(Zi2)
  
  
  ##トピック分布を更新
  #ディレクリ分布のパラメータ
  wsum0 <- matrix(0, nrow=d, ncol=L)
  for(i in 1:d){
    wsum0[i, ] <- Zi2_T[, doc_list[[i]]] %*% doc_vec[[i]]
  }
  wsum <- wsum0 + alpha2 
  
  #ディレクリ分布からパラメータをサンプリング
  theta <- extraDistr::rdirichlet(d, wsum)
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(rp%%keep==0 & rp >= burnin){
    }
    
    if(rp%%disp==0){
      #サンプリング結果を確認
      print(rp)
      print(c(k1, k2, k3))
      print(c(kt1, kt2, kt3))
      print(colSums(Zi12))
      print(round(cbind(theta[1:10, ], thetat1[1:10, ]), 3))
      print(round(rbind(phi2[1:nrow(phi2), 1:20], phit2[, 1:20]), 3))
    }
  }
}

