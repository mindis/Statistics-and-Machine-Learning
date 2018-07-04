#####Tree_Structured_Switching_LDA_model#####
options(warn=0)
library(stringr)
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
#set.seed(2506787)

####データの発生####
##データの設定
s <- 3   #データタイプ数
k0 <- 10   #一般語のトピック数
k11 <- 3   #ユーザーの上位トピック数
k12 <- rtpois(k11, a=1, b=5, 2.5)   #ユーザーの下位トピックス
k21 <- 3   #アイテムの上位トピック数
k22 <- rtpois(k21, a=1, b=5, 2.5)   #アイテムの下位トピック数
hh <- 1000   #レビュアー数
item <- 500   #アイテム数
v1 <- 300   #評価スコアの語彙数
v2 <- 400   #ユーザートピックの語彙数
v3 <- 400   #アイテムトピックの語彙数
v <- v1 + v2 + v3   #総語彙数
spl <- matrix(1:v1, nrow=s, ncol=v1/s, byrow=T)
v1_index <- 1:v1
v2_index <- (v1+1):(v1+v2)
v3_index <- (v2+1):v

##IDと欠損ベクトルの作成
#IDを仮設定
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)

#欠損ベクトルを作成
for(rp in 1:100){
  m_vec <- rep(0, hh*item)
  for(i in 1:item){
    prob <- runif(1, 0.02, 0.12)
    m_vec[item_id0==i] <- rbinom(hh, 1, prob)
  }
  m_index <- which(m_vec==1)
  
  #完全なIDを設定
  user_id <- user_id0[m_index]
  item_id <- item_id0[m_index]
  d <- length(user_id)   #総レビュー数
  
  #すべてのパターンが生成されればbreak
  if(length(unique(user_id))==hh & length(unique(item_id))==item) break
}

#単語数を設定
w <- rpois(d, rgamma(d, 25, 0.5))   #文書あたりの単語数
f <- sum(w)   #総単語数
n_user <- plyr::count(user_id)$freq
n_item <- plyr::count(item_id)$freq

#単語IDを設定
u_id <- rep(user_id, w)
i_id <- rep(item_id, w)
d_id <- rep(1:d, w)


#インデックスを設定
user_index <- list()
item_index <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
}


##パラメータの事前分布を設定
#単語分布の事前分布
alpha11 <- c(rep(0.2, v1), rep(0.001, v2+v3))
alpha12 <- c(rep(0.001, v1), rep(0.1, v2), rep(0.001, v3))
alpha13 <- c(rep(0.001, v1+v2), rep(0.1, v3))

#スイッチング変数の事前分布
alpha2 <- rep(1.5, s)

#木構造トピック分布の事前分布
alpha31 <- c(18.5, 15.0)   #通過確率の事前分布
alpha32 <- rep(0.2, max(c(k12, k22)))   #ノード選択確率の事前分布

###すべての単語が出現するまでデータの生成を続ける
for(rp in 1:1000){
    
  ##パラメータを生成
  #ディリクレ分布から単語分布を生成
  phi0 <- phit0 <- extraDistr::rdirichlet(k0, alpha11)
  phi11 <- phit11 <- extraDistr::rdirichlet(k11, alpha12)
  phi21 <- phit21 <- extraDistr::rdirichlet(k21, alpha13)
  phi12 <- phi22 <- list()
  for(j in 1:k11){
    phi12[[j]] <- extraDistr::rdirichlet(k12[j], alpha12)
  }
  for(j in 1:k21){
    phi22[[j]] <- extraDistr::rdirichlet(k22[j], alpha13)
  }
  phit12 <- phi12; phit22 <- phi22
  
  #ディリクレ分布からスイッチング変数を生成
  lambda <- lambdat <- extraDistr::rdirichlet(hh, alpha2)
  
  ##木構造トピックを生成
  #ベータ分布から通過確率を生成
  gamma1 <- gammat1 <- matrix(rbeta(hh*k11, alpha31[1], alpha31[2]), nrow=hh, ncol=k11)   #ユーザーの通過率を生成
  gamma2 <- gammat2 <- matrix(rbeta(item*k21, alpha31[1], alpha31[2]), nrow=item, ncol=k21)   #アイテムの通過率を生成
  
  #ディリクレ分布からノード生成確率を生成
  theta0 <- thetat0 <- as.numeric(extraDistr::rdirichlet(1, rep(2.5, k0)))
  theta11 <- thetat11 <- extraDistr::rdirichlet(hh, alpha32[1:k11])
  theta21 <- thetat21 <- extraDistr::rdirichlet(item, alpha32[1:k21])
  theta12 <- theta22 <- list()
  for(j in 1:k11){
    theta12[[j]] <- extraDistr::rdirichlet(hh, alpha32[1:k12[j]])
  }
  for(j in 1:k21){
    theta22[[j]] <- extraDistr::rdirichlet(item, alpha32[1:k22[j]])
  }
  thetat12 <- theta12; thetat22 <- theta22
  
  
  ##モデルに基づきデータを生成
  #データの格納用配列
  WX <- matrix(0, nrow=d, ncol=v)
  wd_list <- y_list <- list()
  Z0_list <- Z1_list <- Z2_list <- G_list <- list() 
  options(warn=2)
  for(i in 1:d){
    if(i%%1000==0){
      print(i)
    }
    ##ユーザーとアイテムを抽出
    word <- matrix(0, nrow=w[i], ncol=v)
    u_index <- user_id[i]
    i_index <- item_id[i]
    
    ##多項分布からスイッチング変数を生成
    y <- rmnom(w[i], 1, lambda[u_index, ])
    y_vec <- as.numeric(y %*% 1:s)
    index_y1 <- which(y[, 1]==1); index_y2 <- which(y[, 2]==1); index_y3 <- which(y[, 3]==1)
    
    ##一般語のトピックと単語を生成
    if(length(index_y1) > 0){
      #トピックを生成
      z0 <- matrix(0, nrow=w[i], ncol=k0)
      z0[index_y1, ] <- rmnom(length(index_y1), 1, theta0)
      
      #単語を生成
      word[index_y1, ] <- rmnom(length(index_y1), 1, phi0)
    }
    
    ##ユーザートピックを生成
    #データの格納用配列
    z1 <- matrix(0, nrow=w[i], ncol=2)
    g1 <- rep(0, w[i])
    
    if(length(index_y2) > 0){
  
      #上位階層のトピックを生成
      z11 <- rmnom(length(index_y2), 1, theta11[u_index, ])
      z11_vec <- as.numeric(z11 %*% 1:k11)
      z1[index_y2, 1] <- z11_vec
      
      #通過変数を生成
      g1[index_y2] <- rbinom(length(index_y2), 1, gamma1[u_index, ][z1[, 1]])
      
      #下位階層のトピックを生成
      index_g1 <- which(g1==1)
      if(length(index_g1) > 0){
        for(j in 1:length(index_g1)){
          z12 <- as.numeric(rmnom(1, 1, theta12[[z1[index_g1, 1][j]]][u_index, ]))
          z1[index_g1[j], 2] <- as.numeric(z12 %*% 1:length(z12))
        }
      }
      #トピックから単語を生成
      for(j in 1:length(index_y2)){
        node <- z1[index_y2[j], ]
        if(sum(node > 0)==1){
          word[index_y2[j], ] <- rmnom(1, 1, phi11[node[1], ])
        } else {
          word[index_y2[j], ] <- rmnom(1, 1, phi12[[node[1]]][node[2], ])
        }
      }
    }
    
    ##アイテムトピックを生成
    #データの格納用配列
    z2 <- matrix(0, nrow=w[i], ncol=2)
    g2 <- rep(0, w[i])
    
    if(length(index_y3) > 0){
      #上位階層のトピックを生成
      z21 <- rmnom(length(index_y3), 1, theta21[i_index, ])
      z21_vec <- as.numeric(z21 %*% 1:k21)
      z2[index_y3, 1] <- z21_vec
      
      #通過変数を生成
      g2[index_y3] <- rbinom(length(index_y3), 1, gamma2[i_index, ][z2[, 1]])
      
      #下位階層のトピックを生成
      index_g2 <- which(g2==1)
      if(length(index_g2) > 0){
        for(j in 1:length(index_g2)){
          z22 <- as.numeric(rmnom(1, 1, theta22[[z2[index_g2, 1][j]]][i_index, ]))
          z2[index_g2[j], 2] <- as.numeric(z22 %*% 1:length(z22))
        }
      }
      #トピックから単語を生成
      for(j in 1:length(index_y3)){
        node <- z2[index_y3[j], ]
        if(sum(node > 0)==1){
          word[index_y3[j], ] <- rmnom(1, 1, phi21[node[1], ])
        } else {
          word[index_y3[j], ] <- rmnom(1, 1, phi22[[node[1]]][node[2], ])
        }
      }
    }
    
    ##データを格納
    WX[i, ] <- colSums(word)
    wd_list[[i]] <- as.numeric(word %*% 1:v)
    y_list[[i]] <- y
    Z0_list[[i]] <- as.numeric(z0 %*% 1:k0)
    Z1_list[[i]] <- z1
    Z2_list[[i]] <- z2
    G_list[[i]] <- cbind(g1, g2)
  }
  if(min(colSums(WX)) > 0) break   #break条件
}

##リストを変換
wd <- unlist(wd_list)
y <- do.call(rbind, y_list)
Z0 <- unlist(Z0_list)
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)
G <- do.call(rbind, G_list)
sparse_data <- sparseMatrix(i=1:f, wd, x=rep(1, f), dims=c(f, v))
sparse_data_T <- t(sparse_data)
rm(wd_list); rm(y_list); rm(Z0_list); rm(Z1_list); rm(Z2_list); rm(G_list)
gc(); gc()


####マルコフ連鎖モンテカルロ法でTree Structured Switching LDA modelを推定####
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
burnin <- 1000
disp <- 10

##事前分布の設定
alpha01 <- 0.1
alpha02 <- 0.1
beta01 <- 1.0
beta02 <- 0.5

##パラメータの真値
lambda <- lambdat
gamma1 <- gammat1
gamma2 <- gammat2
theta0 <- thetat0
theta11 <- thetat11
theta12 <- thetat12
theta21 <- thetat21
theta22 <- thetat22
phi0 <- phit0
phi11 <- phit11
phi12 <- phit12
phi21 <- phit21
phi22 <- phit22

##データのインデックスを設定
vec0 <- rep(1, k0); vec11 <- rep(1, k11); vec21 <- rep(1, k21)


####ギブスサンプリングでパラメータをサンプリング####
##期待対数尤度を設定
#一般語の期待対数尤度
Lho0 <- matrix(theta0, nrow=f, ncol=k0, byrow=T) * t(phi0)[wd, ]  
Lis0 <- as.numeric(Lho0 %*% vec0)

#ユーザー関連語の期待対数尤度
Lho11 <- theta11[u_id, ] * t(phi11)[wd, ]
Lis11 <- as.numeric(((1-gamma1[u_id, ]) * Lho11) %*% vec11)
Lho12 <- list(); Lis12 <- matrix(0, nrow=f, ncol=k11)
for(j in 1:k11){
  Lho12[[j]] <- theta12[[j]][u_id, ] * t(phi12[[j]])[wd, ]
  Lis12[, j] <- gamma1[u_id, j] * as.numeric(Lho12[[j]] %*% rep(1, k12[j]))
}

#アイテム関連語の期待対数尤度
Lho21 <- theta21[i_id, ] * t(phi21)[wd, ]
Lis21 <- as.numeric(((1-gamma2[i_id, ]) * Lho21) %*% vec21)
Lho22 <- list(); Lis22 <- matrix(0, nrow=f, ncol=k21)
for(j in 1:k21){
  Lho22[[j]] <- theta22[[j]][i_id, ] * t(phi22[[j]])[wd, ]
  Lis22[, j] <- gamma2[i_id, j] * as.numeric(Lho22[[j]] %*% rep(1, k22[j]))
}

##多項分布からスイッチング変数をサンプリング
#スイッチング変数の割当確率
Lis12_1 <- theta11[u_id, ] * Lis12
Lis22_1 <- theta21[i_id, ] * Lis22
Lis1 <- Lis11 + as.numeric(Lis12_1 %*% vec11)
Lis2 <- Lis21 + as.numeric(Lis22_1 %*% vec21)
Lis_par <- lambda[u_id, ] * cbind(Lis0, Lis1, Lis2)
y_rate <- Lis_par / as.numeric(Lis_par %*% rep(1, s))   #潜在変数の割当確率

#スイッチング変数をサンプリング
yi <- rmnom(f, 1, y_rate)   #多項分布からサンプリング
y_vec <- as.numeric(yi %*% 1:s)
index_user <- which(yi[, 2]==1)


##ユーザートピックの木構造を生成
##上位階層のノードを生成
#トピックノードの割当確率
par <- ((1-gamma1[u_id[index_user], ]) * Lho11[index_user, ]) + Lis12_1[index_user, ]
par_rate <- par / as.numeric(par %*% vec11)   #潜在変数の割当確率

#多項分布からトピックノードを生成
Zi11 <- matrix(0, nrow=f, ncol=k11)
Zi11[index_user, ] <- rmnom(length(index_user), 1, par_rate)
z11_vec <- as.numeric(Zi11 %*% 1:k11)
Zi11_T <- t(Zi11)

##通過変数をサンプリング
#通過確率を設定
g1_rate <- Lis12_1[index_user, ] / par
g1_rate[is.nan(g1_rate)] <- 0

#ベルヌーイ分布から通過変数をサンプリング
g1 <- rep(0, f)
g1[index_user] <- rbinom(length(index_user), 1, as.numeric((g1_rate * Zi11[index_user, ]) %*% vec11))

##下位階層のノードを生成
Zi12_T_list <- list()
z12_vec <- matrix(0, nrow=f, ncol=k11)

for(j in 1:k11){
  #トピックノードの割当確率
  index_node <- which(Zi11[, j]==1 & g1==1)
  par <- Lho12[[j]][index_node, ]
  par_rate <- par / as.numeric(par%*% rep(1, ncol(par)))   #潜在変数の割当確率
  
  #多項分布からトピックノードを生成
  Zi12 <- matrix(0, nrow=f, ncol=k12[j])
  Zi12[index_node, ] <- rmnom(length(index_node), 1, par_rate)  
  Zi12_T_list[[j]] <- t(Zi12)
  z12_vec[, j] <- as.numeric(Zi12 %*% 1:k12[j])
}


##アイテムトピックの木構造を生成
##上位階層のノードを生成
#トピックノードの割当確率
par <- ((1-gamma1[u_id[index_user], ]) * Lho11[index_user, ]) + Lis12_1[index_user, ]
par_rate <- par / as.numeric(par %*% vec11)   #潜在変数の割当確率

#多項分布からトピックノードを生成
Zi11 <- matrix(0, nrow=f, ncol=k11)
Zi11[index_user, ] <- rmnom(length(index_user), 1, par_rate)
z11_vec <- as.numeric(Zi11 %*% 1:k11)
Zi11_T <- t(Zi11)

##通過変数をサンプリング
#通過確率を設定
g1_rate <- Lis12_1[index_user, ] / par
g1_rate[is.nan(g1_rate)] <- 0

#ベルヌーイ分布から通過変数をサンプリング
g1 <- rep(0, f)
g1[index_user] <- rbinom(length(index_user), 1, as.numeric((g1_rate * Zi11[index_user, ]) %*% vec11))

##下位階層のノードを生成
Zi12_T_list <- list()
z12_vec <- matrix(0, nrow=f, ncol=k11)

for(j in 1:k11){
  #トピックノードの割当確率
  index_node <- which(Zi11[, j]==1 & g1==1)
  par <- Lho12[[j]][index_node, ]
  par_rate <- par / as.numeric(par%*% rep(1, ncol(par)))   #潜在変数の割当確率
  
  #多項分布からトピックノードを生成
  Zi12 <- matrix(0, nrow=f, ncol=k12[j])
  Zi12[index_node, ] <- rmnom(length(index_node), 1, par_rate)  
  Zi12_T_list[[j]] <- t(Zi12)
  z12_vec[, j] <- as.numeric(Zi12 %*% 1:k12[j])
}



