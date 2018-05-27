#####Hierarchical Probablistic Automaton topic model#####
options(warn=0)
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
k1 <- 6   #上位階層
k21 <- 10   #トピック数
k22 <- 4   #機能語のBIESオートマトン
k3 <- 4   #トピック語のBIESオートマトン
d <- 3000   #文書数
v1 <- 700   #一般語
v2 <- 400   #独立語
v3 <- 100   #機能語
v <- v1 + v2 + v3   #総語彙数
w <- rpois(d, rgamma(d, 70, 0.4))   #語彙数
f <- sum(w)   #総語彙数
B <- 1; I <- 2; E <- 3; S <- 4   #BIESのインデックス
topic_index <- matrix(1:v1, nrow=k21, ncol=v1/k21, byrow=T)
functional_index <- matrix((v1+1):(v1+v2), nrow=k1-1, ncol=v2/(k1-1), byrow=T)

#IDの設定
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}


##パラメータの事前分布の設定
#初期分布の事前分布
pi01 <- c(rep(30, 2), rep(5, k1-2))

#マルコフ推移行列の事前分布
alpha011 <- c(60.0, rep(15.0, k1-1))   #トピック割当ての条件付き最上位階層の事前分布
alpha012 <- c(10.0, rep(1.0, k1-1))   #機能語割当の最上位階層の事前分布
alpha021 <- rep(0.1, k21)   #トピック分布の事前分布
alpha022 <- c(15, 30)   #BIEとstop wordの事前分布
alpha023 <- c(10, 20)   #機能語のBIESオートマトンの事前分布
alpha031 <- c(15, 30)   #BIEとstop wordの事前分布
alpha032 <- c(10, 20)   #トピック語のBIESオートマトンの事前分布


##全単語が生成されるまでデータの生成を継続
for(rp in 1:1000){
  print(rp)
  
  ##パラメータの生成
  #マルコフ推移行列のパラメータを生成
  eta <- etat <- as.numeric(extraDistr::rdirichlet(1, pi01))
  theta1 <- thetat1 <- rbind(as.numeric(extraDistr::rdirichlet(1, alpha011)), extraDistr::rdirichlet(k1-1, alpha012))
  theta21 <- thetat21 <- extraDistr::rdirichlet(d, alpha021)
  beta21 <- betat21 <- rbeta(k1-1, alpha022[1], alpha022[2])
  beta22 <- betat22 <- rbeta(k1-1, alpha023[1], alpha023[2])
  beta31 <- betat31 <- rbeta(k21, alpha031[1], alpha031[2])
  beta32 <- betat32 <- rbeta(k21, alpha032[1], alpha032[2])
  
  ##トピックの単語分布を生成
  phi <- array(0, dim=c(k21, v, k22))
  
  for(j in 1:k21){
    #事前分布を設定
    alpha11 <- c(rep(0.15, v1), rep(0.001, v2), rep(0.001, v3))   #トピック語のBの事前分布
    alpha12 <- c(rep(0.25, v1), rep(0.001, v2), rep(0.2, v3))   #トピック語のIの事前分布
    alpha13 <- c(rep(0.15, v1), rep(0.001, v2), rep(0.01, v3))   #トピック語のEの事前分布
    alpha14 <- c(rep(0.15, v1), rep(0.001, v2), rep(0.001, v3))   #トピック語のSの事前分布
    
    #トピックに応じて事前分布を設定
    alpha11[as.numeric(t(topic_index[-j, ]))] <- 0.015
    alpha12[as.numeric(t(topic_index[-j, ]))] <- 0.025
    alpha13[as.numeric(t(topic_index[-j, ]))] <- 0.015
    alpha14[as.numeric(t(topic_index[-j, ]))] <- 0.015
    
    #パラメータを生成
    phi[j, , 1] <- extraDistr::rdirichlet(1, alpha11)
    phi[j, , 2] <- extraDistr::rdirichlet(1, alpha12)
    phi[j, , 3] <- extraDistr::rdirichlet(1, alpha13)
    phi[j, , 4] <- extraDistr::rdirichlet(1, alpha14)
  }
  phit <- phi
  
  ##機能語の単語分布を生成
  omega <- array(0, dim=c(k1-1, v, k3))
  
  for(j in 1:(k1-1)){
    #事前分布を設定
    alpha21 <- c(rep(0.001, v1), rep(1.0, v2), rep(0.001, v3))   #独立語のBの事前分布
    alpha22 <- c(rep(0.001, v1), rep(0.75, v2), rep(0.5, v3))   #独立語のIの事前分布
    alpha23 <- c(rep(0.001, v1), rep(1.0, v2), rep(0.1, v3))   #独立語のEの事前分布
    alpha24 <- c(rep(0.001, v1), rep(1.0, v2), rep(0.001, v3))   #独立語のSの事前分布
    
    #トピックに応じて事前分布を設定
    alpha21[as.numeric(t(functional_index[-j, ]))] <- 0.1
    alpha22[as.numeric(t(functional_index[-j, ]))] <- 0.15
    alpha23[as.numeric(t(functional_index[-j, ]))] <- 0.1
    alpha24[as.numeric(t(functional_index[-j, ]))] <- 0.1
   
    #パラメータを生成
    omega[j, , 1] <- extraDistr::rdirichlet(1, alpha21)
    omega[j, , 2] <- extraDistr::rdirichlet(1, alpha22)
    omega[j, , 3] <- extraDistr::rdirichlet(1, alpha23)
    omega[j, , 4] <- extraDistr::rdirichlet(1, alpha24)
  }
  omegat <- omega
  
  #区切り文字のパラメータを生成
  delta <- rbeta(d, 15, 60)
  
  
  ##モデルに基づきデータを生成
  z1_list <- z2_list <- z3_list <- topic_list <- r2_list <- r3_list <- delimiter_list <- wd_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  
  for(i in 1:d){
    if(i%%100==0){
      print(i)
    }
    #データの格納用配列
    z1 <- matrix(0, nrow=w[i], ncol=k1)
    z1_vec <- rep(0, w[i])
    delimiter <- rep(0, w[i])
    topic <- matrix(0, nrow=w[i], ncol=k21)
    topic_vec <- rep(0, w[i])
    r3 <- r2 <- rep(0, w[i])
    z2 <- matrix(0, nrow=w[i], ncol=k22)
    z3 <- matrix(0, nrow=w[i], ncol=k3)
    word <- matrix(0, nrow=w[i], ncol=v)
    
    #1単語ごとに潜在変数と単語を生成
    for(j in 1:w[i]){
      
      #最上位階層を生成
      if(j==1 | sum(delimiter[j-1])==1){
        z1[j, ] <- as.numeric(rmnom(1, 1, eta))
        z1_vec[j] <- which.max(z1[j, ])
      }
      if(j > 1 & sum(delimiter[j-1])==0){
        if(sum(z2[j-1, c(E, S)]) > 0 | sum(z3[j-1, c(E, S)]) > 0){
          z1[j, ] <- rmnom(1, 1, theta1[z1_vec[j-1], ])
          z1_vec[j] <- which.max(z1[j, ])
        }
        if(sum(z2[j-1, c(E, S)])==0 & sum(z3[j-1, c(E, S)])==0){
          z1[j, ] <- z1[j-1, ]
          z1_vec[j] <- which.max(z1[j-1, ])      
        }
      }
      
      #機能語のスイッチング変数とトピックを生成
      if(z1_vec[j] > 1){
        if(j==1 | sum(z2[j-1, c(E, S)]) > 0 | sum(z3[j-1, c(E, S)]) > 0){
          #機能語のスイッチング変数を生成
          r2[j] <- rbinom(1, 1, beta21[z1_vec[j]-1])
          z2[j, S] <- 1-r2[j]
        }
      }
      
      if(z1_vec[j]==1){
        if(j==1 | sum(z2[j-1, c(E, S)]) > 0 | sum(z3[j-1, c(E, S)]) > 0){
          #トピックを生成
          topic[j, ] <- rmnom(1, 1, theta21[i, ])
          topic_vec[j] <- which.max(topic[j, ])
            
          #トピックのスイッチング変数を生成
          r3[j] <- rbinom(1, 1, beta31[topic_vec[j]])
          z3[j, S] <- 1-r3[j]
        }
        if(j > 1 & sum(z2[j-1, c(E, S)])==0 & sum(z3[j-1, c(E, S)])==0){
          topic[j, ] <- topic[j-1, ]
          topic_vec[j] <- topic_vec[j-1]
        }
      }
      
      #機能語とトピック語のBIEを生成
      #初期値Bを生成
      if(z1_vec[j] > 0){
        if(z1_vec[j] > 1 & r2[j]==1){
          z2[j, B] <- 1
        }
        if(z1_vec[j]==1 & r3[j]==1){
          z3[j, B] <- 1
        }
      }
      
      #継続値Iを生成
      if(j > 1 & sum(z2[j-1, B])==1){
        z2[j, I] <- 1
        z1[j, ] <- z1[j-1, ]
        z1_vec[j] <- z1_vec[j-1]
      }
      
      if(j > 1 & sum(z3[j-1, B])==1){
        z3[j, I] <- 1
        z1[j, ] <- z1[j-1, ]
        z1_vec[j] <- z1_vec[j-1]
        topic[j, ] <- topic[j-1, ]
        topic_vec[j] <- topic_vec[j-1]
      }
    
      #BIEが終わるかどうかを生成
      if(sum(z2[j-1, I])==1){
        x <- rbinom(1, 1, beta22[z1_vec[j-1]-1])
        z2[j, I] <- x
        z2[j, E] <- 1-x
      }
      if(sum(z3[j-1, I])==1){
        x <- rbinom(1, 1, beta32[topic_vec[j-1]])
        z3[j, I] <- x
        z3[j, E] <- 1-x
      }
      
      #区切り文字を生成
      if(sum(z2[j, c(E, S)]) > 0 | sum(z3[j, c(E, S)]) > 0){
        delimiter[j] <- rbinom(1, 1, delta[i])
      }
      
      #単語を生成
      if(z1_vec[j] > 1){
        word[j, ] <- rmnom(1, 1, omega[z1_vec[j]-1, , which.max(z2[j, ])])
      }
      if(z1_vec[j]==1){
        word[j, ] <- rmnom(1, 1, phi[topic_vec[j], , which.max(z3[j, ])])
      }
    }
    
    #データを格納
    z1_list[[i]] <- z1
    z2_list[[i]] <- z2
    z3_list[[i]] <- z3
    topic_list[[i]] <- topic
    r2_list[[i]] <- r2
    r3_list[[i]] <- r3
    delimiter_list[[i]] <- delimiter
    wd_list[[i]] <- as.numeric(word %*% 1:v)
    WX[i, ] <- colSums(word)
  }
  if(min(colSums(WX)) > 0) break
}

#リストを変換
Z1 <- do.call(rbind, z1_list)
Z2 <- do.call(rbind, z2_list)
Z3 <- do.call(rbind, z3_list)
topic <- do.call(rbind, topic_list)
r2 <- unlist(r2_list)
r3 <- unlist(r3_list)
delimiter <- unlist(delimiter_list)
wd <- unlist(wd_list)

#スパース行列に変換
sparse_data <- sparseMatrix(1:f, wd, dims=c(f, v))
sparse_data_T <- t(sparse_data)


####マルコフ連鎖モンテカルロ法でHierarchical Probablistic Automaton topic modelを推定####
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
alpha2 <- 1/k21
alpha3 <- c(1, 1) 
alpha4 <- 0.01

#パラメータの真値
eta <- etat
theta1 <- thetat1
theta21 <- thetat21
beta21 <- betat21
beta22 <- betat22
beta31 <- betat31
beta32 <- betat32
phi <- phit
omega <- omegat


##BIESの単語ごとの尤度関数を設定
#トピック語のBIES尤度
Lit_B <- theta21_d * t(phi[, , B])[wd, ]; Lit_I <- theta21_d * t(phi[, , I])[wd, ]
Lit_E <- theta21_d * t(phi[, , E])[wd, ]; Lit_S <- theta21_d * t(phi[, , S])[wd, ]

#機能語後のBIES尤度
Lif_B <- t(omega[, , B])[wd, ]; Lif_I <- t(omega[, , I])[wd, ]
Lif_E <- t(omega[, , E])[wd, ]; Lif_S <- t(omega[, , S])[wd, ]


eta

vec_topic <- rep(1, k21)
vec_functional <- rep(1, k1-1)

Lit_B %*% vec_topic
Lit_S %*% vec_topic

##上位階層の期待尤度を推定
#トピック分布の階層的期待尤度
beta31_d <- matrix(beta31, nrow=f, ncol=k21, byrow=T)
Lit_switch1 <- beta31_d * Lit_B; Lit_switch2 <- (1-beta31_d) * Lit_S 
Lit_switch <- eta[1] * ((Lit_switch1 + Lit_switch2) %*% vec_topic)

Lit_I
theta


#機能語分布の階層的期待尤度
beta21_d <- matrix(beta21 , nrow=f, ncol=k1-1, byrow=T)
Lif_switch1 <- beta21_d * Lif_B; Lif_switch2 <- (1-beta21_d) * Lif_S
Lif_switch <- (matrix(eta[-1], nrow=f, ncol=k1-1, byrow=T) * (Lif_switch1 + Lif_switch2)) %*% vec_functional

Li_switch1 <- cbind(Lit_switch, Lif_switch)



cbind(round(Li_switch1[index_t11, ] / rowSums(Li_switch1[index_t11, ]), 3), Z1[index_t11, 1], rowSums(Z1[index_t11, -1]) > 0)











