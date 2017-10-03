#####CJSモデルによる生存率推定#####
library(MASS)
library(nlme)
library(glmm)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(89467)

####データの発生####

#パラメータの設定
n.occasions <- 6   #観測期間
marked <- 5000   #参入数
phi0 <- 0.7   #生存率
p0 <- 0.4   #訪問率

##訪問履歴データを発生させる
n.occasions <- dim(PHI)[2] + 1
CH <- matrix(0, nrow=marked, ncol=n.occasions)
Z0 <- matrix(0, nrow=marked, ncol=n.occasions) 
CH[, 1] <- 1
Z0[, 1] <- 1
index_z <- 1:marked

##訪問履歴を逐次的に発生
#生存しているかどうかの決定
for(t in 2:n.occasions){
  Z0[index_z, t] <- rbinom(length(index_z), 1, phi0)
  index_z <- subset(1:nrow(Z0), Z0[, t]==1)
  CH[index_z, t] <- rbinom(length(index_z), 1, p0)
}


#発生させたデータの確認
cbind(Z0, CH)


####マルコフ連鎖モンテカルロ法で潜在生存率を推定####
##アルゴリズムの設定
R <- 5000
keep <- 2
sbeta <- 1.5

##事前分布の設定
a1 <- a2 <- 1
b1 <- b2 <- 1

##初期値の設定
phi <- 0.5
p <- 0.5

##サンプリング結果の保存用配列
Z <- array(0, dim=c(sum(marked), n.occasions, R/keep))
Z_rate <- array(0, dim=c(sum(marked), n.occasions, R/keep))
PHI1 <- array(0, R/keep)
P1 <- array(0, R/keep)

#生存の定義
z0 <- matrix(0, nrow=nrow(CH), ncol=ncol(CH))
z0[, 1] <- 1

for(j in 2:ncol(CH)){
  if(j==ncol(CH)){
    z0[, j] <- ifelse(CH[, j] > 0, 1, 0)
  } else {
    z0[, j] <- ifelse(rowSums(CH[, j:ncol(CH)]) > 0, 1, 0)
  }
}

####MCMCで生存率を推定####
for(rp in 1:R){
  z <- matrix(0, nrow=sum(marked), ncol=n.occasions)
  z[, 1] <- 1
  index_z <- subset(1:nrow(z), z[, 1]==1)
  sur <- c()
  obz <- c()
  
  for(j in 2:(n.occasions)){
    
    ##状態プロセスのパラメータを発生
    #最初の訪問時の潜在変数を定義
    index_surv <- subset(1:nrow(z0), z0[, j]==1)   #観測変数から生存の指示変数を設定
    
    #生存しているかどうかの決定
    if(j!=n.occasions){
      z[index_z, j] <- rbinom(length(index_z), 1, phi*(1-p))   #生存しているかどうかを発生
      z[index_surv, j] <- 1   #観測時点は生存
      obz_vec <- CH[z[, j]==1, j]
    } else {
      z[index_z, j] <- rbinom(length(index_z), 1, phi)
      obz_vec <- CH[z[, j]==1 | z0[, j]==1, j]
    }
    
    #パラメータをサンプリング
    sur_vec <- z[index_z, j]
    phi <- rbeta(1, a1 + sum(sur)+sum(sur_vec), b1 + length(sur)-sum(sur) + length(sur_vec)-sum(sur_vec))   
    p <- rbeta(1, a2 + sum(obz)+sum(obz_vec), b2 + length(obz)-sum(obz) + length(obz_vec)-sum(obz_vec)) 
  
    #パラメータを更新
    index_z <- subset(1:nrow(z), z[, j]==1)   #潜在変数zの生存の指示変数を設定
    sur <- c(sur, sur_vec)
    obz <- c(obz, obz_vec)
  }

  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    Z[, , mkeep] <- z
    PHI1[mkeep] <- phi
    P1[mkeep] <- p
    print(rp)
    print(phi)
  }
}


plot(1:length(P1), P1, type="l")
plot(1:length(PHI1), PHI1, type="l")

burnin <- 1000
mean(PHI1[burnin:(R/keep)])
mean(P1[burnin:(R/keep)])

mean(sur_vec)
