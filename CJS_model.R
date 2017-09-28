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
marked <- rep(300, n.occasions-1)   #1年ごとの参入数
phi0 <- rep(0.50, n.occasions-1)   #生存率
p0 <- rep(0.4, n.occasions-1)   #訪問率

#生存率と再捕獲率の行列の定義
PHI <- matrix(phi0, nrow=sum(marked), ncol=n.occasions-1)
P <- matrix(p0, nrow=sum(marked), ncol=n.occasions-1)

##訪問履歴データを発生させる
n.occasions <- dim(PHI)[2] + 1
CH <- matrix(0, nrow=sum(marked), ncol=n.occasions)

#初めて訪問した期間を定義するベクトル
mark.occ <- rep(1:length(marked), marked[1:length(marked)])

#訪問履歴を逐次的に発生
for(i in 1:sum(marked)){
  CH[i, mark.occ[i]] <- 1
  if(mark.occ[i] == n.occasions) next   #観測期間最終なら次の個体へ
  
  for(t in (mark.occ[i]+1):n.occasions){
    #生存しているかどうかを発生
    sur <- rbinom(1, 1, PHI[i, t-1])
    if(sur==0) break   #死んでいたらbreak
    
    #生存しているなら今期訪問したかどうかを発生
    rp <- rbinom(1, 1, P[i, t-1])
    if(rp==1) CH[i, t] <- 1
  }
}

#初めての訪問期を特定
f <- apply(CH, 1, function(x) min(which(x != 0)))
f0 <- 1:(n.occasions-1) 
index_f <- list()
for(i in 1:(n.occasions-1)){
  index_f[[i]] <- subset(1:length(f), f==i)
}


####マルコフ連鎖モンテカルロ法で潜在生存率を推定####
##アルゴリズムの設定
R <- 5000
keep <- 2
sbeta <- 1.5

##事前分布の設定
a1 <- a2 <- 1
b1 <- b2 <- 1

##初期値の設定
phi <- 0.3
p <- 0.5

##サンプリング結果の保存用配列
Z <- array(0, dim=c(sum(marked), n.occasions, R/keep))
Z_rate <- array(0, dim=c(sum(marked), n.occasions, R/keep))
PHI1 <- array(0, R/keep)
P1 <- array(0, R/keep)

#生存の定義
z0 <- matrix(0, nrow=nrow(CH), ncol=ncol(CH))
for(j in 1:ncol(CH)){
  if(j==ncol(CH)){
    z0[, j] <- ifelse(CH[, j] > 0, 1, 0)
  } else {
    z0[, j] <- ifelse(rowSums(CH[, j:ncol(CH)]) > 0, 1, 0)
  }
}


####MCMCで生存率を推定####
for(rp in 1:R){
  z <- matrix(0, nrow=sum(marked), ncol=n.occasions)
  sur <- c()
  obz <- c()
  
  for(j in 1:(n.occasions-1)){
    
    #最初の訪問時の潜在変数を定義
    z[index_f[[j]], j] <- 1
    
    for(t in (f0[j]+1):n.occasions){
    
      ##状態プロセス
      #前期の潜在変数zを条件付けた条件付き観測データを抽出
      index_z <- subset(index_f[[j]], z[index_f[[j]], t-1]==1)
      ch <- CH[index_z, t]
      
      #潜在変数zの混合率を計算
      if(t > f0[j]+1){
        Li1 <- dbinom(ch, 1, phi) * sur_rate
        Li2 <- dbinom(ch, 1, 1-phi) * (1-sur_rate)
      } else {
        Li1 <- dbinom(ch, 1, phi) 
        Li2 <- dbinom(ch, 1, 1-phi) 
      }
      
      z_rate <- Li1 / (Li1 + Li2)
      z_rate[z0[index_z, t]==1] <- 1   #観測されているなら生存率は1

      #生存有無zを発生
      z[index_z, t] <- rbinom(length(z_rate), 1, z_rate)  
      sur_rate <- mean(z[index_z, t])
      
      #サンプリングした潜在変数をもとに生存と観測をベクトル化
      sur <- c(sur, z[index_z, t])
      obz <- c(obz, CH[index_z, t][z[index_z, t]==1])
    }
  }
  
  ##パラメータを更新
  phi <- rbeta(1, a1 + sum(sur), b1 + length(sur)-sum(sur))   #ベータ分布から生存率をサンプリング
  p <- rbeta(1, a2 + sum(obz), b2 + length(obz)-sum(obz))   #ベター分布から観測率をサンプリング
  
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

length(sur)
obz

plot(1:length(P1), P1, type="l")
plot(1:length(PHI1), PHI1, type="l")

burnin <- 1000
mean(PHI1[burnin:(R/keep)])
mean(P1[burnin:(R/keep)])


