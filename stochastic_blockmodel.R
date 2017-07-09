#####共クラスタリング(確率的ブロックモデル)#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(gtools)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

#set.seed(318)

####データの発生####
#データの設定
N <- 1000   #ユーザー数
K <- 150   #アイテム数
sg_n <- 4   #ユーザーのセグメント数
sg_k <- 3   #アイテムのセグメント数

##パラメータとセグメントを決定し観測データを発生させる
#ユーザーセグメントを発生
alpha.s1 <- rep(5, sg_n)
pi.s1 <- rdirichlet(1, alpha.s1)
z.s1 <- t(rmultinom(N, 1, pi.s1))
zno1 <- z.s1 %*% 1:sg_n
z1_cnt <- as.numeric(table(zno1))

#アイテムセグメントの発生
alpha.s2 <- rep(5, sg_k)
pi.s2 <- rdirichlet(1, alpha.s2)
z.s2 <- t(rmultinom(K, 1, pi.s2))
zno2 <- z.s2 %*% 1:sg_k
z2_cnt <- as.numeric(table(zno2))

#観測変数のパラメータの設定
#ユーザーセグメント×アイテムセグメントのベータ事前分布のパラメータを発生
theta.s <- rbeta(sg_n*sg_k, 1.0, 2.0)
round(theta.m <- matrix(theta.s, nrow=sg_n, ncol=sg_k), 3)

#ベルヌーイ分布から観測行列を発生させる
X <- matrix(0, nrow=N, ncol=K)
for(i in 1:N){
  print(i)
  for(j in 1:K){
    X[i, j] <- rbinom(1, 1, theta.m[zno1[i], zno2[j]])
  }
}


####周辺化ギブスサンプリングで共クラスタリングを推定####
##MCMCアルゴリズムの設定
R <- 10000
keep <- 2
sbeta <- 1.5

##ハイパーパラメータの設定
#ディクレリ分布のパラメータ
alpha1 <- rep(5, sg_n)
alpha2 <- rep(5, sg_k)
alpha <- 5

#ベータ分布のパラメータ
alpha_b <- 1.0
theta_b <- 1.5
beta <- 1


####MCMCアルゴリズムのデータの設定####
##十分統計量を更新するためのIDとデータフォーマットを準備
#ユーザーごとに観測データを集計する
X.cnt_u <- apply(X, 1, sum)
V.cnt_u <- apply(1-X, 1, sum)

#アイテムごとに観測データを集計する
X.cnt_i <- apply(X, 2, sum)
V.cnt_i <- apply(1-X, 2, sum)

#ユーザーIDの設定
ID_X1 <- rep(1:N, X.cnt_u)
ID_V1 <- rep(1:N, V.cnt_u)
ID1 <- rep(1:N, rep(K, N))

#アイテムIDの設定
ID_X2 <- rep(1:K, X.cnt_i)
ID_V2 <- rep(1:K, V.cnt_i)
ID2 <- rep(1:K, rep(N, K))


##MCMC推定のためのデータの準備
#観測データをベクトル化
X.vec <- as.numeric(t(X))
V.vec <- as.numeric(t(1-X))

#観測データをユーザー単位で多次元配列化
#ユーザー用の配列
XVu <- array(0, dim=c(sg_k, ncol(X), N))
XZu_M <- array(0, dim=c(sg_n, sg_k, N))
WVu <- array(0, dim=c(sg_k, ncol(X), N))
WZu_M <- array(0, dim=c(sg_n, sg_k, N))

#アイテム用の配列
XVi <- array(0, dim=c(sg_n, nrow(X), K))
XZi_M <- array(0, dim=c(sg_n, sg_k, K))
WVi <- array(0, dim=c(sg_n, nrow(X), K))
WZi_M <- array(0, dim=c(sg_n, sg_k, K))

for(i in 1:N){
  XVu[, , i] <- matrix(X[i, ], nrow=sg_k, ncol=ncol(X), byrow=T)
  WVu[, , i] <- matrix((1-X[i, ]), nrow=sg_k, ncol=ncol(X), byrow=T)
}

#観測データをアイテム単位で多次元配列化
for(i in 1:K){
  XVi[, , i] <- matrix(X[, i], nrow=sg_n, ncol=nrow(X), byrow=T)
  WVi[, , i] <- matrix((1-X[, i]), nrow=sg_n, ncol=nrow(X), byrow=T)
}

##サンプリング結果の保存用配列
Z1_M <- matrix(0, nrow=R/keep, ncol=nrow(X))
Z2_M <- matrix(0, nrow=R/keep, ncol=ncol(X))
P1_M <- array(0, dim=c(N, sg_n, R/keep))
P2_M <- array(0, dim=c(K, sg_k, R/keep))

####周辺化ギブスサンプリングで潜在変数zをサンプリング####
##所属セグメントの初期値を設定
for(t in 1:10000){
  print(t)
  
  #kmeansで初期クラスタを決定
  z1 <- matrix(0, nrow=N, ncol=sg_n)
  clust1 <- kmeans(X, sg_n)$cluster
  for(i in 1:N){
    z1[i, clust1[i]] <- 1
  }
  Z1 <- as.numeric(z1 %*% 1:sg_n)
  
  z2 <- t(rmultinom(K, 1, rep(1/sg_k, sg_k)))
  Z2 <- as.numeric(z2 %*% 1:sg_k)
  
  ##統計量の初期値を計算
  M1k <- as.numeric(table(Z1))   #ユーザーセグメントの統計量   
  M1l <- as.numeric(table(Z2))   #アイテムセグメントの統計量
  
  Nkl <- matrix(0, nrow=sg_n, ncol=sg_k)
  Mkl <- matrix(0, nrow=sg_n, ncol=sg_k)
  Xkl <- array(0, dim=c(N, K, sg_n*sg_k))
  Wkl <- array(0, dim=c(N, K, sg_n*sg_k))
  Zkl <- array(0, dim=c(N, K, sg_n*sg_k))
  
  for(i in 1:sg_n){
    for(j in 1:sg_k){
      r <- (i-1) + j
      Zkl[, , r] <- as.matrix(ifelse(Z1==i, 1, 0), nrow=N, ncol=1) %*% ifelse(as.numeric(Z2)==j, 1, 0)
      Xkl[, , r] <- X * Zkl[, , r]
      Wkl[, , r] <- (1-X) * Zkl[, , r]
      Mkl[i, j] <- sum(Xkl[, , r])
      Nkl[i, j] <- sum(Wkl[, , r])
    }
  }
  
  ##ここからMCMCサンプリングで供クラスタリングを推定
  for(rp in 1:R){
    
    ##ユーザーセグメントをサンプリング
    ##ユーザーの行に対するセグメントの割当を解除
    M1k_hat <- matrix(M1k, nrow=N, ncol=length(M1k), byrow=T) - z1
    
    ##ユーザーの行と列に対するセグメント割当を解除
    #ユーザーごとにMkl、Nklの全体の合計値をそれぞれ配列に格納しておく
    Mkl_u <- array(Mkl, dim=c(sg_n, sg_k, N))
    Nkl_u <- array(Nkl, dim=c(sg_n, sg_k, N))
    
    #列のセグメントを行列形式に変更
    Zkl_ind <- matrix(Z2, nrow=N, ncol=K, byrow=T)
    Zkl_v <- as.numeric(t(Zkl_ind))
    
    #観測データにセグメント行列を割り当てる
    Mkl_v <- Zkl_v * X.vec
    Nkl_v <- Zkl_v * V.vec
    
    #観測データが観測されたベクトルのみ取得
    Mkl_zero <- Mkl_v[Mkl_v > 0]
    Nkl_zero <- Nkl_v[Nkl_v > 0]
    
    #ユーザーごとにセグメント割当を集計
    Mkl_cnt <- as.matrix(t(table(Mkl_zero, ID_X1)))
    Nkl_cnt <- as.matrix(t(table(Nkl_zero, ID_V1)))
    
    #ユーザーごとにセグメント割当を解除したMkl、Nklを格納するための配列
    Mkl_hat <- array(Mkl, dim=c(sg_n, sg_k, N))
    Nkl_hat <- array(Nkl, dim=c(sg_n, sg_k, N))
  
    #ユーザーごとにMkl、Nklのセグメント割当を解除
    for(i in 1:N){
      Mkl_hat[Z1[i], , i] <- Mkl_u[Z1[i], , i] - Mkl_cnt[i, ]
      Nkl_hat[Z1[i], , i] <- Nkl_u[Z1[i], , i] - Nkl_cnt[i, ]
    
      ##ユーザーごとにセグメント割当確率を計算
      #Z1の事後分布のサンプリング式を計算
      #ユーザーごとにセグメント割当を計算
      XZ <- XVu[, , i] * matrix(rep(Z2, sg_k), nrow=sg_k, ncol=ncol(X), byrow=T)
      XZ_sum <- rowSums(t(apply(cbind(XZ, 1:sg_k), 1, function(x) ifelse(x[-length(x)]==x[length(x)], 1, 0))))
      XZu_M[, , i] <- matrix(XZ_sum, nrow=sg_n, ncol=sg_k, byrow=T)
      
      WZ <- WVu[, , i] * matrix(rep(Z2, sg_k), nrow=sg_k, ncol=ncol(X), byrow=T)
      WZ_sum <- rowSums(t(apply(cbind(WZ, 1:sg_k), 1, function(x) ifelse(x[-length(x)]==x[length(x)], 1, 0))))
      WZu_M[, , i] <- matrix(WZ_sum, nrow=sg_n, ncol=sg_k, byrow=T)
    }
    
    #ベータ分布のパラメータ
    alpha_hat <- alpha + Mkl_hat
    beta_hat <- beta + Nkl_hat
    
    
    #Z1のサンプリング式の周辺事後分布
    sgu_M <- array(matrix(as.numeric(table(Z2)), nrow=sg_n, ncol=sg_k, byrow=T), dim=c(sg_n, sg_k, N))   #Z2の割当数
    
    Z1_gamma1 <- lgamma(alpha_hat + beta_hat) - lgamma(alpha_hat) - lgamma(beta_hat)
    Z1_gamma2 <- lgamma(alpha_hat + XZu_M) + lgamma(beta_hat + WZu_M)
    Z1_gamma3 <- lgamma(alpha_hat + beta_hat + sgu_M)
    Z1_gamma <- log(alpha_hat) + Z1_gamma1 + Z1_gamma2 - Z1_gamma3
    
    #数値が桁落ちしないように定数を加える
    M1 <- apply(Z1_gamma, 3, mean)
    Z1_gamma <- exp(Z1_gamma - array(M1[rep(1:N, rep(sg_n*sg_k, N))], dim=c(sg_n, sg_k, N)))
    
  
    #セグメント割当確率を計算して潜在変数を発生
    #ユーザーごとに割当確率を計算
    Pr_z1 <- t(apply(Z1_gamma, c(1, 3), prod))
    Pr1 <- Pr_z1 / rowSums(Pr_z1)
    
    #カテゴリカル分布から潜在変数を発生
    z1_new <- try(t(apply(Pr1, 1, function(x) rmultinom(1, 1, x))), silent=TRUE)
    if(class(z1_new) == "try-error") {break}   #エラー処理
    
    Z1_new <- z1_new %*% 1:sg_n
    Z1 <- Z1_new
    
    ##十分統計量を更新
    M1k <- as.numeric(table(Z1_new))
    
    for(i in 1:sg_n){
      for(j in 1:sg_k){
        r <- (i-1) + j
        Zkl[, , r] <- as.matrix(ifelse(Z1_new==i, 1, 0), nrow=N, ncol=1) %*% ifelse(as.numeric(Z2)==j, 1, 0)
        Xkl[, , r] <- X * Zkl[, , r]
        Wkl[, , r] <- (1-X) * Zkl[, , r]
        Mkl[i, j] <- sum(Xkl[, , r])
        Nkl[i, j] <- sum(Wkl[, , r])
      }
    }
    
    ##アイテムセグメントをサンプリング
    ##ユーザーの行に対するセグメントの割当を解除
    M1l_hat <- matrix(M1l, nrow=K, ncol=sg_k, byrow=T) - z2
    
    ##アイテムの行と列に対するセグメント割当を解除
    #アイテムごとにMkl、Nklの全体の合計値をそれぞれ配列に格納しておく
    Mkl_i <- array(Mkl, dim=c(sg_n, sg_k, K))
    Nkl_i <- array(Nkl, dim=c(sg_n, sg_k, K))
    
    #行のセグメントを行列形式に変更
    Zkl_ind <- matrix(Z1, nrow=N, ncol=K)
    Zkl_v <- as.numeric(Zkl_ind)
    
    #観測データにセグメント行列を割り当てる
    Mkl_v <- Zkl_v * X.vec
    Nkl_v <- Zkl_v * V.vec
    
    #観測データが観測されたベクトルのみ取得
    Mkl_zero <- Mkl_v[Mkl_v > 0]
    Nkl_zero <- Nkl_v[Nkl_v > 0]
    
    #アイテムごとにセグメント割当を集計
    Mkl_cnt <- as.matrix(t(table(Mkl_zero, ID_X2)))
    Nkl_cnt <- as.matrix(t(table(Nkl_zero, ID_V2)))
    
    #アイテムごとにセグメント割当を解除したMkl、Nklを格納するための配列
    Mkl_hat <- array(Mkl, dim=c(sg_n, sg_k, K))
    Nkl_hat <- array(Nkl, dim=c(sg_n, sg_k, K))
    
    #アイテムごとにMkl、Nklのセグメント割当を解除
    for(i in 1:K){
      Mkl_hat[, Z2[i], i] <- Mkl_i[, Z2[i], i] - Mkl_cnt[i, ]
      Nkl_hat[, Z2[i], i] <- Nkl_i[, Z2[i], i] - Nkl_cnt[i, ]
    
      ##アイテムごとにセグメント割当確率を計算
      #Z1の事後分布のサンプリング式を計算
      #ユーザーごとにセグメント割当を計算
      XZ <- XVi[, , i] * matrix(rep(Z1, sg_n), nrow=sg_n, ncol=nrow(X), byrow=T)
      XZ_sum <- rowSums(t(apply(cbind(XZ, 1:sg_n), 1, function(x) ifelse(x[-length(x)]==x[length(x)], 1, 0))))
      XZi_M[, , i] <- matrix(XZ_sum, nrow=sg_n, ncol=sg_k)
      
      WZ <- WVi[, , i] * matrix(rep(Z1, sg_n), nrow=sg_n, ncol=nrow(X), byrow=T)
      WZ_sum <- rowSums(t(apply(cbind(WZ, 1:sg_n), 1, function(x) ifelse(x[-length(x)]==x[length(x)], 1, 0))))
      WZi_M[, , i] <- matrix(WZ_sum, nrow=sg_n, ncol=sg_k)
    }
    
    #ベータ分布のパラメータ
    alpha_hat <- alpha + Mkl_hat
    beta_hat <- beta + Nkl_hat
        
    #Z2のサンプリング式の周辺事後分布
    sgi_M <- array(matrix(as.numeric(table(Z1)), nrow=sg_n, ncol=sg_k), dim=c(sg_n, sg_k, K))   #Z1の割当数
    
    Z2_gamma1 <- lgamma(alpha_hat + beta_hat) - lgamma(alpha_hat) - lgamma(beta_hat)
    Z2_gamma2 <- lgamma(alpha_hat + XZi_M) + lgamma(beta_hat + WZi_M)
    Z2_gamma3 <- lgamma(alpha_hat + beta_hat + sgi_M)
    Z2_gamma <- log(alpha_hat) + Z2_gamma1 + Z2_gamma2 - Z2_gamma3
    
    #数値が桁落ちしないように定数を加える
    M2 <- apply(Z2_gamma, 3, mean)
    Z2_gamma <- exp(Z2_gamma - array(M2[rep(1:K, rep(sg_n*sg_k, K))], dim=c(sg_n, sg_k, K)) )
  
    
    #セグメント割当確率を計算して潜在変数を発生
    #ユーザーごとに割当確率を計算
    Pr_z2 <- t(apply(Z2_gamma, 2:3, prod))
    Pr2 <- Pr_z2 / rowSums(Pr_z2)
  
    #カテゴリカル分布から潜在変数を発生
    z2_new <- try(t(apply(Pr2, 1, function(x) rmultinom(1, 1, x))), silent=TRUE)
    if(class(z2_new) == "try-error") {break}   #エラー処理
    
    Z2_new <- z2_new %*% 1:sg_k
    Z2 <- Z2_new
  
    ##十分統計量を更新
    M1l <- as.numeric(table(Z2_new))
    
    for(i in 1:sg_n){
      for(j in 1:sg_k){
        r <- (i-1) + j
        Zkl[, , r] <- as.matrix(ifelse(Z1==i, 1, 0), nrow=N, ncol=1) %*% ifelse(as.numeric(Z2_new)==j, 1, 0)
        Xkl[, , r] <- X * Zkl[, , r]
        Wkl[, , r] <- (1-X) * Zkl[, , r]
        Mkl[i, j] <- sum(Xkl[, , r])
        Nkl[i, j] <- sum(Wkl[, , r])
      }
    }
    
    ##パラメータを格納
    if(rp%%keep==0){
      print(rp)
      mkeep <- rp/keep
      Z1_M[mkeep, ] <- Z1
      Z2_M[mkeep, ] <- Z2
      P1_M[, , mkeep] <- Pr1
      P2_M[, , mkeep] <- Pr2
    
      print(rbind(M1l, z2_cnt))
      print(rbind(M1k, z1_cnt))
    }
  }
  if(R==rp) {break}
}
