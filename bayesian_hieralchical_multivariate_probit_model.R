#####bayesian hieralchical multivariate probit model#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)


####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #固有値分解で強制的に正定値行列に修正する
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

####データの発生####
##データの設定
s <- 2   #応答変数数
dir <- 150   #ディレクトリ数
item <- 10000   #アイテム数
dir_freq <- rtpois(item, 1.0, a=0, b=5)   #アイテムごとのディレクトリ数
max_dir <- max(dir_freq)   #ディレクトリの最大数
w <- rpois(item, rgamma(item, 10.0, 0.075))   #アイテムあたりのサンプル数
f <- sum(w)   #総サンプル数

#IDの設定
item_id <- rep(1:item, w)
t_id <- as.numeric(unlist(tapply(1:f, item_id, rank)))

#ディレクトリの生成
pr <- runif(dir, 0.1, 3.0)
dir_data <- as.numeric(rmnom(item, 1, pr) %*% 1:dir)
dir_item <- dir_data[item_id]

##説明変数の生成
#アイテムの説明変数
k1 <- 3; k2 <- 4; k3 <- 5
x1 <- matrix(runif(f*k1, 0, 1), nrow=f, ncol=k1)
x2 <- matrix(0, nrow=f, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(f, 1, pr)
}
x3 <- rmnom(f, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
X <- cbind(1, x1, x2, x3)   #データを結合
k <- ncol(X)

##パラメータを生成
#応答変数の誤差相関
Cov <- Covt <- matrix(c(1, 0.5, 0.5, 1.0), nrow=s, ncol=s)  

##応答変数が妥当な数値になるまで繰り返す
for(rp in 1:1000){
  
  #回帰パラメータを生成
  theta1 <- thetat1 <- runif(k, -1.3, 1.1); theta2 <- thetat2 <- runif(k, -1.4, 0.9)   #階層モデルの平均
  tau1 <-taut1 <- runif(k, 0.1, 0.7) * diag(k); tau2 <- taut2 <- runif(k, 0.1, 0.7) * diag(k)   #階層モデルの分散
  beta1 <- betat1 <- mvrnorm(dir, theta1, tau1); beta2 <- betat2 <- mvrnorm(dir, theta2, tau2)   #ディレクトリ別の回帰係数
  
  ##応答変数を生成
  #回帰モデルの平均構造
  
  mu1 <- as.numeric((X * beta1[dir_item, ]) %*% rep(1, k))
  mu2 <- as.numeric((X * beta2[dir_item, ]) %*% rep(1, k))
  mu <- cbind(mu1, mu2)
  
  #多変量正規分布から応答変数を生成
  er <- mvrnorm(f, rep(0, s), Cov)  
  UT <- U <- mu + er   #潜在効用を設定
  y <- ifelse(U > 0, 1, 0)   #応答絵変数を生成
  
  #ストップ判定
  if(max(colMeans(y)) < 0.4 & min(colMeans(y)) > 0.15) break
}

####マルコフ連鎖モンテカルロ法でLatent variable bayesian hieralchical probit modelを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##多変量正規分布の条件付き期待値と条件付き分散を計算する関数
cdMVN <- function(mean, Cov, dependent, U){
  
  #分散共分散行列のブロック行列を定義
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent, drop=FALSE]
  Cov21 <- Cov[-dependent, dependent, drop=FALSE]
  Cov22 <- Cov[-dependent, -dependent]
  
  
  #条件付き分散と条件付き平均を計算
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mean[, dependent] + t(CDinv %*% t(U[, -dependent] - mean[, -dependent]))   #条件付き平均を計算
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #条件付き分散を計算
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##多変量正規分布の密度関数
mvdnorm <- function(u, mu, Cov, s){
  er <- u - mu   #誤差
  Lho <- 1 / (sqrt(2*pi)^s*sqrt(det(Cov))) * exp(-1/2 * as.numeric((er %*% solve(Cov) * er) %*% rep(1, s)))
  return(Lho)
}

##MCMCの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##事前分布の設定
#モデルの事前分布
alpha0 <- 1.0   #ディレクトリ割当の事前分布
nu1 <- s   #逆ウィシャート分布の自由度
V1 <- nu1 * diag(s)   #逆ウィシャート分布の自由度

#階層モデルの事前分布
nu2 <- k   #逆ウィシャート分布の自由度
V2 <- nu2 * diag(k)   #逆ウィシャート分布の自由度
Deltabar <- rep(0, k)   #回帰係数の平均の事前分布
Adelta <- solve(diag(100, k))   #回帰係数の分散の事前分布

##パラメータの真値
theta1 <- thetat1
theta2 <- thetat2
tau1 <- taut1
tau2 <- taut2

diag(c(diag(taut1), diag(taut2)))

beta1 <- betat1
beta2 <- betat2
Cov <- Covt

##初期値を設定
lambda <- matrix(0, nrow=item, ncol=max_dir)
for(i in 1:item){
  if(dir_freq[i]==1){
    lambda[i, 1] <- 1
  } else {
    lambda[i, 1:dir_freq[i]] <- as.numeric(extraDistr::rdirichlet(1, rep(10.0, dir_freq[i])))
  }
}
theta1 <- theta2 <- rep(0, k)
tau1 <- tau2 <- diag(0.5, k) 
beta1 <- mvrnorm(dir, theta1, tau1)
beta2 <- mvrnorm(dir, theta2, tau2)
Cov <- matrix(c(1, 0.3, 0.3, 1), nrow=s, ncol=s)

##インデックスを作成
#ディレクトリのインデックス
dir_index <- list(); dir_list1 <- dir_list2 <- list()
for(i in 1:dir){
  data <- matrix(0, nrow=f, ncol=max_dir)
  for(j in 1:max_dir){
    index <- which(dir_item[, j]==i)
    if(length(index)==0) next
    data[index, j] <- index
  }
  index <- which(rowSums(data > 0) > 0)
  dir_list1[[i]] <- data[index, ]
  dir_list2[[i]] <- rowSums(dir_list1[[i]])
}

for(j in 1:max_dir){
  dir_index[[j]] <- which(dir_data[item_id, j] > 0)
}
multi_index <- which(rowSums(dir_data[item_id, ] > 0) >= 2)

#アイテムのインデックス
item_index <- item_vec <- list()
for(i in 1:item){
  item_index[[i]] <- which(item_id==i)
  item_vec[[i]] <- rep(1, length(item_index[[i]]))
}


##切断正規分布の切断領域を定義
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)
a0 <- ifelse(y==0, -10^-100, 0)
b0 <- ifelse(y==1, 10^-100, 0)


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  ##切断正規分布から潜在効用を生成
  U1 <- U2 <- Lho <- matrix(0, nrow=f, ncol=max_dir)
  
  #効用の平均構造
  for(i in 1:max_dir){
    index <- dir_index[[i]]
    target <- dir_item[, i][index]
    mu1 <- as.numeric((X[index, ] * beta1[target, ]) %*% rep(1, k))
    mu2 <- as.numeric((X[index, ] * beta2[target, ]) %*% rep(1, k))
    mu <- cbind(mu1, mu2)
    
    #潜在効用を生成
    for(j in 1:s){
      MVR <- cdMVN(mu, Cov, j, U[index, ])   #条件付き分布と分散を生成
      MVR.U <- MVR$CDmu
      MVR.S <- sqrt(MVR$CDvar)
      if(j==1){
        U1[index, i] <- rtnorm(MVR.U, MVR.S, a[index, j], b[index, j])   #パラメータを生成
        index_inf <- which(is.infinite(U1[index, i])==TRUE)
        U1[index, i][index_inf] <- 0
      } else {
        U2[index, i] <- rtnorm(MVR.U, MVR.S, a[index, j], b[index, j])
        index_inf <- which(is.infinite(U2[index, i])==TRUE)
        U2[index, i][index_inf] <- 0
      }
    }
    #多変量正規分布の尤度関数
    u <- cbind(U1[index, i], U2[index, i])
    Lho[index, i] <- mvdnorm(u, mu, Cov, s)
  }
  
  ##潜在ディレクトリを生成
  #多項分布から潜在変数をサンプリング
  Zi <- matrix(0, nrow=f, ncol=max_dir); Zi[-multi_index, 1] <- 1
  dir_prob <- Lho[multi_index, ] / as.numeric(Lho[multi_index, ] %*% rep(1, max_dir))   #ディレクトリの割当確率
  Zi[multi_index, ] <- rmnom(length(multi_index), 1, dir_prob)   #多項分布からサンプリング
  Zi_T <- t(Zi)
  
  #ディリクレ分布から潜在変数をサンプリング
  for(i in 1:item){
    if(dir_freq[i]==1) next
    wsum <- (Zi_T[, item_index[[i]]] %*% item_vec[[i]])[1:dir_freq[i]] + alpha0
    lambda[i, 1:dir_freq[i]] <- extraDistr::rdirichlet(1, wsum)
  }
  
  ##ディレクトリ別に回帰係数をサンプリング
  util_mu <- U <- matrix(0, nrow=f, ncol=s)
  for(i in 1:dir){
    #データを抽出
    for(i in 1:dir){
      print(i)
      index <- as.numeric((dir_list1[[i]] * Zi[dir_list2[[i]], ]) %*% rep(1, max_dir)); index <- index[index > 0]
      u <- cbind((U1[index, ]*Zi[index, ]) %*% rep(1, max_dir), (U2[index, ]*Zi[index, ]) %*% rep(1, max_dir))
      
      #回帰係数のパラメータ
      XVX <- matrix(0, nrow=s*k, ncol=s*k)
      XVY <- matrix(0, nrow=s*k, ncol=1)
      for(r in 1:length(index)){
        XV <- DT_T[, , index[r]] %*% inv_cov
        XVX <- XVX + XV %*% DT[, , index[r]]
        XVY <- XVY + XV %*% u[r, ]
      }
    }
    pmv
    
    inv_tau1 <- solve(tau1); inv_tau2 <- solve(tau2)
    inv_xxv1 <- solve(xx + inv_tau1); inv_xxv2 <- solve(xx + inv_tau2)
    beta_mu1 <- inv_xxv1 %*% (xu[, 1] + inv_tau1 %*% theta1)   #回帰パラメータの期待値
    beta_mu2 <- inv_xxv2 %*% (xu[, 2] + inv_tau2 %*% theta2)
    
    #多変量正規分布から回帰係数をサンプリング
    beta1[i, ] <- mvrnorm(1, beta_mu1, inv_xxv1)
    beta2[i, ] <- mvrnorm(1, beta_mu2, inv_xxv2)
    
    #効用の期待値
    beta <- cbind(beta1[i, ], beta2[i, ])
    util_mu[index, ] <- x %*% beta
    U[index, ] <- u
  }
  U1 <- U; util_mu1 <- util_mu
  
  ##モデルの相関行列をサンプリング
  #逆ウィシャート分布のパラメータ
  R_error <- U - util_mu
  IW_R <- solve(V1) + t(R_error) %*% R_error   #パラメータ
  Sn <-  nu1 + f
  
  #逆ウィシャート分布からパラメータをサンプリング
  Cov_hat <- rwishart(Sn, solve(IW_R))$IW
  
  
  ##階層モデルのパラメータをサンプリング
  #回帰係数の階層モデルをサンプリング
  mu1 <- colMeans(beta1); mu2 <- colMeans(beta2)
  theta_par1 <- solve(Adelta + dir*solve(tau1)) %*% (dir*solve(tau1) %*% mu1)
  theta_par2 <- solve(Adelta + dir*solve(tau2)) %*% (dir*solve(tau2) %*% mu2)
  theta1 <- mvrnorm(1, theta_par1, tau1/dir)
  theta2 <- mvrnorm(1, theta_par2, tau2/dir)
  
  #階層モデルの標準偏差の事後分布をサンプリング
  er1 <- beta1 - matrix(theta1, nrow=dir, ncol=k, byrow=T); er2 <- beta2 - matrix(theta2, nrow=dir, ncol=k, byrow=T)
  IW_R1 <- solve(V2) + t(er1) %*% er1; IW_R2 <- solve(V2) + t(er2) %*% er2
  Sn <- nu2 + dir
  tau1 <- rwishart(Sn, solve(IW_R1))$IW   #逆ウィシャート分布からtauをサンプリング
  tau2 <- rwishart(Sn, solve(IW_R2))$IW
  
  ##識別性の問題を回避するために分散共分散行列の対角成分を1を固定する
  gamma <- diag(diag(Cov_hat)^(-1/2))
  Cov <- cov2cor(Cov_hat)
  inv_cov <- solve(Cov)
  #util_mu <- util_mu %*% gamma
  #U <- U %*% gamma
  
  
  #対数尤度関数を計算
  prob <- pnorm(util_mu, 0, 1)
  LL <- sum(y*log(prob) + (1-y)*log(1-prob))
  
  #サンプリング結果を表示
  print(LL)
  print(Cov)
  print(round(rbind(theta1, thetat1), 3))
  print(round(rbind(theta2, thetat2), 3))
  print(round(rbind(diag(tau1), diag(taut1)), 3))
  print(round(rbind(diag(tau2), diag(taut2)), 3))
}

i <- 1
X[item_index[[i]], ]

DT <- array(0, dim=c(s, s*k, f))
DT_T <- array(0, dim=c(s*k, s, f))
for(i in 1:f){
  if(i%%10000==0){
    print(i)
  }
  DT[, , i] <- rbind(as.numeric(t(cbind(X[i, ], 0))), as.numeric(t(cbind(0, X[i, ]))))
  DT_T[, , i] <- t(DT[, , i])
}



xx <- array(0, dim=c(s*k, s*k, f))
for(i in 1:f){
  xx[, , i] <- DT_T[, , i] %*% DT[, , i]
}

for(i in 1:f){
  xx[, , f]
}


XVX <- matrix(0, nrow=s*k, ncol=s*k)
XVY <- matrix(0, nrow=s*k, ncol=1)
for(i in 1:f){
  XV <- DT_T[, , i] %*% inv_cov
  XVX <- XVX + XV %*% DT[, , i]
  XVY <- XVY + XV %*% U[i, ]
}
XVX

t(t(as.numeric(t(X)))) %*% inv_cov[1, ]

solve(XVX) %*% XVY
