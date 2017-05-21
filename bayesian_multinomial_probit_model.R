#####多項プロビットモデル#####
library(MASS)
library(bayesm)
library(condMVNorm)
library(MCMCpack)
library(glmm)
library(lme4)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####多変量正規乱数を発生させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

#多変量回帰モデルの相関行列を作成
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
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
#set.seed(8437)
##データの設定
hh <- 500   #プレイヤー数
choise <- 5   #選択可能数
st <- 5   #基準ブランド
k <- 5   #回帰係数の数

##説明変数の発生
#通常価格の発生
PRICE <- matrix(runif(hh*choise, 0.7, 1), nrow=hh, ncol=choise)   

#ディスカウント率の発生
DISC <- matrix(runif(hh*choise, 0, 0.3), nrow=hh, ncol=choise)

#特別陳列の発生
DISP <- matrix(0, nrow=hh, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hh, 1, r)
}

#特別キャンペーンの発生
CAMP <- matrix(0, nrow=hh, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hh, 1, r)
}

##分散共分散行列の設定
corM <- corrM(col=choise-1, lower=-0.55, upper=0.60)   #相関行列を作成
Sigma <- covmatrix(col=choise-1, corM=corM, lower=2, upper=4)   #分散共分散行列
Cov <- Sigma$covariance

##パラメータの設定
beta1 <- -7.0   #価格のパラメータ
beta2 <- 8.2   #割引率のパラメータ
beta3 <- 2.0   #特別陳列のパラメータ
beta4 <- 1.8   #キャンペーンのパラメータ
beta0 <- c(0.5, 1.1, 1.4, 2.2)   #ブランド1～4の相対ベース販売力
betat <- c(beta0, beta1, beta2, beta3, beta4)

##相対効用を発生し、選択されたブランドを決定
#基準ブランドとの相対説明変数
PRICE.r <- PRICE[, -5] - PRICE[, 5]
DISC.r <- DISC[, -5] - DISC[, 5]
DISP.r <- DISP[, -5] - DISP[, 5]
CAMP.r <- CAMP[, -5] - CAMP[, 5]

#相対効用の発生
U.mean <- matrix(0, nrow=hh, ncol=choise-1)
for(b in 1:length(beta0)){
  U.mean[, b] <- beta0[b] + PRICE.r[, b]*beta1 + DISC.r[, b]*beta2 + DISP.r[, b]*beta3 + CAMP.r[, b]*beta4
}
U <- U.mean + mvrnorm(n=hh, rep(0, 4), Cov)   #誤差構造を加えた効用
round(cbind(U.mean, U), 2)
U <- t(apply(U.mean, 1, function(x) mvrnorm(1, x, Cov)))

#効用最大化原理に基づき選択ブランドを決定
Y <- apply(U, 1, function(x) ifelse(max(x) < 0, 5, which.max(x)))

#購買を0、1行列に変更
BUY <- matrix(0, hh, choise)
for(i in 1:hh){
  BUY[i, Y[i]] <- 1
}

table(Y)   #選択ブランドの集計
round(data.frame(Y, U=U), 1)   #効用と選択ブランドを比較


####マルコフ連鎖モンテカルロ法で多項プロビットモデルを推定####
####MCMC推定のための推定準備####

##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##多変量正規分布の条件付き期待値と分散を計算する関数
MVR.E <- function(Cov, U.mean, U, choise){
  UE <- matrix(0, hh, choise-1)
  Sig.E <- c()
  
  for(b in 1:(choise-1)){
    #分散共分散行列のブロック行列を定義
    Cov11 <- Cov[b, b]
    Cov22 <- Cov[-b, -b]
    Cov12 <- Cov[b, -b]
    Cov21 <- Cov[-b, b]
    
    #条件付き分散を計算
    Sig.E <- c(Sig.E, Cov11 - Cov12 %*% solve(Cov22) %*% Cov21)
    
    #条件付き期待値を計算
    CovE <- Cov12 %*% solve(Cov22)
    UE.b <- U.mean[, b] - t(CovE %*% t((U[, -b] - U.mean[, -b])))

    UE[, b] <- UE.b
  }
  val <- list(UE=UE, Sig.E=Sig.E)
  return(val)
}
UE <- MVR.E(Cov=Cov, U.mean=U.mean, U=U, choise=choise)$UE
round(data.frame(E=UE, M=U.mean, U=U), 1)   #条件付き期待値、効用の平均構造、効用を比較

##アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 2
llike <- array(0, dim=c(R/keep))   #対数尤度の保存用

##回帰モデルを推定するために説明変数をベクトル形式に変更設定
#切片の設定
p <- c(1, rep(0, choise-1))
bp <- matrix(p, nrow=hh*choise, ncol=choise-1, byrow=T)
BP <- subset(bp, rowSums(bp) > 0)

#説明変数の設定
PRICE.v <- as.numeric(t(PRICE.r))
DISC.v <- as.numeric(t(DISC.r))
DISP.v <- as.numeric(t(DISP.r))
CAMP.v <- as.numeric(t(CAMP.r))

X <- data.frame(BP=BP, PRICE.v, DISC.v, DISP.v, CAMP.v)   #データの結合
XM <- as.matrix(X)

##事前分布の設定
nu <- 5   #逆ウィシャート分布の自由度
V <- df*diag(choise-1) + 10   #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(X))  #回帰係数の平均の事前分布
Adelta <- 100 * diag(rep(1, ncol(X)))   #回帰係数の事前分布の分散

##サンプリング結果の保存用配列
Util <- array(0, dim=c(hh, choise-1, R/keep))
BETA <- matrix(0, nrow=R/keep, length(beta0)+k-1)
SIGMA <- array(0, dim=c(choise-1, choise-1, R/keep))

##初期値の設定
#回帰係数の初期値
oldbeta <- c(runif(choise-1, 0, 3), -3.0, 3.0, runif(2, 0, 2))   

#分散共分散行列の初期値
corM.f <- corrM(col=choise-1, lower=-0.6, upper=0.6)   #相関行列を作成
Sigma.f <- covmatrix(col=choise-1, corM=corM.f, lower=4, upper=9)   #分散共分散行列
oldcov <- Sigma.f$covariance

#効用の平均構造の初期値
b <- oldbeta[choise:length(oldbeta)]
brand_power <- matrix(oldbeta[1:(choise-1)], hh, choise-1, byrow=T)
old.utilm <- brand_power + PRICE.r*b[1] + DISC.r*b[2] + DISP.r*b[3]  + CAMP.r*b[4]

#効用の初期値
old.util <- old.utilm + mvrnorm(nrow(old.utilm), rep(0, 4), oldcov)
old.util <- t(apply(old.utilm, 1, function(x) mvrnorm(1, x, oldcov)))

cov2cor(oldcov)
cov2cor(var(old.util - old.utilm))

####マルコフ連鎖モンテカルロ法で多項プロビットモデルを推定####
for(rp in 1:R){
  
  ##選択結果と整合的な潜在効用を発生させる
  #条件付き期待値と条件付き分散を計算
  MVR <- MVR.E(Cov=oldcov, U.mean=old.utilm, U=old.util, choise=choise)   #条件付き分布を計算
  UM <- MVR$UE   #条件付き期待値を取り出す
  S <- MVR$Sig.E   #条件付き分散を取り出す
  SD <- sqrt(S)
  
  #切断正規分布より潜在効用の発生
  util.M <- matrix(0, nrow=hh, ncol=choise-1)   #潜在効用の格納用行列
  
  #サンプル数分潜在変数を発生させる
  for(i in 1:hh){
    u <- old.util[i, ]
    
    for(b in 1:(choise-1)){
      if(Y[i]!=choise){
        
        #選択ブランドが5以外の場合の潜在効用の発生
        if(Y[i]!=b){
          u[b] <- rtrun(UM[i, b], SD[b], -100, max(c(u[-b], 0)))
          } else {
          u[b] <- rtrun(UM[i, b], SD[b], max(c(u[-b], 0)), 100)
          }
      
        #選択ブランドが5の場合の潜在効用の発生
      } else {
        c <- rtrun(UM[i, b], SD[b], -100, 0)
        u[b] <- c
      }
    }
    util.M[i, ] <- u
  }
  
  ##回帰モデルのギブスサンプリングでbetaとsigmaを推定
  u.vec <- as.numeric(t(util.M))   #潜在効用をベクトルに変更
  ##betaのギブスサンプリング
  oldcovi <- solve(oldcov)
  I <- diag(hh)
  SIGMA.B <- kronecker(I, oldcovi)
  
  
  #回帰係数の平均構造
  B <- solve(t(XM) %*% SIGMA.B %*% XM) %*% t(XM) %*% SIGMA.B %*% u.vec   #回帰係数の最小二乗推定量
  XVX <- t(XM) %*% SIGMA.B %*% XM
  BETA.M <- solve(XVX + solve(Adelta)) %*% (XVX %*% B + solve(Adelta) %*% Deltabar)
  
  #回帰係数の分散共分散行列
  BETA.SIG <- solve(XVX + solve(Adelta))
  
  #多変量正規分布から回帰係数をサンプリング
  oldbeta <- mvrnorm(1, as.numeric(BETA.M), BETA.SIG)
  
  ##sigmaのギブスサンプリング
  #逆ウィシャート分布の自由度を計算
  Sn <- nu + hh
  
  #逆ウィシャート分布のパラメータを計算
  #二乗誤差を計算して和を取る
  Vi <- solve(V) 
  EE <- matrix(0, nrow=choise-1, ncol=choise-1)

  for(i in 1:hh){
    r <- (i-1)*(choise-1)
    u.vech <- u.vec[(r+1):(r+4)]
    XMh <- XM[(r+1):(r+4), ]
    ee <- (u.vech - XMh %*% oldbeta) %*% t(u.vech - XMh %*% oldbeta)
    EE <- EE + ee
  }
  R <- solve(EE + Vi)
  
  #逆ウィシャート乱数を発生
  oldcov <- rwishart(Sn, R)$IW
  
  sigma11 <- oldcov[1, 1]
  
  ##潜在効用とパラメータを更新
  #潜在効用を更新
  old.util <- util.M
  
  #潜在効用の平均構造を更新
  b <- oldbeta[choise:length(oldbeta)]
  brand_power <- matrix(oldbeta[1:(choise-1)], hh, choise-1, byrow=T)
  old.utilm <- brand_power + PRICE.r*b[1] + DISC.r*b[2] + DISP.r*b[3]  + CAMP.r*b[4]
  
  print(round(c(rp, B), 2))
  
  ##サンプリング結果を保存
  mkeep <- rp/keep
  if(rp%%keep==0){
    Util[, , mkeep] <- util.M
    BETA[mkeep, ] <- oldbeta
    SIGMA[, , mkeep] <- oldcov
    #print(round(THETA[mkeep, 1:20], 2))
  }
}
round(cbind(Y, old.utilm, old.util, U), 2)

cov2cor(Cov)
cov2cor(oldcov)

library(MNP)
mnp
