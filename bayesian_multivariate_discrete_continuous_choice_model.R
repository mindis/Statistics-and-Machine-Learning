#####誘導型多変量離散・連続モデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(psych)
library(gtools)
library(MNP)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
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
hh <- 1500
items <- 15   #観測ブランド数


####説明変数の発生####
#連続変数の発生
cont <- 2
X.cont <- matrix(rnorm(hh*cont, 0, 1), nrow=hh, ncol=cont)

#二値変数の発生
bin <- 2
X.bin <- matrix(0, nrow=hh, ncol=bin)
for(i in 1:bin){
  X.bin[, i] <- rbinom(hh, 1, runif(1, 0.35, 0.65))
}

#価格とプロモーション説明変数
Price <- matrix(0, nrow=hh, ncol=items)
Disp <- matrix(0, nrow=hh, ncol=items)
freq <- rpois(hh, 5)
freq <- ifelse(freq < 1, 1, freq)

for(i in 1:items){
  lower <- runif(1, 0.4, 0.7)
  upper <- runif(1, 0.8, 1.0)
  Price[, i] <- runif(hh, lower, upper) 
  Disp[, i] <- rbinom(hh, freq, runif(1, 0.2, 0.5))/freq
}


##説明変数をベクトル化
#切片と説明変数の設定
BP.vec <- matrix(0, nrow=hh*items, ncol=items)
X1.vec <- matrix(0, nrow=hh*items, ncol=items*cont)
X2.vec <- matrix(0, nrow=hh*items, ncol=items*bin)

for(i in 1:hh){
  r <- ((i-1)*items+1):((i-1)*items+items)
  BP.vec[r, ] <- diag(items) 
  X1.vec[r, ] <- cbind(diag(X.cont[i, 1], items), diag(X.cont[i, cont], items))
  X2.vec[r, ] <- cbind(diag(X.bin[i, 1], items), diag(X.bin[i, bin], items))
}

#価格とプロモーションの設定
Price.v <- as.numeric(t(Price))
Disp.v <- as.numeric(t(Disp))


#説明変数の結合
X.vec <- data.frame(bp=BP.vec, cont=X1.vec, bin=X2.vec, price=Price.v, disp=Disp.v)
XM.vec <- as.matrix(X.vec)
k1 <- ncol(XM.vec)

##IDの設定
id.v <- rep(1:hh, rep(items, hh))
pd <- rep(1:items, hh)
ID.vec <- data.frame(no=1:(hh*items), id=id.v, pd=pd)


####購買有無を多変量プロビットモデルで発生させる####
#相関行列の設定
corM <- corrM(col=items, lower=-0.55, upper=0.85)   #相関行列を作成
Sigma <- covmatrix(col=items, corM=corM, lower=1, upper=1)   #分散共分散行列
Cov <- Sigma$covariance

#妥当な購買率が出るまでパラメータの設定を繰り返す
for(i in 1:10000){
  print(i)
  
  #パラメータの設定
  beta0 <- runif(items, -0.7, 2.0)
  beta1 <- runif(items*cont, 0, 0.9)
  beta2 <- runif(items*bin, -0.8, 1.0)
  beta3 <- runif(1, -3.5, -2.5)
  beta4 <- runif(1, 2.0, 2.5)
  betat <- c(beta0, beta1, beta2, beta3, beta4)
  
  ##効用関数の計算と応答変数の発生
  #効用関数の計算
  U_mean.vec <- cbind(XM.vec) %*% betat
  error.vec  <- as.numeric(t(mvrnorm(hh, rep(0, items), Cov)))
  U.vec <- U_mean.vec + error.vec
  U.mx <- matrix(U.vec, nrow=hh, ncol=items, byrow=T)
  
  #応答変数の発生
  y1 <- ifelse(U.vec > 0, 1, 0)
  Y1 <- matrix(y1, nrow=hh, ncol=items, byrow=T)
 
  #妥当な購買率が出るまで繰り返す
  if(min(colMeans(Y1)) > 0.1 & max(colMeans(Y1)) < 0.75) break   
}

#発生させた応答変数を集計
colMeans(Y1); colSums(Y1)


####購買があった場合、ポアソン分布から購買数を発生####
U.pois <- ifelse(U.mx < 0, 0, U.mx)   #ポアソン分布の効用を定義

##購買確率と効用の関係
mu <- c()
for(i in 1:items){
  mu <- c(mu, mean(U.pois[U.pois[, i] > 0, i]))
}

round(rbind(u=mu, pr=colMeans(Y1)), 3)


##ポアソン分布の説明変数
#購買個数は家族人数と収入に依存
family <- rpois(hh, 2.5)
income <- rpois(hh, 5.3)

#ゼロを1に置き換え
family <- ifelse(family < 1, 1, family)
income <- ifelse(income < 1, 1, income)

##説明変数をベクトル形式に変更
family.v <- as.numeric(t(matrix(family, nrow=hh, ncol=items)))
income.v <- as.numeric(t(matrix(income, nrow=hh, ncol=items)))
upois.v <- as.numeric(t(U.pois))

#購買が発生していないならincomeとfamilyは0にしておく
income.vec <- income.v * y1
family.vec <- family.v * y1

#説明変数を結合
Z <- data.frame(u=upois.v, income=income.vec, family=family.vec)
ZM <- as.matrix(Z)
k2 <- ncol(ZM)


##ポアソン分布から購入数を発生
#パラメータの設定
theta1 <- runif(1, 0.4, 0.8)
theta2 <- runif(1, 0.05, 0.08)
theta3 <- runif(1, 0.08, 0.13)
thetat <- c(theta1, theta2, theta3)

#lambdaを計算
lambda <- exp(ZM %*% thetat) * y1

#ポアソン分布化から応答変数を発生
y2 <- rep(0, length(lambda))

for(i in 1:length(lambda)){
  print(i)
  if(lambda[i]==0) next   #lambdaが0なら次へ
  
  for(j in 1:1000){
    y2[i] <- rpois(1, lambda[i])
    if(y2[i]==0) {next} else {break}
  }
}

Y2 <- matrix(y2, nrow=hh, ncol=items, byrow=T)   #行列形式に変更


##結果の集計と確認
summary(Y2)

Y2.mean <- c()
for(i in 1:items){
  Y2.mean <- c(Y2.mean, mean(Y2[Y2[, i] > 0, i]))
}

round(rbind(u=mu, pr=colMeans(Y1), buy=Y2.mean), 3)   #効用、購入確率、購買数量の関係


####マルコフ連鎖モンテカルロ法で多変量離散連続モデルを推定####
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

##ポアソン回帰モデルの対数尤度
##変量効果ポアソン回帰モデルの対数尤度
loglike <- function(theta, y, X){
  
  #尤度を定義する
  lambda <- exp(X %*% theta)   #平均構造
  LLi <- y*log(lambda)-lambda - lfactorial(y)   #対数尤度
  LL <- sum(LLi)   #対数尤度の和
  
  #結果を返す
  LL.val <- list(LLi=LLi, LL=LL)
  return(LL.val)
}


####MCMCアルゴリズムの設定####
R <- 20000
sbeta <- 1.5
keep <- 4

##事前分布の設定
nu <- items   #逆ウィシャート分布の自由度
V <- nu*diag(items) + 5   #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(XM.vec))   #多変量プロビットモデルの回帰係数の事前分布の平均
Adelta <- solve(100 * diag(rep(1, k1)))   #多変量プロビットモデルの回帰係数の事前分布の分散
Zetabar <- rep(0, k2)   #ポアソン回帰モデルの回帰係数の事前分布
Azeta <- solve(100 * diag(rep(1, k2)))   #ポアソン回帰モデルの回帰係数の事前分布

##サンプリング結果の保存用配列
Util <- matrix(0, nrow=mcmc/keep, ncol=items)
BETA <- matrix(0, nrow=mcmc/keep, k1)
SIGMA <- matrix(0, nrow=mcmc/keep, ncol=items^2)
THETA <- matrix(0, nrow=R/keep, ncol=k2)

##データの設定
#デザイン行列を多次元配列に変更
X.array <- array(0, dim=c(items, k1, hh))
for(i in 1:hh){
  X.array[, , i] <- XM.vec[ID.vec$id==i, ]
}
YX.array <- array(0, dim=c(items, k1+1, hh))


#インデックスの作成
id_r <- matrix(1:(hh*items), nrow=hh, ncol=items, byrow=T)

##計算用パラメータの格納用
Z <- matrix(0, nrow=hh, ncol=items)   #効用関数の格納用
MVR.U <- matrix(0, nrow=hh, ncol=items)

##初期値の設定
#回帰係数の初期値
##プロビットモデルの対数尤度の定義
probit_LL <- function(x, Y, X){
  #パラメータの設定
  b0 <- x[1]
  b1 <- x[2:(ncol(X)+1)]
  
  #効用関数の定義
  U <- b0 + as.matrix(X) %*% b1
  
  #対数尤度を計算
  Pr <- pnorm(U)   #確率の計算
  LLi <- Y*log(Pr) + (1-Y)*log(1-Pr)
  LL <- sum(LLi)
  return(LL)
}

#応答変数ごとに独立にプロビットモデルを当てはめ初期値を設定
first_beta <- matrix(0, nrow=items, ncol=(k1-2)/items+2)

for(b in 1:items){
  print(b)
  for(i in 1:1000){
    #初期パラメータの設定
    X <- cbind(X.cont, X.bin, Price[, b], Disp[, b])
    x <- c(runif(1, -0.5, 1.0), runif(cont, 0, 1), runif(bin, -0.5, 0.5), runif(1, -3.5, -2.5), runif(1, 2.0, 2.5))

    #準ニュートン法で最大化
    res <- try(optim(x, probit_LL, Y=Y1[, b], X=X, method="BFGS", hessian=FALSE, 
                     control=list(fnscale=-1)), silent=TRUE)
    
    #エラー処理
    if(class(res) == "try-error") {
      next
    } else {
      first_beta[b, ] <- res$par
      break
    }   
  }
}

oldbeta <- c(as.numeric(first_beta[, 1:(1+cont+bin)]), mean(first_beta[, 2+cont+bin]), mean(first_beta[, 3+cont+bin]))
betaf <- oldbeta

#ポアソン回帰の初期値
oldtheta <- c(runif(1, 0.4, 0.8), runif(2, 0, 0.1))
thetaf <- oldtheta

#分散共分散行列の初期値
corf <- corrM(col=items, lower=0, upper=0)   #相関行列を作成
Sigmaf <- covmatrix(col=items, corM=corf, lower=1, upper=1)   #分散共分散行列
oldcov <- Sigmaf$covariance

#潜在効用の初期値
old.utilm<- matrix(XM.vec %*% oldbeta, nrow=hh, ncol=items, byrow=T)   #潜在効用の平均構造
Z <- old.utilm + mvrnorm(hh, rep(0, items), oldcov)   #平均構造+誤差

#切断正規分布の切断領域を定義
a <- ifelse(Y1==0, -100, 0)
b <- ifelse(Y1==1, 100, 0)


####マルコフ連鎖モンテカルロ法で多変量離散連続モデルを推定####
for(rp in 1:R){

  ##データ拡大法で潜在効用を発生
  #選択結果と整合的な潜在効用を発生させる
  #潜在効用を計算
  old.utilm<- matrix(XM.vec %*% oldbeta, nrow=hh, ncol=items, byrow=T)   #潜在効用の平均構造
  
  #切断正規分布より潜在効用を発生
  MVR.S <- c()
  for(j in 1:items){
    MVR <- cdMVN(old.utilm, oldcov, j, Z)   #多変量正規分布の条件付き分布を計算
    MVR.U[, j] <- MVR$CDmu   #条件付き平均を抽出
    MVR.S <- c(MVR.S, MVR$CDvar)   #条件付き分散を抽出
    Z[, j] <- rtnorm(MVR.U[, j], sqrt(MVR.S[j]), a[, j], b[, j])   #切断正規分布より潜在変数を抽出
  }
  
  Z[is.infinite(Z)] <- 0
  Z[is.nan(Z)] <- 0
  Zold <- Z
  
  ##betaの分布のパラメータの計算とmcmcサンプリング
  #z.vecとX.vecを結合して多次元配列に変更
  z.vec <- as.numeric(t(Zold))   #潜在効用をベクトルに変更
  YX.bind <- cbind(z.vec, XM.vec)
  
  for(i in 1:hh){
    YX.array[, , i] <- YX.bind[id_r[i, ], ]
  }
  
  #betaの平均パラメータを計算
  invcov <- solve(oldcov)
  xvx.vec <- rowSums(apply(X.array, 3, function(x) t(x) %*% invcov %*% x))
  XVX <- matrix(xvx.vec, nrow=ncol(XM.vec), ncol=ncol(XM.vec), byrow=T)
  XVY <- rowSums(apply(YX.array, 3, function(x) t(x[, -1]) %*% invcov %*% x[, 1]))
  
  #betaの分布の分散共分散行列パラメータ
  inv_XVX <- solve(XVX + Adelta)
  
  #betaの分布の平均パラメータ
  B <- inv_XVX %*% (XVY + Adelta %*% Deltabar)   #betaの平均
  
  #多変量正規分布からbetaをサンプリング
  oldbeta <- mvrnorm(1, as.numeric(B), inv_XVX)
  
  ##Covの分布のパラメータの計算とmcmcサンプリング
  #逆ウィシャート分布のパラメータを計算
  R.error <- matrix(z.vec - XM.vec %*% oldbeta, nrow=hh, ncol=items, byrow=T)
  R <- solve(V) + matrix(rowSums(apply(R.error, 1, function(x) x %*% t(x))), nrow=items, ncol=items, byrow=T)
  
  #逆ウィシャート分布の自由度を計算
  Sn <- nu + hh
  
  #逆ウィシャート分布からCovをサンプリング
  Cov_hat <- try(rwishart(Sn, solve(R))$IW, silent=TRUE)
  if(class(Cov_hat)=="try-error") {Cov_hat <- oldcov}
  

  #分散共分散行列の対角成分を1に固定する
  oldcov <- cov2cor(Cov_hat)
  
  ##効用関数から購買量をベイジアンポアソン回帰モデルで結びつける
  #ランダムウォークサンプリング
  newtheta <- oldtheta + rnorm(ncol(ZM), 0, 0.025) * c(0.5, 0.05, 0.1)
  ZM[, 1] <- z.vec
  
  
  #対数尤度と対数事前分布を計算
  lognew <- loglike(theta=newtheta, y=y2, X=ZM)$LL
  logold <- loglike(theta=oldtheta, y=y2, X=ZM)$LL
  logpnew <- lndMvn(newtheta, Zetabar, Azeta)
  logpold <- lndMvn(oldtheta, Zetabar, Azeta)
  
  #MHサンプリング
  alpha <- min(1, exp(lognew + logpnew - logold - logpold))
  if(alpha == "NAN") alpha <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しい固定効果betaを採択
  if(u < alpha){
    oldtheta <- newtheta
    logl <- lognew
    
    #そうでないなら固定効果betaを更新しない
  } else {
    logl <- logold
  }
  
  ##サンプリング結果を保存
  mkeep <- rp/keep
  if(rp%%keep==0){
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    THETA[mkeep, ] <- oldtheta
    print(rp)
    print(round(rbind(oldbeta, betat), 2))
    print(round(rbind(oldtheta, thetat), 3))
    print(round(cbind(oldcov[1:8, 1:8], Cov[1:8, 1:8]), 2))
    print(round(logl, 3))
  }
}

