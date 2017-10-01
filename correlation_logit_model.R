#####相関のある多項ロジットモデル#####
library(MASS)
library(mlogit)
library(MCMCpack)
library(bayesm)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho <- ifelse(abs(rho) < 0.1, 0, rho)
  rho[upper.tri(rho)] <- 0
  
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  X.Sigma <- eigen(Sigma)
  Lambda <- diag(X.Sigma$values)
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda)
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
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
hh <- 6000   #サンプル数
choise <- 5   #選択可能数
st <- 5   #基準ブランド
k <- 5   #説明変数の数

##説明変数の発生
#通常価格の発生
PRICE <- matrix(runif(hh*choise, 0.6, 1), nrow=hh, ncol=choise) - 1   

#ディスカウント率の発生
DISC <- matrix(runif(hh*choise, 0, 0.4), nrow=hh, ncol=choise) 

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

#カテゴリーロイヤルティ
ROYL <- matrix(rnorm(hh), nrow=hh, ncol=1)

##説明変数をベクトル化
#切片のベクトル化
BP <- matrix(as.numeric(diag(choise)), nrow=hh*choise, ncol=choise, byrow=T)[, -choise]

#ブランドロイヤルティのベクトル化
Royl_vec <- matrix(0, nrow=hh*choise, ncol=choise)
for(i in 1:hh){
  r <- ((i-1)*choise+1):((i-1)*choise+choise)
  Royl_vec[r, ] <- diag(ROYL[i], choise)
}
Royl_vec <- Royl_vec[, -choise]

#その他の変数のベクトル化
Price_vec <- as.numeric(t(PRICE))
Disc_vec <- as.numeric(t(DISC))
Disp_vec <- as.numeric(t(DISP))
Camp_vec <- as.numeric(t(CAMP))

#データを結合
round(X <- data.frame(BP=BP, PRICE=Price_vec, DISC=Disc_vec, DISP=Disp_vec, CAMP=Camp_vec, ROYL=Royl_vec), 2)
XM <- as.matrix(X)

##IDの設定
id <- rep(1:hh, rep(choise, hh))
c <- rep(1:choise, hh)
ID <- data.frame(no=1:(hh*choise), id=id, c=c)


####応答変数を発生####
##分散共分散行列の設定
corM <- corrM(col=choise, lower=0, upper=0.9, eigen_lower=0.01, eigen_upper=0.2)   #相関行列を作成
Sigma <- covmatrix(col=choise, corM=corM, lower=1, upper=1)   #分散共分散行列
Cov <- Sigma$covariance

##パラメータの設定
beta1 <- -5.5   #価格のパラメータ
beta2 <- 5.7   #割引率のパラメータ
beta3 <- 2.0   #特別陳列のパラメータ
beta4 <- 1.8   #キャンペーンのパラメータ
beta5 <- c(1.1, 0.6, -0.5, 0.3)   #カテゴリーロイヤルティのパラメータ
beta0 <- c(0.5, 1.1, 1.4, 2.2)   #ブランド1～4の相対ベース販売力
betat <- c(beta0, beta1, beta2, beta3, beta4, beta5)

##ロジットと確率の計算
logit0 <- matrix(XM %*% betat, nrow=hh, ncol=choise, byrow=T)
logit <- logit0 + mvrnorm(hh, rep(0, choise), Cov)
Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=choise)

##カテゴリカル分布から選択肢を発生
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
colMeans(Y); colSums(Y)


####マルコフ連鎖モンテカルロ法で相関構造のある混合型ロジスティック回帰モデルを推定####
##多項ロジットモデルの対数尤度関数
LLike <- function(logit, X, Y, hh, choise){
  
  d <- rowSums(exp(logit))
  LLl <- rowSums(Y * logit) - log(d)
  LL <- sum(LLl)
  LL_val <- list(LL=LL, LLl=LLl)
  return(LL_val)
}

##初期値の設定用の多項ロジットの対数尤度
loglike <- function(theta, X, Y, hh, choise){
  matrix(XM %*% runif(ncol(XM)), nrow=hh, ncol=choise, byrow=T)
  
  logit <- matrix(XM %*% theta, nrow=hh, ncol=choise, byrow=T)
  
  d <- rowSums(exp(logit))
  LLl <- rowSums(Y * logit) - log(d)
  LL <- sum(LLl)
  return(LL)
}


##MCMCアルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- c()   #対数尤度の保存用

##説明変数を多次元配列化
X.array <- array(0, dim=c(choise, ncol(XM), hh))
for(i in 1:hh){
  X.array[, , i] <- XM[ID$id==i, ]
}
YX.array <- array(0, dim=c(choise, ncol(XM)+1, hh))
id_r <- matrix(1:(hh*choise), nrow=hh, ncol=choise, byrow=T)


##事前分布の設定
nu <- ncol(XM)-2   #逆ウィシャート分布の自由度
V <- nu*diag(choise)   #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(XM))  #回帰係数の平均の事前分布
Adelta <- 0.01 * diag(rep(1, ncol(XM)))   #回帰係数の事前分布の分散

##サンプリング結果の保存用
THETA <- array(0, dim=c(hh, choise, R/keep))
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
SIGMA <- matrix(0, nrow=R/keep, ncol=choise^2)

##初期値の設定
#回帰係数の初期値
beta00 <- rep(0, ncol(XM))
res <- optim(beta00, loglike, gr=NULL, XM, Y, hh, choise, method="BFGS", hessian=TRUE, control=list(fnscale=-1))
oldbeta <- res$par

#分散共分散行列の初期値
oldcov <- diag(choise)
oldcov0 <- oldcov
cov_inv <- solve(oldcov)

#効用の平均構造の初期値
theta_mu <- matrix(XM %*% oldbeta, nrow=hh, ncol=choise, byrow=T)
oldtheta <- matrix(XM %*% oldbeta, nrow=hh, ncol=choise, byrow=T) + mvrnorm(hh, rep(0, choise), oldcov)


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){

  ##MH法でthetaをサンプリング
  thetad <- oldtheta
  thetan <- thetad + 0.5 * mvrnorm(hh, rep(0, choise), diag(choise))
  
  #事前分布の誤差を計算
  er_new <- thetan - theta_mu
  er_old <- thetad - theta_mu
  
  #対数尤度と対数事前分布を計算
  lognew <- LLike(logit=thetan, X=XM, Y=Y, hh=hh, choise=choise)$LLl
  logold <- LLike(logit=thetad, X=XM, Y=Y, hh=hh, choise=choise)$LLl
  logpnew <- apply(er_new, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  logpold <- apply(er_old, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  
  #サンプリングを採択するかどうかを決定
  rand <- runif(hh)   #一様乱数から乱数を発生
  LLind.diff <- exp(lognew + logpnew - logold - logpold)   #採択率を計算
  alpha <- ifelse(LLind.diff > 1, 1, LLind.diff)
  alpha <- matrix(alpha, nrow=hh, ncol=choise)
  
  #alphaに基づきbetaを採択
  oldtheta.r <- ifelse(alpha > rand, thetan, oldtheta)   #alphaがrandを上回っていたら採択
  adopt <- sum(oldtheta[, 1]!=oldtheta.r[, 1])/hh   #採択率
  oldtheta <- oldtheta.r   #パラメータを更新
  
  
  ##多変量回帰モデルのパラメータをサンプリング
  util.v <- as.numeric(t(oldtheta))
  
  ##betaの分布のパラメータの計算とmcmcサンプリング
  #z.vecとX.vecを結合して多次元配列に変更
  YX.bind <- cbind(util.v, XM)
  for(i in 1:hh){
    YX.array[, , i] <- YX.bind[id_r[i, ], ]
  }
  
  ##回帰モデルのギブスサンプリングでbetaとsigmaを推定
  #betaのギブスサンプリング
  XVX <- rowSums(apply(X.array, 3, function(x) t(x) %*% cov_inv %*% x))
  XVY <- rowSums(apply(YX.array, 3, function(x) t(x[, -1]) %*% cov_inv %*% x[, 1]))
  
  #betaの分布の共分散行列
  inv_XVX <- solve(XVX + Adelta)
  
  #betaの分布の平均パラメータ
  B <- inv_XVX %*% (XVY + Adelta %*% Deltabar)   #betaの平均
  b1 <- as.numeric(B)
  oldbeta <- mvrnorm(1, b1, inv_XVX)   #多変量正規分布からbetaをサンプリング
  oldbeta <- betat
  theta_mu <- matrix(XM %*% oldbeta, nrow=hh, ncol=choise, byrow=T)
  
  ##Covの分布のパラメータの計算とmcmcサンプリング
  #逆ウィシャート分布のパラメータを計算
  R.error <- oldtheta - theta_mu
  IW.R <- V + t(R.error) %*% R.error
  
  #逆ウィシャート分布の自由度を計算
  Sn <- nu + hh
  
  #逆ウィシャート分布からCovをサンプリング
  Cov_hat <- rwishart(Sn, solve(IW.R))$IW
  cov_inv <- solve(oldcov)
  
  #識別性のためパラメータに制約をかける
  oldcov <- cov2cor(Cov_hat)
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep   
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- oldcov
    
    print(rp)
    print(round(c(sum(lognew), res$value), 2))
    print(round(rbind(oldbeta, res$par, betat), 3))
    print(round(cbind(oldcov, Cov), 3))
    print(round(adopt, 3))
  }
}

matplot(BETA[, 1:4], type="l")
matplot(SIGMA[, 1:5], type="l")
