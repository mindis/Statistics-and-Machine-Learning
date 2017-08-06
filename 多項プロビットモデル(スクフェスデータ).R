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

####データクレンジング####
##データの読み込み
sf_data <- read.csv("スクフェスデータセット.csv")
hh <- nrow(sf_data)   #サンプル数


##基準メンバーの設定
sf_data$パートナー <-  as.character(sf_data$パートナー)
table(sf_data$パートナー)
unique(sf_data$パートナー)

#mobを定義する
sf_data$パートナー <- ifelse(sf_data$パートナー %in% c("あんじゅ", "いるか", "ココ", "るう", "文絵", "咲良", 
                                             "聖闘士", "遊宇", "ココ", "いるか", "小雪"), "モブ", sf_data$パートナー)
#aqoursメンバーを定義
sf_data$パートナー <- ifelse(sf_data$パートナー %in% c("ダイヤ", "果南", "花丸", "ルビィ", "鞠莉", "善子",
                                             "千歌", "曜", "梨子"), "aqours", sf_data$パートナー)

#基準メンバーを設定
choise <- matrix(as.numeric(table(1:nrow(sf_data), sf_data$パートナー)), 
                            nrow=nrow(sf_data), ncol=length(unique(sf_data$パートナー)))
colnames(choise) <- colnames(table(1:nrow(sf_data), sf_data$パートナー))
choise.m <- choise[, c("穂乃果", "ことり", "海未", "凛", "花陽", "真姫", "絵里", "希", "にこ", "aqours", "モブ")]

index.choise <- subset(1:nrow(choise.m), choise.m[, ncol(choise.m)]==1)
choise.m <- choise.m[-index.choise, -ncol(choise.m)]


#応答変数の設定
Y <- choise.m
y <- as.numeric(Y %*% 1:ncol(Y))
member <- length(unique(as.numeric(y)))
hh <- length(y)

##回帰モデルを推定するために説明変数をベクトル形式に変更設定
#切片の設定
p <- c(1, rep(0, member-1))
Pop <- matrix(p, nrow=hh*member, ncol=member-1, byrow=T)
POP <- subset(Pop, rowSums(Pop) > 0)

#レベルと順位の説明変数の設定
LV <- scale(sf_data$lv[-index.choise])
SCORE <- -scale(sf_data$順位[-index.choise])

#ベクトル形式に変更
LV.v <- matrix(0, hh*(member-1), ncol=member-1)
SCORE.v <- matrix(0, hh*(member-1), ncol=member-1)

for(i in 1:hh){
  index.v <- ((i-1)*(member-1)+1):((i-1)*(member-1)+member-1)
  v.lv <- diag(LV[i, ], member-1)
  v.score <- diag(SCORE[i, ], member-1)
  LV.v[index.v, ] <- v.lv 
  SCORE.v[index.v, ] <- v.score
}
X <- data.frame(POP, lv=LV.v, score=SCORE.v)   #説明変数
XM <- as.matrix(X)

#IDの設定
c <- rep(1:(member-1), hh)
id <- rep(1:hh, rep(member-1, hh))
ID <- data.frame(no=1:nrow(XM), id=id, c=c)
id_r <- matrix(1:nrow(XM), nrow=hh, ncol=member-1, byrow=T)


####マルコフ連鎖モンテカルロ法で多項プロビットモデルを推定####
####MCMC推定のための推定準備####

##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##多変量正規分布の条件付き期待値と分散を計算する関数
cdMVN <- function(mu, Cov, dependent, U){
  
  #分散共分散行列のブロック行列を定義
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent]
  Cov21 <- Cov[-dependent, dependent]
  Cov22 <- Cov[-dependent, -dependent]
  
  #条件付き分散と条件付き平均を計算
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mu[, dependent] + t(CDinv %*% t(U[, -dependent] - mu[, -dependent]))   #条件付き平均を計算
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #条件付き分散を計算
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##アルゴリズムの設定
R <- 100000
sbeta <- 1.5
keep <- 10
llike <- array(0, dim=c(R/keep))   #対数尤度の保存用


##説明変数を多次元配列化
X.array <- array(0, dim=c(member-1, ncol(XM), hh))
for(i in 1:hh){
  X.array[, , i] <- XM[ID$id==i, ]
}
YX.array <- array(0, dim=c(member-1, ncol(XM)+1, hh))


#推定プロセスの格納配列
UM <- matrix(0, nrow=hh, ncol=member-1)
util.M <- matrix(0, nrow=hh, ncol=member-1)   

##事前分布の設定
nu <- member   #逆ウィシャート分布の自由度
V <- solve((1/10)*diag(member-1))    #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(XM))  #回帰係数の平均の事前分布
Adelta <- solve(100 * diag(rep(1, ncol(XM))))   #回帰係数の事前分布の分散

##サンプリング結果の保存用配列
Util <- array(0, dim=c(hh, member-1, R/keep))
BETA <- matrix(0, nrow=R/keep, ncol(XM))
SIGMA <- matrix(0, nrow=R/keep, ncol=(member-1)^2)

##初期値の設定
#回帰係数の初期値
oldbeta <- c(colSums(Y[, -member])/sum(Y[, -member])*10, runif(2*(member-1), -1, 1))

#分散共分散行列の初期値
corM.f <- corrM(col=member-1, lower=0, upper=0)   #相関行列を作成
Sigma.f <- covmatrix(col=member-1, corM=corM.f, lower=1, upper=1)   #分散共分散行列
oldcov <- Sigma.f$covariance

#効用の平均構造の初期値
old.utilm <- matrix(XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)


#効用の初期値
old.util <- old.utilm + mvrnorm(nrow(old.utilm), rep(0, member-1), oldcov)


####マルコフ連鎖モンテカルロ法で多項プロビットモデルを推定####
for(rp in 1:R){
  
  ##選択結果と整合的な潜在効用を発生させる
  #条件付き期待値と条件付き分散を計算
  S <- rep(0, member-1)
  
  for(j in 1:(member-1)){
    MVR <- cdMVN(mu=old.utilm, Cov=oldcov, dependent=j, U=old.util)   #条件付き分布を計算
    UM[, j] <- MVR$CDmu   #条件付き期待値を取り出す
    S[j] <- sqrt(MVR$CDvar)    #条件付き分散を取り出す
    
    #潜在変数を発生させる
    #切断領域の設定
    max.u <- apply(cbind(old.util[, -j], 0), 1, max)
    max.u <- ifelse(y==member, 0, max.u)

    #切断正規分布より潜在変数を発生
    old.util[, j] <- ifelse(y==j, rtnorm(mu=UM[, j], sigma=S[j], a=max.u, b=100), 
                            rtnorm(mu=UM[, j], sigma=S[j], a=-100, b=max.u))
    old.util[, j] <- ifelse(is.infinite(old.util[, j]), ifelse(y==j, max.u + runif(1), max.u - runif(1)), old.util[, j])
  }
  util.v <- as.numeric(t(old.util))
  
  ##betaの分布のパラメータの計算とmcmcサンプリング
  #z.vecとX.vecを結合して多次元配列に変更
  YX.bind <- cbind(util.v, XM)
  for(i in 1:hh){
    YX.array[, , i] <- YX.bind[id_r[i, ], ]
  }
  
  ##回帰モデルのギブスサンプリングでbetaとsigmaを推定
  #betaのギブスサンプリング
  invcov <- solve(oldcov)
  xvx.vec <- rowSums(apply(X.array, 3, function(x) t(x) %*% invcov %*% x))
  XVX <- matrix(xvx.vec, nrow=ncol(XM), ncol=ncol(XM), byrow=T)
  XVY <- rowSums(apply(YX.array, 3, function(x) t(x[, -1]) %*% invcov %*% x[, 1]))
  
  #betaの分布の分散共分散行列のパラメータ
  inv_XVX <- solve(XVX + Adelta)
  
  #betaの分布の平均パラメータ
  B <- inv_XVX %*% (XVY + Adelta %*% Deltabar)   #betaの平均
  b1 <- as.numeric(B)
  
  #多変量正規分布から回帰係数をサンプリング
  oldbeta <- mvrnorm(1, b1, inv_XVX)
  
  ##Covの分布のパラメータの計算とmcmcサンプリング
  #逆ウィシャート分布のパラメータを計算
  R.error <- matrix(util.v - XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)
  IW.R <- V + matrix(rowSums(apply(R.error, 1, function(x) x %*% t(x))), nrow=member-1, ncol=member-1, byrow=T)
  
  #逆ウィシャート分布の自由度を計算
  Sn <- nu + hh
  
  #逆ウィシャート分布からCovをサンプリング
  Cov_hat <- rwishart(Sn, solve(IW.R))$IW
  oldcov <- cov2cor(Cov_hat)
  
  ##潜在効用とパラメータを更新
  #潜在効用と潜在効用の平均を更新
  old.utilm <- matrix(XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)

  ##サンプリング結果を保存
  if(rp%%keep==0){
    print(rp)
    mkeep <- rp/keep
    Util[, , mkeep] <- old.util
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    print(round(oldcov, 2))
    print(round(oldbeta, 2))
  }
}

round(matrix(colMeans(SIGMA[burnin:nrow(SIGMA), ]), nrow=member-1, ncol=member-1), 3)

c("穂乃果", "ことり", "海未", "凛", "花陽", "真姫", "絵里", "希", "にこ", "モブ")

####関数で推定####
Data1 <- list(p=member, y=y, X=XM)
Mcmc1 <- list(R=200000, keep=10)

#多項プロビットモデルを推定
out <- rmnpGibbs(Data=Data1,Mcmc=Mcmc1)
BETA.out <- out$betadraw
SIGMA.out <- out$sigmadraw


####推定結果の要約と適合度の確認####
burnin <- 5000/keep   #バーンイン期間

##サンプリング結果を可視化
#回帰係数のプロット
matplot(BETA[, 1:3], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 4:6], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 7:9], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 10:12], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 13:15], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 16:18], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 19:21], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 22:24], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 25:27], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")

#分散供分散行列の可視化
matplot(SIGMA[, 1:4], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 5:8], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 9:12], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 13:16], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")

sigma11 <- matrix(SIGMA.out[, 1], nrow=nrow(SIGMA.out), ncol=ncol(SIGMA.out))
matplot(SIGMA.out[, 1:4]/sigma11[, 1:4], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA.out[, 5:8]/sigma11[, 5:8], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA.out[, 9:12]/sigma11[, 9:12], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA.out[, 13:16]/sigma11[, 13:16], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")

sigma <- matrix(0, nrow=nrow(SIGMA.out), ncol=ncol(SIGMA.out))
for(i in 1:nrow(SIGMA.out)){
  print(i)
  sigma[i, ] <- as.numeric(cov2cor(matrix(SIGMA.out[i, ], nrow=member-1, ncol=member-1)))
}
matplot(sigma[, 1:4], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(sigma[, 5:8], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(sigma[, 9:12], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(sigma[, 13:16], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(sigma[, 17:20], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")



##推定値の事後平均の比較
#betaの要約統計量
round(colMeans(BETA.out[burnin:nrow(BETA.out), ] / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(関数推定)の事後平均
round(colMeans(BETA[burnin:nrow(BETA), ]), 3)   #betaの事後平均
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(BETA[burnin:nrow(BETA), ], 2, sd), 2)   #事後標準偏差

#sigmaの要約統計量
round(colMeans(SIGMA.out[burnin:nrow(SIGMA.out), ]  / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(関数推定)の事後平均
round(colMeans(SIGMA[burnin:nrow(SIGMA), ]), 3)   #betaの事後平均
round(as.numeric(Cov), 3)   #真の値
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, sd), 2) #事後標準偏差

