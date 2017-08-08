#####ランクプロビットモデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(gtools)
library(MNP)
library(reshape2)
library(dplyr)
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

#factor型を文字列に変更
sf_data$X1位 <- as.character(sf_data$X1位)
sf_data$X2位 <- as.character(sf_data$X2位)
sf_data$X3位 <- as.character(sf_data$X3位)
r <- c(sf_data$X1位, sf_data$X2位, sf_data$X3位)

table(c(sf_data$X1位, sf_data$X2位, sf_data$X3位))
unique(c(sf_data$X1位, sf_data$X2位, sf_data$X3位))

##モブとaqoursを1つの変数にまとめる
r1 <- ifelse(r %in% c("いるか", "姫乃", "咲", "咲良", "ココ", "かさね", "遊宇", "理亜", "紗菜", "あきる", "ななか",
                      "咲夜", "仁美", "千鶴子", "小雪", "瑞希", "優理", "クリスティーナ"), "モブ", r)

r2 <- ifelse(r1 %in% c("ダイヤ", "鞠莉", "果南", "ルビィ"), "aqours", r1)
rank <- matrix(r2, nrow=nrow(sf_data), ncol=3)


##モブとaqoursが同時に選ばれているサンプルは除去
#aqoursの重複を削除
index1 <- subset(1:nrow(rank), rowSums(rank=="aqours") > 1)

if(length(index1)==0) {
  print("重複なし") } else {
    rank <- rank[-index1, ]
    sf_data1 <- sf_data[-index1, ]
  }

#モブの重複を削除
index2 <- subset(1:nrow(rank), rowSums(rank=="モブ") > 1)

if(length(index2)==0) {
  print("重複なし") } else {
    rank <- rank[-index2, ]
    sf_data1 <- sf_data[-index2, ]
  }

##メンバーごとに整理番号をつける
t(t(table(rank)))
Rank <- rank   #新しいランクデータ
lovelive <- c("穂乃果", "ことり", "海未", "凛", "花陽", "真姫", "にこ", "希", "絵里", "千歌", "曜", "梨子", 
              "花丸", "善子", "aqours", "モブ")

#メンバーごとに番号をつける
for(i in 1:length(lovelive)){
  Rank[Rank==lovelive[i]] <- i
}

Rank <- apply(Rank, 2, as.numeric)   #文字列になっている番号を数値型に変更
data.frame(Rank, rank)   #データを確認

member <- length(lovelive)   #メンバー数
hh <- nrow(Rank)   #サンプル数


##説明変数をベクトル形式のデータフォーマットに変更
#IDを設定
id <- rep(1:hh, rep(member-1, hh))
m <- rep(1:(member-1), hh)
ID <- data.frame(no=1:length(id), id=id, m=m)

#切片の設定
p <- c(1, rep(0, member-1))
Pop <- matrix(p, nrow=hh*member, ncol=member-1, byrow=T)
Pop <- subset(Pop, rowSums(Pop) > 0)

#レベルと順位の説明変数の設定
LV <- scale(sf_data1$lv)
SCORE <- -scale(sf_data1$順位)

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

##データを結合
X <- data.frame(pop=Pop, lv=LV.v)
XM <- as.matrix(X)

####マルコフ連鎖モンテカルロ法でランクプロビットモデルを推定####
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
R <- 40000
sbeta <- 1.5
keep <- 4
llike <- array(0, dim=c(R/keep))   #対数尤度の保存用


##データの設定
#説明変数を多次元配列化
X.array <- array(0, dim=c(member-1, ncol(XM), hh))
for(i in 1:hh){
  X.array[, , i] <- XM[ID[, 2]==i, ]
}
YX.array <- array(0, dim=c(member-1, ncol(XM)+1, hh))

#IDの設定
id_r <- matrix(1:nrow(XM), nrow=hh, ncol=member-1, byrow=T)

##推定プロセスの格納配列
UM <- matrix(0, nrow=hh, ncol=member-1)
util.M <- matrix(0, nrow=hh, ncol=member-1)   

##事前分布の設定
nu <- member   #逆ウィシャート分布の自由度
V <- solve(0.1*diag(member-1))    #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(XM))  #回帰係数の平均の事前分布
Adelta <- solve(100 * diag(rep(1, ncol(XM))))   #回帰係数の事前分布の分散
inv.facov <- solve(diag(100, factors))
alpha_d <- 1
beta_d <- 100

##サンプリング結果の保存用配列
Util <- array(0, dim=c(hh, member-1, R/keep))
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
SIGMA <- matrix(0, nrow=R/keep, ncol=(member-1)^2)
FA.A <- matrix(0, nrow=R/keep, ncol=(member-1)*factors)
FA.D <- matrix(0, nrow=R/keep, ncol=member-1)
FA.F <- array(0, dim=c(hh, factors, R/keep))


##初期値の設定
#回帰係数の初期値
oldbeta <- c((table(Rank)/sum(table(Rank))*10)[-member], runif(member-1, -0.3, 0.3))

#分散共分散行列の初期値
corM.f <- corrM(col=member-1, lower=0, upper=0)   #相関行列を作成
Sigma.f <- covmatrix(col=member-1, corM=corM.f, lower=1, upper=1)   #分散共分散行列
oldcov <- Sigma.f$covariance

#効用の平均構造の初期値
old.utilm <- matrix(XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)

#効用の初期値
old.util <- old.utilm + mvrnorm(nrow(old.utilm), rep(0, member-1), oldcov)

#因子負荷量と独自因子の初期値
A <- matrix(runif((member-1)*factors, -1, 1), nrow=member-1, ncol=factors)
D <- diag(runif(member-1, 0, 0.5))


####マルコフ連鎖モンテカルロ法でランクプロビットモデルを推定####
for(rp in 20336:R){
  
  ##順位選択結果と整合的な潜在効用を発生させる
  #条件付き期待値と条件付き分散を計算
  S <- rep(0, member-1)
  
  for(j in 1:(member-1)){
    MVR <- cdMVN(mu=old.utilm, Cov=oldcov, dependent=j, U=old.util)
    UM[, j] <- MVR$CDmu   #条件付き期待値を取り出す
    S[j] <- sqrt(MVR$CDvar)   #条件付き分散を取り出す
    
    #潜在変数を発生させる
    #切断領域の設定
    rank.u <- t(apply(cbind(old.util[, -j], 0), 1, function(x) sort(x, decreasing=TRUE)))[, 1:3]
    rank.u <- ifelse(Rank==member, 0, rank.u)
    
    #切断正規分布より潜在変数を発生
    old.util[, j] <- ifelse(Rank[, 1]==j, rtnorm(mu=UM[, j], S[j], a=rank.u[, 1], b=150), 
                            ifelse(Rank[, 2]==j, rtnorm(mu=UM[, j], S[j], a=rank.u[, 2], b=rank.u[, 1]),
                                   ifelse(Rank[, 3]==j, rtnorm(mu=UM[, j], S[j], a=rank.u[, 3], b=rank.u[, 2]),
                                          rtnorm(mu=UM[, j], sigma=S[j], a=-150, b=rank.u[, 3]))))
  }
  
  util.v <- as.numeric(t(old.util))   #発生させた潜在効用をベクトルに置き換える
  
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
  Z <- old.util - old.utilm 
  Z <- scale(Z)
  
  ##潜在効用の誤差項から因子分析モデルを推定
  #多変量正規分布から潜在変数f(共通因子)をサンプリング
  ADA <- t(A) %*% solve(A %*% t(A) + D)
  F_mean <- Z %*% t(ADA)   #共通因子の平均
  F_var <- diag(factors) - ADA %*% A    #共通因子の分散共分散行列
  Fi <- t(apply(F_mean, 1, function(x) mvrnorm(1, x, F_var)))   #多変量正規分布から共通因子をサンプリング 
  
  
  #ガンマ分布から独自因子dをサンプリング
  Z.error <- Z - Fi %*% t(A)
  Zv.R <- matrix(rowSums(apply(Z.error, 1, function(x) x %*% t(x))), nrow=member-1, ncol=member-1)
  
  gamma_alpha <- (hh + alpha_d)/2   #alphaを計算
  gamma_beta <- (diag(Zv.R) + beta_d)/2   #betaを計算
  D <- diag(rgamma(length(gamma_beta), gamma_alpha, gamma_beta))   #ガンマ分布から独自因子をサンプリング
  
  #多変量正規分布から因子負荷量Aをサンプリング
  FF <- t(Fi) %*% Fi
  FZ <- t(Fi) %*% Z
  d_sigma <- list()
  
  for(i in 1:(choise-1)){
    d_sigma[[i]]  <- 1/diag(D)[i] * inv.facov
    A_mu <- solve(d_sigma[[i]] + FF) %*% FZ[, i]
    A_cov <- solve(inv.facov + diag(D)[i]*FF) 
    A[i, ] <- mvrnorm(1, A_mu, A_cov)
  }
  
  
  ##サンプリング結果を保存
  if(rp%%keep==0){
    print(rp)
    mkeep <- rp/keep
    Util[, , mkeep] <- old.util
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    FA.A[mkeep, ] <- as.numeric(A)
    FA.D[mkeep, ] <- diag(D)
    FA.F[, , mkeep] <- Fi
    print(round(oldcov, 2))
    print(round(oldbeta[1:member-1], 2))
    print(round(oldbeta[member:length(oldbeta)], 2))
    print(round(t(A), 2))
  }
}

lovelive <- c("穂乃果", "ことり", "海未", "凛", "花陽", "真姫", "にこ", "希", "絵里", "千歌", "曜", "梨子", 
              "花丸", "善子", "aqours", "モブ")

####推定結果の要約と適合度の確認####
burnin <- 8000/keep   #バーンイン期間

##サンプリング結果を可視化
#回帰係数のプロット
matplot(BETA[, 1:4], type="l", main="メンバーごとの人気のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 5:9], type="l", main="メンバーごとの人気のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 10:11], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 12:15], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 16:20], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 21:24], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 25:29], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")

#分散供分散行列の可視化
matplot(SIGMA[, 1:4], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 5:9], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 10:13], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 14:18], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 19:22], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 23:27], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 28:31], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 32:36], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 37:40], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 41:45], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 46:49], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 50:54], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 55:58], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 59:63], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 64:67], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 68:72], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 73:76], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 77:81], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")


##推定値の事後平均の比較
#betaの要約統計量
#メンバーの切片の要約推定量
round(colMeans(BETA[burnin:, 1:(member-1)]), 3)   #betaの事後平均
round(apply(BETA[burnin:nrow(BETA), 1:(member-1)], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(BETA[burnin:nrow(BETA), 1:(member-1)], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(BETA[burnin:nrow(BETA), 1:(member-1)], 2, sd), 2)   #事後標準偏差

#説明変数の要約推定量
round(colMeans(BETA[burnin:nrow(BETA), member:ncol(BETA)]), 3)   #betaの事後平均
round(apply(BETA[burnin:nrow(BETA), member:ncol(BETA)], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(BETA[burnin:nrow(BETA), member:ncol(BETA)], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(BETA[burnin:nrow(BETA), member:ncol(BETA)], 2, sd), 2)   #事後標準偏差


#sigmaの要約統計量
round(matrix(colMeans(SIGMA[burnin:nrow(SIGMA), ]), nrow=member-1, ncol=member-1), 2)   #sigmaの事後平均
round(matrix(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.05)), nrow=member-1, ncol=member-1), 2)   #5％分位点
round(matrix(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.95)), nrow=member-1, ncol=member-1), 2)   #95％分位点
round(matrix(apply(SIGMA[burnin:nrow(SIGMA), ], 2, sd), nrow=member-1, ncol=member-1), 2) #事後標準偏差


