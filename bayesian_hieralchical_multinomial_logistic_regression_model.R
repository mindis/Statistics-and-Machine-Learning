#####階層ベイズ多項ロジットモデル#####
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

####データの発生####
#set.seed(8437)
##データの設定
hh <- 500   #サンプル数
pt <- rpois(hh, 15); pt <- ifelse(pt==0, 1, pt)   #購買機会(購買機会数が0なら1に置き換え)
hhpt <- sum(pt)
choise <- 5   #選択可能数
st <- 5   #基準ブランド
k <- 5   #説明変数の数
c <- 4   #条件付き説明変数の個数
m <- 2   #多項型説明変数の個数

##個体内説明変数の発生
#IDの設定
id <- rep(1:hh, pt)
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}
ID <- data.frame(no=1:hhpt, id, t)

#通常価格の発生
PRICE <- matrix(runif(hhpt*choise, 0.6, 1), nrow=hhpt, ncol=choise, byrow=T)   

#ディスカウント率の発生
DISC <- matrix(runif(hhpt*choise, 0, 0.5), nrow=hhpt, ncol=choise, byrow=T)

#特別陳列の発生
DISP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hhpt, 1, r)
}

#特別キャンペーンの発生
CAMP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hhpt, 1, r)
}

#カテゴリーロイヤルティ
ROYL <- matrix(0, nrow=hhpt, ncol=1)
for(i in 1:hh){
  ROYL[ID[, 2]==i] <- runif(1, 0, 1)
}

##個体内説明変数をベクトル形式に変換
#IDの設定
id.v <- rep(1:hh, pt*choise)
brand <- rep(1:choise, hhpt)
time.v <- c()
for(i in 1:hh){
  time.v <- c(time.v, rep(1:pt[i], rep(choise, pt[i])))
}
ID.v <- data.frame(id=id.v, t=time.v, b=brand)   #データの結合

#ブランド力の説明変数の設定
bv <- c(1, rep(0, choise))
bp <- matrix(bv, nrow=hhpt*(choise+1), ncol=choise, byrow=T)
BP <- subset(bp, rowSums(bp) > 0)
BP <- BP[, -choise]

#カテゴリロイヤルティの設定
index.royl <- rep(1:hhpt, rep(choise, hhpt))
ROYL.v <- matrix(0, nrow=hhpt*choise, ncol=choise)
for(i in 1:hhpt){
  ROYL.v[index.royl==i, ] <- diag(ROYL[i], choise)
}
ROYL.v <- ROYL.v[, -choise]

#その他の説明変数をベクトル化
PRICE.v <- as.numeric(t(PRICE))
DISC.v <- as.numeric(t(DISC))
DISP.v <- as.numeric(t(DISP))
CAMP.v <- as.numeric(t(CAMP))

#説明変数の結合
X <- data.frame(b=BP, PRICE=PRICE.v, DISC=DISC.v, DISP=DISP.v, CAMP=CAMP.v, ROY=ROYL.v)
XM <- as.matrix(X)


##個体間説明変数の発生
#連続変数の発生
cont <- 4
X.cont <- matrix(runif(hh*cont, 0, 1), nrow=hh, ncol=cont)

#二値変数の発生
bin <- 3
X.bin <- matrix(0, nrow=hh, ncol=bin)
for(i in 1:bin){
  X.bin[, i] <- rbinom(hh, 1, runif(1, 0.3, 0.7))
}

#多値変数の発生
multi <- 4
p.multi <- runif(multi)
X.multi <- t(rmultinom(hh, 1, p.multi))
X.multi <- X.multi[, -which.min(colSums(X.multi))]

#データの結合
XH <- data.frame(cont=X.cont, bin=X.bin, multi=X.multi)
XHi <- as.matrix(data.frame(i=1, XH))


##パラメータの設定
##個体間回帰係数の設定
#妥当な反応変数が出来るまで回帰係数の設定を繰り返す
for(t in 1:1000){
  print(t)
  len <- c + m*(choise-1)
  theta0 <- matrix(runif(len, -0.5, 2.0), nrow=1, ncol=len)   
  thetac <- matrix(runif(len*cont, -1.5, 1.5), nrow=cont, ncol=len)   
  thetab <- matrix(runif(len*bin, -1.3, 1.3), nrow=bin, ncol=len)   
  thetam <- matrix(runif(len*(multi-1), -1.4, 1.4), nrow=multi-1, ncol=len)   
  
  #回帰係数行列を作成
  THETAT <- rbind(theta0, thetac, thetab, thetam)   
  
  ##個体内回帰係数の設定
  #個体内回帰係数の誤差を決定
  Cov <- diag(runif(ncol(THETAT), 0.15, 0.4))
  er <- mvrnorm(hh, rep(0, ncol(THETAT)), Cov)
  
  #個体間回帰係数を線形結合で決定する
  BETAM <- as.matrix(XHi) %*% THETAT 
  BETA <- BETAM + er
  
  #個体内モデルの効用関数を計算
  U <- matrix(0, nrow=hhpt, ncol=choise)
  for(i in 1:hh){
    util <- XM[ID.v$id==i, ] %*% BETA[i, ]
    U.ind <- matrix(util, nrow=pt[i], ncol=choise, byrow=T)
    U[ID$id==i, ] <- U.ind
  }
  
  #効用関数から選択確率を計算して、選択ブランドを決定
  Pr <- exp(U) / rowSums(exp(U))
  Y <- t(apply(Pr, 1, function(x) t(rmultinom(1, 1, x))))
  
  #ブランド選択確率が妥当な数値ならbreak
  if(min(colMeans(Y)) > 0.05 & max(colMeans(Y) < 0.6)) {break} else {next} 
}
BETAT <- BETA

#ブランド選択結果を確認
round(colMeans(Y), 3)
colSums(Y)


####マルコフ連鎖モンテカルロ法で階層ベイズ多項混合ロジットモデルを推定####
##多項ロジットモデルの対数尤度を設定
loglike <- function(beta, Y, X, h, choise){
  #効用関数の設定
  util <- X %*% beta
  Util <- matrix(util, nrow=h, ncol=choise, byrow=T)
  
  #確率の計算
  d <- rowSums(exp(Util))
  LLl <- rowSums(Y * Util) - log(d)
  LL <- sum(LLl)
  LL.val <- list(LLl=LLl, LL=LL)
  return(LL.val)
}

#関数の最大化用対数尤度
LLfunc <- function(beta, Y, X, h, choise){
  #効用関数の設定
  util <- X %*% beta
  Util <- matrix(util, nrow=h, ncol=choise, byrow=T)
  
  #確率の計算
  d <- rowSums(exp(Util))
  LLl <- rowSums(Y * Util) - log(d)
  LL <- sum(LLl)
  return(LL)
}


##対数尤度を最大化する
b0 <- runif(ncol(XM), -1, 1)
res <- optim(b0, LLfunc, Y=Y, X=XM, h=hhpt, choise=choise, method="BFGS", hessian=FALSE, control=list(fnscale=-1))
beta_first <- res$par


##MCMCアルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4

##データの設定
#データをリスト化しておく
Y.list <- list()
X.list <- list()
for(i in 1:hh){
  Y.list[[i]] <- Y[ID$id==i, ]
  X.list[[i]] <- XM[ID.v$id==i, ]
}

#対数尤度の保存用
lognew <- matrix(0, nrow=hh, ncol=1)
logold <- matrix(0, nrow=hh, ncol=1)
logpnew <- matrix(0, nrow=hh, ncol=1)
logpold <- matrix(0, nrow=hh, ncol=1)


##事前分布の設定
Deltabar <- matrix(rep(0, ncol(XHi)*ncol(X)), nrow=ncol(XHi), ncol(X))   #階層モデルの回帰係数の平均の事前分布
Adelta <- 0.01 * diag(rep(1, ncol(XHi)))   #階層モデルの回帰係数の分散の事前分布
nu <- (ncol(X))+choise   #逆ウィシャート分布の自由度
V <- choise * diag(ncol(X))

##サンプリング結果の保存用配列
#パラメータの保存用配列
BETA <- array(0, dim=c(hh, ncol(X), R/keep))
THETA <- matrix(0, nrow=R/keep, ncol=ncol(X)*ncol(XHi))
SIGMA <- matrix(0, nrow=R/keep, ncol=ncol(X)*ncol(X))

#棄却率と対数尤度の保存用配列
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##初期値の設定
tau <- mvrnorm(hh, rep(0, ncol(X)), diag(0.3, ncol(X)))
oldbetas <- matrix(beta_first, nrow=hh, ncol=ncol(X), byrow=T) + tau
oldDelta <- solve(t(XHi) %*% XHi) %*% t(XHi) %*% oldbetas
oldVbetai <- solve(1/hh * (t(oldbetas - XHi %*% oldDelta) %*% (oldbetas - XHi %*% oldDelta)))


####MCMCで階層ベイズ多項混合ロジットモデルのパラメータをサンプリング####
for(rp in 1:R){
  
  ##MH法で個人別にbetaをサンプリング
  rw <- matrix(rnorm(length(oldbetas), 0, 0.15), nrow=nrow(oldbetas), ncol=ncol(oldbetas))
  betad <- oldbetas
  betan <- oldbetas + rw
  
  #パラメータの事前分布との誤差を計算
  er_new <- betan - XHi %*% oldDelta
  er_old <- betad - XHi %*% oldDelta
  
  
  #ID別に対数尤度と対数事前分布を計算
  for(i in 1:hh){
    lognew[i, ] <- loglike(betan[i, ], Y.list[[i]], X.list[[i]], pt[i], choise)$LL
    logold[i, ] <- loglike(betad[i, ], Y.list[[i]], X.list[[i]], pt[i], choise)$LL
    logpnew[i, ] <- -0.5 * (er_new[i, ] %*% oldVbetai %*% er_new[i, ])
    logpold[i, ] <- -0.5 * (er_old[i, ] %*% oldVbetai %*% er_old[i, ])
  }
  
  ##MHサンプリング
  #サンプリングを採択するかどうかを決定
  rand <- matrix(runif(hh*ncol(oldbetas)), nrow=hh, ncol=ncol(oldbetas))   #一様乱数から乱数を発生
  LLind.diff <- exp(lognew + logpnew - logold - logpold)   #採択率を計算
  alpha <- matrix(ifelse(LLind.diff > 1, 1, LLind.diff), nrow=hh, ncol=ncol(oldbetas))   
  
  #alphaに基づきbetaを採択
  oldbetas.r <- ifelse(alpha > rand, oldbetas.r <- betan, oldbetas.r <- betad)   #alphaがrandを上回っていたら採択
  logl <- ifelse(alpha[, 1] > rand[, 1], lognew, logold)
  
  adopt <- sum(oldbetas[, 1]==oldbetas.r[, 1])/hh   #採択率
  LLho <- sum(logl)   #対数尤度の総和
  oldbetas <- oldbetas.r   #パラメータを更新
  
  
  ##多変量回帰モデルによる階層モデルのサンプリング
  out <- rmultireg(Y=oldbetas, X=XHi, Bbar=Deltabar, A=Adelta, nu=nu, V=V)
  oldDelta <- out$B
  sig <- out$Sigma
  oldVbetai <- solve(sig)
  
  
  ##サンプリング結果を保存
  if(rp%%keep==0){
    print(rp)
    mkeep <- rp/keep
    BETA[, , mkeep] <- oldbetas
    THETA[mkeep, ] <- as.numeric(oldDelta)
    SIGMA[mkeep, ] <- as.numeric(sig)
    llike[mkeep] <- LLho
    print(adopt)
    print(LLho)
    print(round(cbind(oldDelta, THETAT), 1))
  }
}


round(oldDelta - THETAT, 2)
round(oldDelta, 3)
round(THETAT, 3)



