####相関のある多項ロジットモデル#####
library(MASS)
library(mlogit)
library(MCMCpack)
library(bayesm)
library(caret)
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


####データの発生####
#set.seed(8437)
##データの設定
hh <- 1000   #サンプル数
choise <- 5   #選択可能数
st <- 5   #基準ブランド
k <- 5   #説明変数の数

##説明変数の発生
#通常価格の発生
PRICE <- matrix(runif(hh*choise, 0.6, 1), nrow=hh, ncol=choise)   

#ディスカウント率の発生
DISC <- matrix(runif(hh*choise, 0, 0.5), nrow=hh, ncol=choise)

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
ROYL <- matrix(runif(hh, 0, 1), nrow=hh, ncol=1)

##分散共分散行列の設定
MU <- matrix(c(1, 1, 0, 1, 1), nrow=hh, ncol=choise, byrow=T)
Cov <- matrix(c(3.5, 0, 0, 1), nrow=2, ncol=2)
RMV <- mvrnorm(hh, rep(0, 2), Cov)   #分散成分
MU <- cbind(RMV[, 1], RMV[, 1], 0, RMV[, 2], RMV[, 2])


##パラメータの設定
beta1 <- -5.8   #価格のパラメータ
beta2 <- 5.5   #割引率のパラメータ
beta3 <- 2.0   #特別陳列のパラメータ
beta4 <- 1.8   #キャンペーンのパラメータ
betat <- c(beta1, beta2, beta3, beta4)

##階層モデルの設定
b1 <- c(1.1, 0.6, -0.7, -0.3)   #カテゴリーロイヤルティのパラメータ
b0 <- c(0.5, 0.8, 1.2, 2.0)   #ブランド1〜4の相対ベース販売力
beta0 <- matrix(b0, nrow=hh, ncol=choise-1, byrow=T)

##効用を発生させ、選択されたブランドを決定
#多変量正規分布からロジットを発生
logit <- matrix(0, nrow=hh, ncol=st)
for(i in 1:(st-1)){
  logit[, i] <- beta0[, i] + beta1*PRICE[, i]  + beta2*DISC[, i] + beta3*DISP[, i] + beta4*CAMP[, i] + MU[, i]
}
#基準変数のロジットを計算
logit[, st] <- beta1*PRICE[, st]  + beta2*DISC[, st] + beta3*DISP[, st] + beta4*CAMP[, st] + MU[, st]

##発生させたロジットから選択ブランドを決定
#ブランド選択確率を計算
Pr <- exp(logit)/rowSums(exp(logit))
colMeans(Pr); apply(Pr, 2, summary)

round(cbind(Pr, Pr1), 3)

#選択ブランドを発生
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
colMeans(Y); apply(Y, 2, table)

round(cbind(Y %*% 1:choise, Pr), 3)

####マルコフ連鎖モンテカルロ法で相関構造のある混合型ロジスティック回帰モデルを推定####
##回帰モデルを推定するために説明変数をベクトル形式に変更設定
#IDの設定
ID <- rep(1:hh, rep(choise, hh))

#切片の設定
p <- c(1, rep(0, choise))
bp <- matrix(p, nrow=hh*(choise+1), ncol=choise, byrow=T)
BP <- subset(bp, rowSums(bp) > 0)
BP <- BP[, -st] 

#カテゴリロイヤルティの設定
ROYL.v <- matrix(0, nrow=hh*choise, ncol=choise)
for(i in 1:hh){
  ROYL.v[ID==i, ] <- diag(c(rep(ROYL[i, ], choise-1), 0))
}
ROYL.v <- ROYL.v[, -st]

#説明変数の設定
PRICE.v <- as.numeric(t(PRICE))
DISC.v <- as.numeric(t(DISC))
DISP.v <- as.numeric(t(DISP))
CAMP.v <- as.numeric(t(CAMP))

round(X <- data.frame(BP, PRICE=PRICE.v, DISC=DISC.v, DISP=DISP.v, CAMP=CAMP.v), 2)   #データの結合
XM <- as.matrix(X)

#IDの設定
brand <- rep(1:choise, hh)
id <- rep(1:hh, rep(choise, hh))
ID <- data.frame(id, brand)

#Zの設定
z <- cbind(c(1, 1, 0, 0, 0), c(rep(0, choise)), c(0, 0, 0, 1, 1))
Z <- kronecker(diag(hh), z)
index.z <- subset(1:ncol(Z), colSums(Z)==0)
Z <- Z[, -index.z]

##MCMCアルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 2
llike <- c()   #対数尤度の保存用

##事前分布の設定
#固定効果の事前分布
betas.fix <- rep(0, ncol(XM))  #回帰係数の平均の事前分布
sigma.fix <- diag(rep(0.01, ncol(XM)))   #回帰係数の事前分布の分散

#変量効果の事前分布
Deltabar <- rep(0, hh*choise)
Adelta <- 0.01*diag(2)
nu <- 2   #逆ウィシャート分布の自由度
V <- nu * diag(rep(1, 2))
beta.random <- matrix(0, nrow=hh, ncol=2)   #変量効果の事前分布の平均を0に固定


##サンプリング結果の保存用
Util <- array(0, dim=c(hh, choise-1, R/keep))
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
BETA0 <- matrix(0, R/keep, ncol=2)
SIGMA <- matrix(0, nrow=R/keep, ncol=2^2)

##初期値の設定
#回帰係数の初期値
oldbeta.f <- c(c(runif(choise-1, 0, 3)), -3.0, 3.0, runif(2, 0, 2))  

#変量効果の初期値
cov.random <- diag(runif(2, 0.1, 1))
oldbeta.r <- matrix(mvrnorm(hh, rep(0, 2), cov.random), nrow=hh, ncol=2, byrow=T)   #変量効果の初期値

#階層モデルの初期値
beta.random <- matrix(0, nrow=hh, ncol=2)

####マルコフ連鎖モンテカルロ法で推定####
##mixed logitモデルの対数尤度
LLike <- function(beta, b, X, Z, Y, hh, choise){
  d <- rowSums(exp(matrix(X %*% beta + Z %*% b, nrow=hh, ncol=choise, byrow=T)))
  LLl <- rowSums(Y * matrix(X %*% beta + Z %*% b, nrow=hh, ncol=choise, byrow=T)) - log(d)
  LL <- sum(LLl)
  LL.val<- list(LLl=LLl, LL=LL)
  return(LL.val)
}


for(rp in 1:R){
  ##MHサンプリングで固定効果betaのサンプリング
  oldbeta.rv <- as.numeric(t(oldbeta.r))
  betad.f <- oldbeta.f
  betan.f <- betad.f + rnorm(length(betad.f), 0, 0.05)   #ランダムウォークサンプリング
  
  #対数尤度と対数事前分布を計算
  lognew.f <- LLike(beta=betan.f, b=oldbeta.rv, X=XM, Z=Z, Y=Y, hh=hh, choise=choise)$LL
  logold.f <- LLike(beta=betad.f, b=oldbeta.rv, X=XM, Z=Z, Y=Y, hh=hh, choise=choise)$LL
  logpnew.f <- lndMvn(betan.f, betas.fix, sigma.fix)
  logpold.f <- lndMvn(betad.f, betas.fix, sigma.fix)
  
  #MHサンプリング
  alpha.f <- min(1, exp(lognew.f + logpnew.f - logold.f - logpold.f))
  if(alpha.f == "NAN") alpha.f <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しい固定効果betaを採択
  if(u < alpha.f){
    oldbeta.f <- betan.f
    logl.f <- lognew.f
    
    #そうでないなら固定効果betaを更新しない
  } else {
    logl.f <- logold.f
  }
  
  
  ##MHサンプリングで個人別に変量効果betaをサンプリング
  betad.random <- oldbeta.r 
  rw <- t(1.25 * chol(cov.random) %*% t(matrix(rnorm(hh*(2)), nrow=hh, ncol=2)))
  
  betan.random <- betad.random + rw
  betad.r <- as.numeric(t(betad.random))
  betan.r <- as.numeric(t(betan.random))
  
  #対数尤度と対数事前分布を計算
  lognew.r <- LLike(beta=oldbeta.f, b=betan.r, X=XM, Z=Z, Y=Y, hh=hh, choise=choise)$LLl
  logold.r <- LLike(beta=oldbeta.f, b=betad.r, X=XM, Z=Z, Y=Y, hh=hh, choise=choise)$LLl
  logpnew.r <- apply((betan.random - beta.random), 1, function(x) -0.5 * x %*% solve(cov.random) %*% x)
  logpold.r <- apply((betad.random - beta.random), 1, function(x) -0.5 * x %*% solve(cov.random) %*% x)
  
  #MHサンプリング
  rand <- matrix(runif(hh), nrow=hh, ncol=2)
  LLind.diff <- exp(lognew.r + logpnew.r - logold.r - logpold.r)   #棄却率を計算
  alpha <- matrix(ifelse(LLind.diff > 1, 1, LLind.diff), nrow=hh, ncol=2)      
  
  oldbeta.r <- ifelse(alpha > rand, oldbeta.r <- betan.random, oldbeta.r <- betad.random)   #alphaがrandを上回っていたら採択
  logl <- ifelse(alpha[, 1] > rand[, 1], logl <- lognew.r, logl <- logold.r)
  
 
  ##正規分布からbeta0をサンプリング
  #beta0 <- colMeans(oldbeta.r)
  #ohm <- cov.random/hh
  #beta0.mv <-  mvrnorm(1, beta0, ohm)
  #beta0.mv + t(chol(cov.random/hh))*rnorm(choise-1)
  #beta.random <- matrix(beta0.mv, nrow=hh, ncol=choise-1, byrow=T)
  
  ##逆ウィシャート分布からsigmaをサンプリング
  V <- var(oldbeta.r)
  VK <- 2 * diag(2) + hh * V
  nu1 <- hh + nu - 1 
  
  cov.random <- rwishart(nu1, solve(VK))$IW   #逆ウィシャート分布から分散共分散行列を発生

  
  ##ギブスサンプリングでDeltaをサンプリング
  #M <- matrix(c(1, 0), hh, 2, byrow=T)   #仮想的な0の説明変数を作成
  #DeltaM <- matrix(Deltabar, 2, 2, byrow=T)   #仮想的な0の回帰係数の事前分布を作成
  
  #out <- rmultireg(oldbeta.r, M, DeltaM, Adelta, nu, V)   #多変量回帰モデルのギブスサンプラー
  #cov.random <- out$Sigma
  
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta.f
    BETA0[mkeep, ] <- beta.random[1, ]
    SIGMA[mkeep, ] <- as.numeric(cov.random)
    Util[, , mkeep] <- oldbeta.r
    
    print(sum(logl))
    print(rp)
    print(round(mean(alpha), 3)); print(round(alpha.f, 3))
    print(round(rbind(oldbeta.f, c(b0, betat)), 3))
    print(round(cbind(cov.random, Cov), 3))
  }
}

burnin <- 5000

matplot(BETA[, 1:4], type="l")
matplot(BETA[, 5:8], type="l")
matplot(SIGMA[, c(1, 4)], type="l")

colMeans(SIGMA[, c(1, 4)]) - pi^2/6



