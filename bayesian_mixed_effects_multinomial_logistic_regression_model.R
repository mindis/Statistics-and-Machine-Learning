#####変量効果混合ロジットモデル#####
library(MASS)
library(mlogit)
library(MCMCpack)
library(bayesm)
library(caret)
library(reshape2)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)


####データの発生####
#set.seed(8437)
##データの設定
hh <- 250   #サンプル数
pt <- rpois(hh, 15); pt <- ifelse(pt==0, 1, pt)   #購買機会(購買機会数が0なら1に置き換え)
hhpt <- sum(pt)
choise <- 5   #選択可能数
st <- 5   #基準ブランド
k <- 5   #説明変数の数

##説明変数の発生
#IDの設定
id <- rep(1:hh, pt)
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}
ID <- cbind(no=1:hhpt, id, t)

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
ROYL <- matrix(runif(hhpt, 0, 1), nrow=hhpt, ncol=1)

##パラメータの設定
beta1 <- -5.8   #価格のパラメータ
beta2 <- 5.5   #割引率のパラメータ
beta3 <- 2.0   #特別陳列のパラメータ
beta4 <- 1.8   #キャンペーンのパラメータ
b1 <- c(1.1, 0.6, -0.7, -0.3)   #カテゴリーロイヤルティのパラメータ
b0 <- c(0.5, 0.9, 1.4, 2.1)   #ブランド1～4の相対ベース販売力
betat <- c(b0, beta1, beta2, beta3, beta4)

##変量効果の設定
k.random <- 5   #変量効果の変数数
b0.random <- matrix(b0, nrow=hh, ncol=choise-1, byrow=T) + 
  matrix(rnorm(hh*(choise-1), 0, 1), nrow=hh, ncol=choise-1, byrow=T)
beta1.random <- matrix(beta1, nrow=hh, ncol=1, byrow=T) + 
  matrix(rnorm(hh, 0, sqrt(3.5)), nrow=hh, ncol=1, byrow=T)


##効用を発生させ、選択されたブランドを決定
#ロジットの発生
logit <- matrix(0, nrow=hhpt, ncol=st)
for(i in 1:hh){
  r <- subset(1:hhpt, ID[, 2]==i)
  for(j in 1:(st-1)){
    logit[r, j] <- b0.random[i, j] + beta1.random[i]*PRICE[r, j]  + beta2*DISC[r, j] + beta3*DISP[r, j] + beta4*CAMP[r, j] 
  }
}

#基準変数のロジットを計算
for(i in 1:hh){
  r <- subset(1:hhpt, ID[, 2]==i)
  logit[r, st] <- beta1.random[i]*PRICE[r, st] + beta2*DISC[r, st] + beta3*DISP[r, st] + beta4*CAMP[r, st] 
}

##発生させたロジットから選択ブランドを決定
#ブランド選択確率を計算
Pr <- exp(logit)/rowSums(exp(logit))
colMeans(Pr); apply(Pr, 2, summary)

#選択ブランドを発生
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
colMeans(Y); apply(Y, 2, table)

round(cbind(Y %*% 1:choise, Pr), 3)   #選択結果と選択確率

####マルコフ連鎖モンテカルロ法で変量効果ロジスティック回帰モデルを推定####
##回帰モデルを推定するために説明変数をベクトル形式に変更設定
#idを設定
id.v <- c()
for(i in 1:hh){
  id.v <- c(id.v, rep(ID[ID[, 2]==i, 2], choise))
}

#切片の設定
p <- c(1, rep(0, choise))
bp <- matrix(p, nrow=hhpt*(choise+1), ncol=choise, byrow=T)
BP <- subset(bp, rowSums(bp) > 0)
BP <- BP[, -st] 

#カテゴリロイヤルティの設定
index.royl <- rep(1:hhpt, rep(choise, hhpt))
ROYL.v <- matrix(0, nrow=hhpt*choise, ncol=choise)

for(i in 1:hhpt){
  ROYL.v[index.royl==i, ] <- diag(c(rep(ROYL[i, ], choise-1), 0))
}
ROYL.v <- ROYL.v[, -st]

#説明変数の設定
PRICE.v <- as.numeric(t(PRICE))
DISC.v <- as.numeric(t(DISC))
DISP.v <- as.numeric(t(DISP))
CAMP.v <- as.numeric(t(CAMP))

round(X <- data.frame(b=BP, PRICE=PRICE.v, DISC=DISC.v, DISP=DISP.v, CAMP=CAMP.v), 2)   #データの結合
XM <- as.matrix(X)


#Zの設定
Z <- matrix(0, nrow=nrow(XM), ncol=k.random*hh)
for(i in 1:hh){
  r <- subset(1:nrow(XM), id.v==i)
  c <- ((i-1)*k.random+1):((i-1)*k.random+k.random)
  Z[r, c] <- XM[id.v==i, 1:k.random]
}


##MCMCアルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- c()   #対数尤度の保存用

##事前分布の設定
#固定効果の事前分布
betas.fix <- rep(0, ncol(XM))  #回帰係数の平均の事前分布
sigma.fix <- diag(rep(0.01, ncol(XM)))   #回帰係数の事前分布の分散

#変量効果の事前分布
Deltabar <- rep(0, hh*choise)
Adelta <- 0.01*diag(k.random)
nu <- k.random   #逆ウィシャート分布の自由度
V <- nu * diag(rep(1, k.random))
beta.random <- matrix(0, nrow=hh, ncol=k.random)   #変量効果の事前分布の平均を0に固定


##サンプリング結果の保存用
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
B.random <- array(0, dim=c(hh, k.random, R/keep))
SIGMA <- matrix(0, nrow=R/keep, ncol=k.random^2)

##初期値の設定
#回帰係数の初期値
oldbeta.f <- c(c(runif(choise-1, 0, 3)), -5.0, 5.0, runif(2, 1, 2))  

#変量効果の初期値
cov.random <- diag(runif(k.random, 0.5, 2))
oldbeta.r <- mvrnorm(hh, c(rep(0, k.random-1), -5), cov.random)   #変量効果の初期値


#階層モデルの初期値
beta.random <- matrix(0, nrow=hh, ncol=k.random)

####マルコフ連鎖モンテカルロ法で推定####
##mixed logitモデルの対数尤度
LLike <- function(beta, b, X, Z, Y, hh, choise){
  d <- rowSums(exp(matrix(X %*% beta + Z %*% b, nrow=hh, ncol=choise, byrow=T)))
  LLl <- rowSums(Y * matrix(X %*% beta + Z %*% b, nrow=hh, ncol=choise, byrow=T)) - log(d)
  LL <- sum(LLl)
  LL.val<- list(LLl=LLl, LL=LL)
  return(LL.val)
}

##マルコフ連鎖モンテカルロ法でパラメータをサンプリング
for(rp in 1:R) {
  ##MHサンプリングで固定効果betaのサンプリング
  oldbeta.rv <- as.numeric(t(oldbeta.r))
  betad.f <- oldbeta.f
  betan.f <- betad.f + rnorm(length(betad.f), 0, 0.025)
  
  #ランダムウォークサンプリング
  #対数尤度と対数事前分布を計算
  lognew.f <- LLike(beta=betan.f, b=oldbeta.rv, X=XM, Z=Z, Y=Y, hh=hhpt, choise=choise)$LL
  logold.f <- LLike(beta=betad.f, b=oldbeta.rv, X=XM, Z=Z, Y=Y, hh=hhpt, choise=choise)$LL
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
  rw <- t(0.5 * chol(cov.random) %*% t(matrix(rnorm(hh*(k.random)), nrow=hh, ncol=k.random)))
  
  betan.random <- betad.random + rw
  betad.r <- as.numeric(t(betad.random))
  betan.r <- as.numeric(t(betan.random))
  inv.cov <- solve(cov.random)
  
  #事前分布の誤差を計算
  er.new <- betan.random - beta.random
  er.old <- betad.random - beta.random
  
  #対数尤度と対数事前分布を計算
  lognew.r <- LLike(beta=oldbeta.f, b=betan.r, X=XM, Z=Z, Y=Y, hh=hhpt, choise=choise)$LLl
  logold.r <- LLike(beta=oldbeta.f, b=betad.r, X=XM, Z=Z, Y=Y, hh=hhpt, choise=choise)$LLl
  logpnew.r <- apply(er.new, 1, function(x) -0.5 * x %*% inv.cov %*% x)
  logpold.r <- apply(er.old, 1, function(x) -0.5 * x %*% inv.cov %*% x)
  
  #ID別に対数尤度の和を取る
  lognew.rind <- data.frame(logl=lognew.r, id=ID[, 2]) %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(sum=sum(logl))
  
  logold.rind <- data.frame(logl=logold.r, id=ID[, 2]) %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(sum=sum(logl))
  
  #MHサンプリング
  rand <- matrix(runif(hh), nrow=hh, ncol=k.random)
  LLind.diff <- exp(lognew.rind$sum + logpnew.r - logold.rind$sum - logpold.r)   #棄却率を計算
  alpha <- matrix(ifelse(LLind.diff > 1, 1, LLind.diff), nrow=hh, ncol=k.random)      
  
  oldbeta.r <- ifelse(alpha > rand, oldbeta.r <- betan.random, oldbeta.r <- betad.random)   #alphaがrandを上回っていたら採択
  logl <- ifelse(alpha[, 1] > rand[, 1], logl <- lognew.r, logl <- logold.r)
  
  
  ##逆ウィシャート分布からsigmaをサンプリング
  #逆ウィシャート分布のパラメータ
  V <- var(oldbeta.r)
  VK <- k.random * diag(k.random) + hh * V
  nu1 <- hh + nu - 1 
  
  #逆ウィシャート分布から分散共分散行列を発生
  cov.random <- rwishart(nu1, solve(VK))$IW   
  
  
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta.f
    B.random[, , mkeep] <- oldbeta.r
    SIGMA[mkeep, ] <- as.numeric(cov.random)
    
    print(sum(logl))
    print(rp)
    print(round(mean(alpha), 3)); print(round(alpha.f, 3))
    print(round(rbind(oldbeta.f, betat), 3))
    print(round(cov.random, 3))
  }
}

####推定結果と要約####
burnin <- 2500
i <- 6

matplot(BETA[, 1:4], type="l")
matplot(BETA[, 5:8], type="l")
matplot(SIGMA[, c(1, 7, 13)], type="l")
matplot(SIGMA[, c(19, 25)], type="l")
matplot(t(B.random[i, 1:2, ]), type="l")
matplot(t(B.random[i, 3:4, ]), type="l")
plot(1:5000, t(B.random[i, 5, ]), type="l")

#変量効果の分散の事後平均
colMeans(SIGMA[burnin:nrow(SIGMA), c(1, 7, 13, 19, 25)]) 
colMeans(t(B.random[i, , burnin:nrow(SIGMA)])) 
c(b0.random[i, ]-b0, beta1.random[i, ]-beta1)
