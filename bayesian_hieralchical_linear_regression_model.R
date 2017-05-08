#####階層ベイズ線形回帰モデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(nlme)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(13294)
#データの設定
hh <- 300   #サンプル人数
cnt <- 60   #ガチャ機会の最大値
#1人あたりのガチャ回数
for(br in 1:1000){
  pt <- c()
  for(i in 1:hh){
    p.cnt <- runif(1, 0.25, 0.9)
    t <- rbinom(1, cnt, p.cnt)
    pt <- c(pt, t)
  }
  if(min(pt)>2) break
  print(br)
}
table(pt)
hist(pt, col="grey")

hhpt <- sum(pt)   #総サンプル数
coefi <- 14   #個体内回帰係数の個数
coefh <- 9   #個体間回帰係数の個数

##データの発生
##個体間モデルの説明変数の発生
#連続変数
h.cont <- 5
Xh.cont <- matrix(runif(hh*h.cont, 0, 1), nrow=hh, ncol=h.cont)

#二値変数
h.bin <- 4
Xh.bin <- matrix(0, nrow=hh, ncol=h.bin)
for(i in 1:h.bin){
  runi <- runif(1, 0.3, 0.7)
  Xh.bin[, i] <- rbinom(hh, 1, runi)
}

#データの結合
Xh <- data.frame(hc=Xh.cont, hb=Xh.bin)


##個体内モデルの説明変数の発生
id <- rep(1:hh, rep(cnt, hh))
t <- rep(1:cnt, hh)
no. <- 1:hh*cnt
ID <- data.frame(id, t, no.)

#連続変数
cont <- 3  
x.cont <- matrix(runif(cnt*cont, 0, 1), nrow=cnt, cont) 
X.cont <- matrix(x.cont, nrow=hh*cnt, ncol=cont, byrow=T)

#二値変数
bin <- 3
x.bin <- matrix(0, nrow=cnt, ncol=bin)
for(i in 1:bin){
  runi <- runif(1, 0.3, 0.7)
  x.bin[, i] <- rbinom(cnt, 1, runi)
}
X.bin <- matrix(x.bin, nrow=hh*cnt, ncol=bin, byrow=T)

#多値変数(誰のガチャ回だったか？)
m <- 9
p <- rep(1/m, m)
r <- 10000
for(i in 1:r){
  x.multi <- t(rmultinom(cnt, 2, p))
  if(max(x.multi)==1 & min(colSums(x.multi))!=0) break
  print(i)
}
x.multiw <- x.multi[, -which.min(colSums(x.multi))]
X.multi <- matrix(t(x.multiw), nrow=hh*cnt, ncol=(m-1), byrow=T)

#ガチャした回を特定して抽出
X.k <- data.frame(c=X.cont, b=X.bin, m=X.multi)   #データの結合

index <- c()
for(i in 1:hh){
  irand <- sort(sample(1:cnt, pt[i]))   #IDごとにガチャした回数分ランダムに全ガチャ機会より抽出
  rt <- (i-1)*cnt+irand
  index <- c(index, rt)
}

ID.d <- ID[index, -3]   #ガチャした回のみ抽出
ID <- data.frame(ID.d, no=1:nrow(ID.d))   #IDを再構成
rownames(ID) <- 1:nrow(ID)   #行番号を付け直す

X <- X.k[index, ]   #個体内説明変数を再構成
rownames(X) <- 1:nrow(X)   #行番号を付け直す

Xc <- cbind(ID, X)   #IDと個体内説明変数を結合

##回帰係数と応答変数の発生
#個体間回帰モデルの回帰係数
for(rp in 1:1000){
  thetah1 <- matrix(c(runif(7, 0.2, 0.3), runif(7*h.cont, -0.2, 0.3), runif(7*h.bin, -0.3, 0.25)),
                    nrow=(1+coefh), ncol=7, byrow=T) 
  thetah2 <- matrix(c(runif(8, -0.1, 0.65), runif(8*h.cont, -0.2, 0.25), runif(8*h.bin, -0.3, 0.3)),
                    nrow=(1+coefh), ncol=8, byrow=T) 
  THETA.h <- cbind(thetah1, thetah2)
  round(THETA.h, 3)
  
  #個体内回帰モデルの回帰係数
  BETAt.h <- cbind(1, as.matrix(Xh)) %*% THETA.h + matrix(rnorm((coefi+1)*hh, 0, 0.1), hh, (coefi+1))
  
  #応答変数(ガチャ回数)の発生
  Y <- c()
  for(i in 1:hh){
    y <- BETAt.h[i, 1] + as.matrix(X[ID$id==i, ]) %*% BETAt.h[i, -1] + rnorm(pt[i], 0, 0.1)
    y <- round(exp(y), 1)
    Y <- c(Y, y)
  }
  if(max(Y) <= 200 & max(Y) >= 100) break
  print(max(Y))
}

#応答変数の要約
summary(Y)
hist(Y, col="grey", breaks=30)   #結果をプロット
data.frame(freq=names(table(Y)), y=as.numeric(table(Y)))

round(BETAt.h, 2)   #個人別の回帰係数
round(Xcomp <- cbind(ID, Y, Ylog=log(Y), X), 1)   #すべてのデータを結合


####マルコフ連鎖モンテカルロ法で階層回帰モデルを推定####
##アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
Xi <- as.matrix(cbind(1, X))
Zi <- as.matrix(cbind(1, Xh))

#事前分布の設定
sigmapr <- 0.01*diag(coefi+1)   #betaの標準偏差の事前分布
nu <- 0.01  
s0 <- 0.01
thetapr <- matrix(rep(0, (coefh+1)*(coefi+1)), nrow=(coefh+1), ncol=(coefi+1))   #階層モデルの回帰係数の事前分布
deltapr <- 0.01*diag(coefh+1)   #階層モデルの分散の事前分布
df <- 15   #逆ウィシャート分布の自由度
V <- df*diag(ncol(X)+1)   #逆ウィシャート分布のパラメータ

#サンプリング結果の保存用
oldbetas <- matrix(0, hh, (coefi+1))
BETA <- array(0, dim=c(hh, (coefi+1), R/keep))
SIGMA <- c()
THETA <- matrix(0, R/keep, (ncol(X)+1)*ncol(Zi))
VSIG <- list()

#初期値の設定
thetapr <- matrix(runif((coefh+1)*(coefi+1), -0.3, 0.5), nrow=(coefh+1), ncol=(coefi+1)) 
betapr <- Zi %*% thetapr   #betaの事前推定量
oldsigma <- rep(0.1, hh)


####ギブスサンプリングで階層回帰モデルで推定値をサンプリング####
##個体内回帰係数と標準偏差をサンプリング
for(rp in 1:R){
  sigmasq.s <- c()
  for(i in 1:hh){
    #1人分のサンプルと回帰係数を取得
    Xind <- Xi[ID$id==i, ]
    yind <- log(Y[ID$id==i])
    betapr_h <- betapr[i, ]
    
    #回帰係数の事後分布のサンプリング
    sigma_m <- solve(oldsigma[i] * t(Xind) %*% Xind + solve(sigmapr))
    beta_m <- sigma_m %*% (oldsigma[i] * t(Xind) %*% yind + solve(sigmapr) %*% betapr_h)
    sigma_diag <- diag(diag(sigma_m))   #分散共分散行列を対角化
    betan <- mvrnorm(n=1, beta_m, sigma_diag)   #多変量正規乱数からbetaをサンプリング
    
    #分散の事後分布のサンプリング
    nu1 <- nu+length(yind)
    s1 <- s0 + t(yind - Xind %*% betan) %*% (yind - Xind %*% betan)
    sigma_sq <- 1/(rgamma(1, nu1/2, s1/2))   #逆ガンマ分布からsigma^2をサンプリング
    
    oldbetas[i, ] <- betan
    sigmasq.s <- c(sigmasq.s, sigma_sq)
  }
  oldsigma <- c()
  oldsigma <- sigmasq.s
  
  ##多変量回帰モデルによる階層モデルのギブスサンプリング
  out <- rmultireg(Y=oldbetas, X=Zi, Bbar=thetapr, A=deltapr, nu=df, V=V)
  thetapr <- out$B
  sigmapr <- out$Sigma
  sigmapri <- solve(sigmapr)
  betapr <- Zi %*% thetapr
  print(rp)
  
  ##サンプリング結果を保存
  mkeep <- rp/keep
  if(rp%%keep==0){
    THETA[mkeep, ] <- as.vector(thetapr)
    VSIG[[mkeep]] <- sigmapr
    BETA[, , mkeep] <- oldbetas
    VSIG[[mkeep]] <- sigmasq.s
    #print(round(THETA[mkeep, 1:20], 2))
  }
}

####サンプリング結果の確認と適合度の確認####
burnin <- 2000   #バーンイン期間(8000サンプルまで)
RS <- R/keep 

#サンプリングされたパラメータをプロット
matplot(THETA[1:RS, 1:3], type="l", ylab="parameter")
matplot(t(BETA[1, 1:3, 1:RS]), type="l", ylab="parameter")

##個人別のパラメータ
i <- 155; sum(ID$id==i)   #個人idを抽出
round(rowMeans(BETA[i, , burnin:RS]), 3)   #個人別のパラメータ推定値の事後平均
round(BETAt.h[i, ], 3)   #個人別の真のパラメータの値
apply(BETA[i, , burnin:RS], 1, summary)   #個人別のパラメータ推定値の要約統計
apply(BETA[i, , burnin:RS], 1, function(x) quantile(x, c(0.05, 0.95)))   #個人別のパラメータ推定値の事後信用区間

hist(BETA[i, 10, burnin:RS], col="grey", xlab="beta", main="betaの個人内の事後分布", breaks=20)
hist(BETA[, 10, 5000], col="grey", xlab="beta", main="betaの個人別の事後分布", breaks=20)

##階層モデルのパラメータ
round(colMeans(THETA[burnin:RS, ]), 2)   #階層モデルのパラメータ推定値   
round(as.vector(THETA.h), 2)   #階層モデルの真のパラメータの値

##事後予測分布でガチャ回数を予測
Y.pre <- exp(BETA[i, 1, burnin:RS] + t(as.matrix(X[ID$id==i, ]) %*% BETA[i, 2:dim(BETA)[2], burnin:RS]))
index.c <- as.numeric(colnames(Y.pre)[1])
summary(Y.pre)
apply(Y.pre, 2, function(x) round(quantile(x, c(0.05, 0.95)), 3))
hist(Y.pre[, 1], col="grey", xlab="予測値", main="個人別の事後予測分布", breaks=25)
Y[index.c]   #真のガチャ回数