#####トピックモデル#####
library(MASS)
library(lda)
library(RMeCab)
library(gtools)
library(reshape2)
library(plyr)
library(ggplot2)


####データの発生####
#set.seed(423943)
#データの設定
k <- 5   #トピック数
d <- 200   #文書数
v <- 100   #語彙数
w <- 300   #1文書あたりの単語数

#パラメータの設定
alpha0 <- round(runif(k, 1.5, 2.5), 3)   #文書のディレクリ事前分布のパラメータ
alpha1 <- rep(0.5, v)   #単語のディレクリ事前分布のパラメータ

#ディレクリ乱数の発生
theta0 <- rdirichlet(d, alpha0)   #文書のトピック分布をディレクリ乱数から発生
phi0 <- rdirichlet(k, alpha1)   #単語のトピック分布をディレクリ乱数から発生

#多項分布の乱数からデータを発生
WX <- matrix(0, d, v)
ZS <- list()
for(i in 1:d){
  z <- t(rmultinom(w, 1, theta0[i, ]))   #文書のトピック分布を発生
  
  zn <- z %*% c(1:k)   #0,1を数値に置き換える
  zdn <- cbind(zn, z)   #apply関数で使えるように行列にしておく
  
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi0[x[1], ])))   #文書のトピックから単語を生成
  wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
  WX[i, ] <- wdn  
  ZS[[i]] <- zdn[, 1]
  print(i)
}
barplot(colSums(z), names.arg=c("seg1", "seg2", "seg3", "seg4", "seg5"))
round(colSums(WX)/sum(WX), 3)   #単語の出現頻度


####マルコフ連鎖モンテカルロ法でトピックモデルを推定####
##マルコフ連鎖モンテカルロ法の設定
#アルゴリズムの設定
R <- 20000   #サンプリング回数
keep <- 4   #4回に1回の割合でサンプリング結果を利用
iter <- 0

#ハイパーパラメータの事前分布の設定
alpha <- rep(1.5, k)
beta <- rep(100/v, v)

#パラメータの初期値
theta.ini <- matrix(runif(d*k), nrow=d, ncol=k)
phi.ini <- matrix(runif(v*k), nrow=k, ncol=v)
theta <- theta.ini/rowSums(theta.ini)   #文書トピックのパラメータの初期値
phi <- phi.ini/rowSums(phi.ini)   #単語トピックのパラメータの初期値

#パラメータの格納用配列
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
W.SEG <- matrix(0, nrow=d*v, ncol=R/keep)


####ギブスサンプリングでトピックモデルを推定####
id1 <- rep(1:d, v)
id2 <- rep(1:v, rep(d, v))

##単語のトピックzをサンプリング
for(rp in 1:R){
  
  #トピックzをサンプリングするための配列
  theta.seg <- array(0, dim=c(d, v, k))
  phi.seg <- array(0, dim=c(d, v, k))
  burden.seg <- array(0, dim=c(d, v, k))
  burden.sum <- matrix(0, nrow=d, ncol=v) 
  zpt <- array(0, dim=c(d, v, k))
  
  #潜在確率のパラメータを計算
  for(i in 1:k){
    theta.seg[, , i] <- matrix(theta[, i], nrow=d, ncol=v)
    phi.seg[, , i] <- matrix(phi[i, ], nrow=d, ncol=v, byrow=T)
    burden.seg[, , i] <- theta.seg[, , i] * phi.seg[, , i]
    burden.sum <- burden.sum + burden.seg[, , i]
  }

  
  #潜在確率を計算
  for(i in 1:k){
    zpt[, , i] <- burden.seg[, , i] / burden.sum
  }
  ZP <- matrix(zpt, nrow=v*d, ncol=k)   #潜在確率を単語順に配列
 
  #単語ごとにトピックzをサンプリング
  Z <- t(apply(ZP, 1, function(x) rmultinom(1, 1, x)))
  
  
  ##文書トピック分布thetaをサンプリング
  #ディレクリ分布のパラメータを決定
  dir.theta <- matrix(0, nrow=d, ncol=k) 
  for(i in 1:d){
   dir.theta[i, ] <- colSums(Z[id1==i, ]) + alpha 
  }
  theta <- t(apply(dir.theta, 1, function(x) rdirichlet(1, x)))   #ディレクリ分布からthetaをサンプリング
  
  ##単語トピック分布phiをサンプリング
  #ディレクリ分布のパラメータを決定
  dir.phi <- matrix(0, nrow=k, ncol=v)
  for(i in 1:v){
   dir.phi[, i] <- colSums(matrix(WX[, i], nrow=d, ncol=k) * Z[id2==i, ]) + beta[i]
  }
  phi <- t(apply(dir.phi, 1, function(x) rdirichlet(1, x)))   #ディレクリ分布からphiをサンプリング
  
  ##サンプリング結果の保存
  #サンプリングを保存する回数ならbetaを書き込む
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    W.SEG[, mkeep] <- Z %*% 1:5
    print(rp)
  }
}

Z.cnt <- t(apply(W.SEG[, 2500:500], 1, function(x) xtabs(~ x)))

ZS.cnt
for(i in 1:v){
  for(j in 1:d){
    ZS[[j]][i]
  }
}
ZS[[1]]


colMeans(rdirichlet(1000, dir.phi[1, ]))


cbind(round(theta0, 3), round(THETA[, , 2000], 3))

round(rowMeans(THETA[1:100, 1, 1000:2000]), 3)
round(theta0[1:100, 1], 3)

plot(THETA[1, 1, 1:2000], type="l")

t(THETA[1, , 1000:2000])

rowMeans(THETA[1:10, 1, 1000:2000])
theta[1:10, 1]

colSums(matrix(WX[, i], nrow=d, ncol=k) * Z[id2==i, ])

rowMeans(THETA[1:100, 1, 2500:5000])
theta0[1:100, 3]
