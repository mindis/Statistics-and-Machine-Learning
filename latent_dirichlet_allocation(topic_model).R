#####トピックモデル(Latent Dirichlet Allocation Model)#####
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
w <- 200   #1文書あたりの単語数 

#パラメータの設定
alpha0 <- runif(k, 0.1, 1.25)   #文書のディレクリ事前分布のパラメータ
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
R <- 5000   #サンプリング回数
keep <- 2   #4回に1回の割合でサンプリング結果を利用
iter <- 0

#ハイパーパラメータの事前分布の設定
alpha <- rep(1.5, k)
beta <- rep(100/v, v)

#パラメータの初期値
theta.ini <- runif(k, 0.5, 2)
phi.ini <- runif(v, 0.5, 1)
theta <-    rdirichlet(d, theta.ini)   #文書トピックのパラメータの初期値
phi <- rdirichlet(k, phi.ini)   #単語トピックのパラメータの初期値


#パラメータの格納用配列
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
W.SEG <- matrix(0, nrow=d*w, ncol=R/(keep*2))


####ギブスサンプリングでトピックモデルを推定####
id1 <- rep(1:d, v)
id2 <- rep(1:v, rep(d, v))

##単語頻度行列をベクトルに変換
d.vec <- c() 
v.vec <- c()
for(i in 1:v){
  for(j in 1:d){
    d.vec <- c(d.vec, rep(j, WX[j, i]))
    v.vec <- c(v.vec, rep(i, WX[j, i]))
  }
}

index.vec <- cbind(v.vec, d.vec)
id.w <- cbind(id2, id1)

#抽出する潜在確率の行を取得
index.z <- apply(index.vec, 1, function(x) subset(1:nrow(id.w), id.w[, 1]==x[1] & id.w[, 2]==x[2]))

#トピックzをサンプリングするための配列
theta.seg <- array(0, dim=c(d, v, k))
phi.seg <- array(0, dim=c(d, v, k))
burden.seg <- array(0, dim=c(d, v, k))
zpt <- array(0, dim=c(d, v, k))

#パラメータをサンプリングするための配列
dir.theta <- matrix(0, nrow=d, ncol=k) 
dir.phi <- matrix(0, nrow=k, ncol=v)

####ギブスサンプリング####
##単語のトピックzをサンプリング
for(rp in 1:R){
  
  burden.sum <- matrix(0, nrow=d, ncol=v) 
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
  ZW <- ZP[index.z, ]

  #単語ごとにトピックzをサンプリング
  Z <- t(apply(ZW, 1, function(x) rmultinom(1, 1, x)))
  
  ##文書トピック分布thetaをサンプリング
  #ディレクリ分布のパラメータを決定
  for(i in 1:d){
   dir.theta[i, ] <- colSums(Z[d.vec==i, ]) + alpha 
  }
  
  theta <- t(apply(dir.theta, 1, function(x) rdirichlet(1, x)))   #ディレクリ分布からthetaをサンプリング

  ##単語トピック分布phiをサンプリング
  #ディレクリ分布のパラメータを決定
  for(i in 1:v){
   dir.phi[, i] <- colSums(Z[v.vec==i, ]) + beta[i]
  }
  
  phi <- t(apply(dir.phi, 1, function(x) rdirichlet(1, x)))   #ディレクリ分布からphiをサンプリング
  
  ##ハイパーパラメータの更新
  
  
  ##サンプリング結果の保存
  #サンプリングを保存する回数ならbetaを書き込む
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
  }
  if(rp%%(keep*2)==0){
    W.SEG[, mkeep/2] <- Z %*% 1:k
    print(rp)
  }
}


####推定結果と適合度####
burnin <- 1500
Rkeep <- R/keep

##サンプリング結果のプロット
matplot(t(THETA[1:3, 1, ]), type="l", lty=1, ylab="value")
matplot(t(THETA[4:6, 1, ]), type="l", lty=1, ylab="value")
matplot(t(PHI[1, 1:3, ]), type="l", lty=1, ylab="value")
matplot(t(PHI[2, 1:3, ]), type="l", lty=1, ylab="value")

#thetaの事後平均と真のthetaの比較
round(THETA1 <- cbind(theta0, rowMeans(THETA[, 1, burnin:Rkeep]), rowMeans(THETA[, 2, burnin:Rkeep]), 
      rowMeans(THETA[, 3, burnin:Rkeep]), rowMeans(THETA[, 4, burnin:Rkeep]), 
      rowMeans(THETA[, 5, burnin:Rkeep])), 3)

#phiの事後平均と真のphiの比較
round(PHI1 <- rbind(phi0, rowMeans(PHI[1, , burnin:Rkeep]), rowMeans(PHI[2, , burnin:Rkeep]),
      rowMeans(PHI[3, , burnin:Rkeep]), rowMeans(PHI[4, , burnin:Rkeep]), 
      rowMeans(PHI[5, , burnin:Rkeep])), 3)

#トピックの事後分布
segment <- apply(W.SEG[, (burnin/2):(Rkeep/2)], 1, table)


