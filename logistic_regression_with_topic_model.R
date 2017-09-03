#####トピックモデルを含んだロジスティック回帰モデル#####
library(MASS)
library(lda)
library(bayesm)
library(RMeCab)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(317209)

####データの発生####
#データの設定
k <- 5   #トピック数
d <- 1000   #ユーザー数
v <- 50   #アイテム数
w <- rpois(d, 50)   #購買数

####トピックモデルからアイテム購買を生成####
##トピックモデルの設定
#パラメータの設定
alpha0 <- runif(k, 0.1, 0.8)   #ユーザーのディレクリ事前分布のパラメータ
alpha1 <- rep(0.25, v)   #アイテムのディレクリ事前分布のパラメータ

#ディレクリ乱数の発生
theta0 <- rdirichlet(d, alpha0)   #ユーザーのトピック分布をディレクリ乱数から発生
phi0 <- rdirichlet(k, alpha1)   #アイテムのトピック分布をディレクリ乱数から発生

#多項分布からデータを発生
WX <- matrix(0, nrow=d, ncol=v)
Z <- matrix(0, nrow=d, ncol=k)
Rate <- matrix(0, nrow=d, ncol=k)
ZS <- list()

for(i in 1:d){
  #トピックを生成
  z <- t(rmultinom(w[i], 1, theta0[i, ]))   #ユーザーのアイテムのトピックを生成
  Z[i, ] <- colSums(z)
  Rate[i, ] <- Z[i, ] / w[i]
  zn <- z %*% c(1:k)   #0,1を数値に置き換える
  zdn <- cbind(zn, z)   #apply関数で使えるように行列にしておく
 
  #トピックから応答変数を生成
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi0[x[1], ])))   #ユーザーのトピックからアイテムを生成
  
  wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
  WX[i, ] <- wdn
  ZS[[i]] <- cbind(rep(i, w[i]), zdn[, 1])
  print(i)
}

#リストを行列方式に変更
ZS <- do.call(rbind, ZS)

#トピックの単純集計
z_table <- table(ZS[, 2])
z_ratet <- z_table/sum(z_table)


####発生させたトピックから応答変数を発生####
##補助変数を発生
#連続変数の発生
k1 <- 3
X1 <- matrix(rnorm(d*k1, 0, 1), nrow=d, ncol=k1)

#二値変数の発生
k2 <- 3
X2 <- matrix(0, d, k2)
for(i in 1:k2){
  X2[, i] <- rbinom(d, 1, runif(1, 0.3, 0.7))
}

#トピックの不要変数を削除
ColSums(Rate1)
Rate1 <- Rate[, -which.min(colSums(Rate))]

#データを結合
X <- cbind(1, X1, X2)
XZ <- cbind(X, Rate1)


##ロジスティック回帰モデルから応答変数を発生
#パラメータの設定
b0 <- runif(1, -1.2, -0.5)
b1 <- c(runif(k1, 0, 0.9), runif(k2, -0.9, 1.1))
b2 <- runif(k-1, -2.5, 4.0)
b <- c(b0, b1, b2)

#確率の計算
logit <- XZ %*% b   #ロジット
Pr <- exp(logit)/(1+exp(logit))

#二項分布から応答変数を発生
y <- rbinom(d, 1, Pr)

loglt <- as.numeric(logLik(glm(y ~ ., data.frame(XZ[, -1]), family=binomial)))

res1 <- glm(y ~ ., data.frame(XZ[, -1]), family=binomial)
res2 <- glm(y ~ ., data.frame(XZ[, -c(1, 8:11)]), family=binomial)
summary(res1)
summary(res2)


####マルコフ連鎖モンテカルロ法でトピックモデル+ロジスティック回帰モデルを推定####
##アルゴリズムの設定
R <- 10000   #サンプリング回数 
keep <- 2
iter <- 0

##事前分布のパラメータの設定
#ハイパーパラメータの事前分布のパラメータ
alpha01 <- alpha0   #文書のディレクリ事前分布のパラメータ
beta01 <- alpha1[1]   #単語のディレクリ事前分布のパラメータ

#ロジステック回帰モデルの事前分布のパラメータ
betas <- rep(0, ncol(XZ))  #回帰係数の初期値
betaBi <- rep(0, ncol(XZ))   #回帰係数の事前分布の平均
rootBi <- 0.01*diag(ncol(XZ))   #事前分布の精度

##パラメータの初期値の設定
theta.ini <- runif(k, 0.3, 1.5)
theta <- rdirichlet(d, theta.ini)   #ユーザートピックのパラメータの初期値
phi.ini <- runif(v, 0.5, 1)
phi <- rdirichlet(k, phi.ini)   #アイテムトピックのパラメータの初期値

##パラメータの格納用配列
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
W.SEG <- matrix(0, nrow=sum(w), R/keep)
RATE <- matrix(0, nrow=R/keep, ncol=k)
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XZ))


####MCMC推定のためのデータの準備####
#IDを作成
d.id <- rep(1:d, rep(v, d))
w.id <- rep(1:v, d)
ID <- data.frame(d.id, w.id)

#インデックスを作成
index_w <- list()
index_v <- list()

for(i in 1:d) {index_w[[i]] <- subset(1:length(d.id), d.id==i)}
for(i in 1:v) {index_v[[i]] <- subset(1:length(w.id), w.id==i)}


##トピック割当の初期値を生成
Zx <- matrix(0, nrow=nrow(ID), ncol=k)

#ユーザーごとにアイテムのトピックを生成
for(i in 1:d){
  theta.m <- matrix(theta[i, ], nrow=k, ncol=v)   #thetaを行列形式に変更 
  z.rate <- t(phi * theta.m) / matrix(rowSums(t(phi * theta.m)), nrow=v, ncol=k)   #混合率を計算
  Zx[index_w[[i]], ] <- t(apply(cbind(WX[i, ], z.rate), 1, function(x) rmultinom(1, x[1], x[-1])))   #単語のトピック割当を決定
}

#トピック割当数の初期値を計算
#全体でのトピック割当
k_sum <- colSums(Zx) 

#ユーザーごとのトピック割当
kw_sum <- matrix(0, nrow=d, ncol=k)
for(i in 1:d) {kw_sum[i, ] <- colSums(Zx[index_w[[i]], ])}

#アイテムのトピック割当
kv_sum <- matrix(0, nrow=v, ncol=k)
for(i in 1:v) {kv_sum[i, ] <- colSums(Zx[index_v[[i]], ])}


#トピック割当をベクトル形式に変更
seg_vec <- unlist(apply(Zx, 1, function(x) rep(1:k, x)))   

#行列形式に変更
seg_mx <- matrix(0, nrow=length(seg_vec), ncol=k)
for(i in 1:nrow(seg_mx)) {seg_mx[i, seg_vec[i]] <- 1}


#トピック割当ベクトルのIDを作成
id_vec11 <- rep(1:d, w)
id_vec12 <- c()
for(i in 1:d) {id_vec12 <- c(id_vec12, rep(1:v, rowSums(Zx[index_w[[i]], ])))}

Z1 <- matrix(0, nrow=length(id_vec11), ncol=k)


#ギブスサンプリング用のインデックスを作成
index_user <- list()
for(i in 1:d) {index_user[[i]] <- subset(1:length(id_vec11), id_vec11==i)}


##ロジスティック回帰モデルの対数尤度を定義
loglike <- function(b, X, Y){
  
  #尤度を定義して合計する
  logit <- X %*% b 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}


####マルコフ連鎖モンテカルロ法でトピックモデル+ロジスティック回帰モデルを推定####
for(rp in 1:R){
  
  ##トピックをサンプリング
  for(i in 1:d){
    ##単語のトピックをサンプリング
    for(wd in 1:length(index_user[[i]])){
      index <- index_user[[i]][wd]   #単語のインデックス
      
      #トピック生成する単語のトピックを取り除く
      mx <- seg_mx[index, ]
      k1 <- k_sum - mx
      kw <- kw_sum[i, ] - mx
      kv <- kv_sum[id_vec12[index], ] - mx
    id_vec12
    kv_sum
    
      #単語のトピック割当確率を計算
      z_sums <- (kw + alpha01) * (kv + beta01) / (k1 + beta01*v)
      z_rate <- z_sums / sum(z_sums)
      
      #トピックをサンプリング
      Z1 <- t(rmultinom(1, 1, z_rate))
      
      #データを更新
      k_sum <- k1 + Z1
      kw_sum[i, ] <- kw + Z1
      kv_sum[id_vec12[index], ] <- kv + Z1
      seg_mx[index, ] <- Z1
      
    }
  }
  
  ##ロジスティック回帰モデルのパラメータをサンプリング
  #トピックを集計
  ZL <- matrix(0, nrow=d, ncol=k)
  
  for(i in 1:d){
    z_ind <- colSums(seg_mx[index_user[[i]], ])
    ZL[i, ] <- z_ind / sum(z_ind)
  }
  XZ <- cbind(X, ZL[, -k])   #データを結合
  
  #betaのサンプリング
  if(rp %in% c(1, 10, 100, 1000)){
    res <- glm(y ~ ., data=data.frame(XZ[, -1]), family=binomial)
    invBi <- diag(summary(res)[[12]][, 2]^2)
  }
  
  #ランダムウォークサンプリング
  betad <- betas
  betan <- betad + 0.1 * mvrnorm(1, rep(0, length(betad)), invBi)   #新しいbetaをランダムウォークでサンプリング
  
  #対数尤度と対数事前分布の計算
  lognew <- loglike(betan, XZ, y)
  logold <- loglike(betad, XZ, y)
  logpnew <- lndMvn(betan, betaBi, rootBi)
  logpold <- lndMvn(betad, betaBi, rootBi)
  
  #MHサンプリング
  alpha <- min(1, exp(lognew + logpnew - logold - logpold))
  if(alpha == "NAN") alpha <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しいbetaを採択
  if(u < alpha){
    betas <- betan
    logl <- lognew
    
    #そうでないならbetaを更新しない
  } else {
    logl <- logold
  } 
  
  ##サンプリング結果を保存
  mkeep <- rp/keep
  if(rp%%keep==0){
    
    #混合率の計算
    rate <- colSums(seg_mx)/nrow(seg_mx)
    
    #サンプリング結果を保存
    W.SEG[, mkeep] <- seg_mx %*% 1:k 
    RATE[mkeep, ] <- rate
    
    #サンプリング状況を表示
    print(rp)
    print(round(rbind(rate1=rate, ratet=z_ratet), 3))
    print(round(c(logl, loglt), 2))
    print(round(rbind(betas, b), 3))
  }
}



