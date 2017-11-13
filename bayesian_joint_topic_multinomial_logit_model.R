#####対応トピックモデル#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
detach("package:gtools", unload=TRUE)
detach("package:bayesm", unload=TRUE)
library(extraDistr)
library(monomvn)
library(glmnet)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

set.seed(8079)

####データの発生####
#set.seed(423943)
#文書データの設定
k <- 10   #トピック数
d <- 7000   #文書数
v <- 300   #語彙数
w <-  rpois(d, rgamma(d, 30, 0.5))   #1文書あたりの単語数
a <- 100   #補助変数数
x <- rpois(d, 20)
select <- 7   #選択肢数

#パラメータの設定
alpha0 <- rep(0.2, k) #文書のディレクリ事前分布のパラメータ
alpha1 <- rep(0.15, v)   #単語のディレクリ事前分布のパラメータ
alpha2 <- rep(0.25, a)   #補助情報のディクレリ事前分布のパラメータ

#ディレクリ乱数の発生
theta0 <- theta <- extraDistr::rdirichlet(d, alpha0)   #文書のトピック分布をディレクリ乱数から発生
phi0 <- phi <- extraDistr::rdirichlet(k, alpha1)   #単語のトピック分布をディレクリ乱数から発生
omega0 <- omega <- extraDistr::rdirichlet(k, alpha2) #補助データのトピック分布をディクレリ乱数から発生

#多項分布の乱数からデータを発生
WX <- matrix(0, nrow=d, ncol=v)
AX <- matrix(0, nrow=d, ncol=a)
Z1 <- list()
Z2 <- list()

for(i in 1:d){
  print(i)
  
  #文書のトピック分布を発生
  z1 <- t(rmultinom(w[i], 1, theta[i, ]))   #文書のトピック分布を発生
  
  #文書のトピック分布から単語を発生
  zn <- z1 %*% c(1:k)   #0,1を数値に置き換える
  zdn <- cbind(zn, z1)   #apply関数で使えるように行列にしておく
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi[x[1], ])))   #文書のトピックから単語を生成
  wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
  WX[i, ] <- wdn  
  
  #文書のトピック分布から補助変数を発生
  z2 <- t(rmultinom(x[i], 1, theta[i, ]))
  zx <- z2 %*% 1:k
  zax <- cbind(zx, z2)
  an <- t(apply(zax, 1, function(x) rmultinom(1, 1, omega[x[1], ])))
  adn <- colSums(an)
  AX[i, ] <- adn
  
  #文書トピックおよび補助情報トピックを格納
  Z1[[i]] <- z1
  Z2[[i]] <- z2
}

#データ行列を整数型行列に変更
storage.mode(WX) <- "integer"
storage.mode(AX) <- "integer"

####応答変数の発生####
#応答変数の格納用配列
y <- matrix(0, nrow=d, ncol=select)
Pr <- matrix(0, nrow=d, ncol=select)
Pr0 <- matrix(0, nrow=d, ncol=select)

##妥当な応答変数が発生するまで反復させる
for(j in 1:5000){
  ##パラメータの設定
  #トピックモデルのパラメータ
  sparse1 <- matrix(rbinom((select-1)*k, 1, 0.4), nrow=k, ncol=select-1)   #パラメータのスパース行列
  b00 <- runif(select-1, -0.5, 0.5)
  b01 <- (matrix(runif((select-1)*k, -4.0, 4.0), nrow=k, ncol=select-1)) * sparse1
  b02 <- (b01 + mvrnorm(k, rep(0, select-1), diag(0.2, select-1))) * sparse1
  b0 <- rbind(b00, b01, b02)
  rownames(b0) <- NULL
  
  ##文書ごとに確率と応答変数を発生
  for(i in 1:d){
    #潜在トピックzを集計
    z1 <- colSums(Z1[[i]])
    z1 <- ifelse(z1==0, 1, z1)
    z2 <- colSums(Z2[[i]])
    z2 <- ifelse(z2==0, 1, z2)
    
    #多項分布から応答変数を発生
    logit <- c(c(1, log(z1), log(z2)) %*% b0, 0)
    Pr[i, ] <- exp(logit) / sum(exp(logit))
    y[i, ] <- t(rmultinom(1, 1, Pr[i, ]))
  }
  
  t1 <- sum(apply(Pr, 1, which.max)==y %*% 1:select)/d
  print(round(c(t1, min(colSums(y))), 3))
  
  if(t1 > 0.85 &  min(colSums(y)) > 250) break
}



####マルコフ連鎖モンテカルロ法で結合トピックモデル+L1正則化ロジットモデルを推定####
####トピックモデルのためのデータと関数の準備####
##それぞれの文書中の単語の出現および補助情報の出現をベクトルに並べる
##データ推定用IDを作成
ID1_list <- list()
wd_list <- list()
ID2_list <- list()
ad_list <- list()

#求人ごとに求人IDおよび単語IDを作成
for(i in 1:nrow(WX)){
  print(i)
  
  #単語のIDベクトルを作成
  ID1_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
  num2 <- subset(num1, num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
  
  #補助情報のIDベクトルを作成
  ID2_list[[i]] <- rep(i, x[i])
  num1 <- (AX[i, ] > 0) * (1:a)
  num2 <- subset(num1, num1 > 0)
  A1 <- AX[i, (AX[i, ] > 0)]
  number <- rep(num2, A1)
  ad_list[[i]] <- number
}

#リストをベクトルに変換
ID1_d <- unlist(ID1_list)
ID2_d <- unlist(ID2_list)
wd <- unlist(wd_list)
ad <- unlist(ad_list)

##インデックスを作成
doc1_list <- list()
word_list <- list()
doc2_list <- list()
aux_list <- list()

for(i in 1:length(unique(ID1_d))) {doc1_list[[i]] <- which(ID1_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- which(wd==i)}
for(i in 1:length(unique(ID2_d))) {doc2_list[[i]] <- which(ID2_d==i)}
for(i in 1:length(unique(ad))) {aux_list[[i]] <- which(ad==i)}
gc(); gc()


####マルコフ連鎖モンテカルロ法で対応トピックモデルを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #負担係数の格納用
  for(kk in 1:k){
    #負担係数を計算
    Bi <- rep(theta[, kk], w) * phi[kk, c(wd)]   #尤度
    Bur[, kk] <- Bi   
  }
  
  Br <- Bur / rowSums(Bur)   #負担率の計算
  r <- colSums(Br) / sum(Br)   #混合率の計算
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##多項ロジットモデルの対数尤度
fr <- function(beta, y, x, hh, select){
  
  #ロジットと確率の計算
  logit <- t(x %*% t(beta))
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #対数尤度を設定
  LLi <- rowSums(y*log(Pr)) 
  return(LLi)
}


##アルゴリズムの設定
R <- 20000   #サンプリング回数
keep <- 4   #4回に1回の割合でサンプリング結果を格納
iter <- 0

##学習データとテストデータに分割
index_test <- sample(1:d, 1000)
n1 <- d-length(index_test)
n2 <- length(index_test)
y_train <- y[-index_test, ]
y_test <- y[index_test, ]
y_vec <- y_test %*% 1:select


##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- rep(1.0, k)
beta0 <- rep(0.5, v)
gamma0 <- rep(0.5, a)
alpha01m <- matrix(alpha01, nrow=d, ncol=k, byrow=T)
beta0m <- matrix(beta0, nrow=v, ncol=k)
gamma0m <- matrix(gamma0, nrow=a, ncol=k)


##パラメータの初期値
#トピックモデルの初期値
theta.ini <- runif(k, 0.5, 2)
phi.ini <- runif(v, 0.5, 1)
omega.ini <- runif(a, 0.5, 1)
theta <- rdirichlet(d, theta.ini)   #文書トピックのパラメータの初期値
phi <- rdirichlet(k, phi.ini)   #単語トピックのパラメータの初期値
omega <- rdirichlet(k, omega.ini)   #補助情報トピックのパラメータの初期値

#多項ロジットモデルの初期値
par <- 2*k+1
beta0 <- scale(colSums(y_train))
oldbeta <- mvrnorm(n1, beta0[-select], diag(0.2, select-1))
oldtheta <- matrix(0, nrow=par, ncol=select-1)
oldcov <- diag(0.1, select-1)
inv_cov <- solve(oldcov)
lambda <- c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01)


##パラメータの格納用配列
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
OMEGA <- array(0, dim=c(k, a, R/keep))
W_SEG <- matrix(0, nrow=R/(keep*10), ncol=sum(w))
A_SEG <- matrix(0, nrow=R/(keep*10), ncol=sum(x))
storage.mode(W_SEG) <- "integer"
storage.mode(A_SEG) <- "integer"
gc(); gc()

##MCMC推定用配列
wsum0 <- matrix(0, nrow=d, ncol=k)
vf0 <- matrix(0, nrow=v, ncol=k)
af0 <- matrix(0, nrow=a, ncol=k)
asum0 <- matrix(0, nrow=d, ncol=k)
x_diag <- diag(select)[, -select]


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##単語トピックをサンプリング
  #単語ごとにトピックの出現率を計算
  word_rate <- burden_fr(theta, phi, wd, w, k)$Br
  
  #多項分布から単語トピックをサンプリング
  Zi1 <- rmnom(nrow(word_rate), 1, word_rate) 
  
  ##文書トピックのパラメータを更新
  #ディクレリ分布からthetaをサンプリング
  for(i in 1:d){
    wsum0[i, ] <- colSums(Zi1[doc1_list[[i]], ]) 
  }
  wsum <- wsum0 + alpha01m   #ディクレリ分布のパラメータ
  theta <- extraDistr::rdirichlet(d, wsum)   #ディクレリ分布からthetaをサンプリング
  
  ##単語トピックのパラメータを更新
  #ディクレリ分布からphiをサンプリング
  for(i in 1:v){
    vf0[i, ] <- colSums(Zi1[word_list[[i]], ])
  }
  vf <- t(vf0 + beta0m)   #ディクレリ分布のパラメータ
  phi <- extraDistr::rdirichlet(k, vf)   #ディクレリ分布からphiをサンプリング

  
  ##補助情報トピックをサンプリング
  #文書トピックを用いて補助情報ごとにトピックの出現率を計算
  aux_rate <- burden_fr(theta, omega, ad, x, k)$Br
  
  #多項分布から補助情報トピックをサンプリング
  Zi2 <- rmnom(nrow(aux_rate), 1, aux_rate) 
  
  
  ##補助情報トピックのパラメータを更新
  #ディクレリ分布からomegaをサンプリング
  for(i in 1:a){
    af0[i, ] <- colSums(Zi2[aux_list[[i]], ])
  }
  af <- t(af0 + gamma0m)   #ディクレリ分布のパラメータ
  omega <- extraDistr::rdirichlet(k, af)   #ディクレリ分布からomegaをサンプリング
  
  
  ####ベイジアンL1正則化ロジットモデルでbetaをサンプリング####
  ##データの設定
  #文書トピックの説明変数の設定
  wsum0[wsum0==0] <- 1
  Data01 <- log(wsum0)
  
  #補助情報トピックの説明変数の設定
  for(i in 1:d){
    asum0[i, ] <- colSums(Zi2[doc2_list[[i]], ])
  }
  asum0[asum0==0] <- 1
  Data02 <- log(asum0)
  
  #学習データとテストデータに分割
  Data <- cbind(1, Data01, Data02)
  Data1 <- Data[-index_test, ]
  Data2 <- Data[index_test, ]

  
  ##多項ロジットモデルのパラメータをサンプリング
  #新しいパラメータをサンプリング
  betad <- oldbeta
  betan <- betad + mvrnorm(n1, rep(0, select-1), diag(0.25, select-1))
  
  #誤差を設定
  er_new <- betan - mu
  er_old <- betad - mu

  #対数尤度と対数事前分布を計算
  lognew <- fr(betan, y_train, x_diag, n1, select)
  logold <- fr(betad, y_train, x_diag, n1, select)
  logpnew <- -0.5 * rowSums(er_new %*% inv_cov * er_new)
  logpold <- -0.5 * rowSums(er_old %*% inv_cov * er_old)  
  
  #メトロポリスヘイスティング法でパラメータの採択を決定
  rand <- runif(n1)   #一様分布から乱数を発生
  LLind_diff <- exp(lognew + logpnew - logold - logpold)   #採択率を計算
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=n1, ncol=select-1)
  oldbeta <- flag*betan + (1-flag)*betad   #alphaがrandを上回っていたら採択
  
  ##lassoで階層モデルの回帰パラメータをサンプリング
  for(j in 1:(select-1)){
    
    #lambdaが一定以上なら強制的に小さいlambdaにする
    if(lambda[j] > 0.5){
      lambda[j] <- 0.001
    }
    
    #ベイジアンlassoで階層モデルの回帰係数を更新
    res <- blasso(X=Data1[, -1], y=oldbeta[, j], beta=oldtheta[-1, j], lambda2=lambda[j], s2=diag(oldcov)[j], 
                  normalize=TRUE, T=2)
    oldtheta[, j] <- c(res$mu[2], res$beta[2, ])
    lambda[j] <- res$lambda2[2]
  }
  
  ##階層モデルの分散共分散行列を推定
  mu <- Data1 %*% oldtheta
  oldcov <- var(oldbeta - mu)
  inv_cov <- solve(oldcov)
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    mkeep1 <- rp/keep
    THETA[, , mkeep1] <- theta
    PHI[, , mkeep1] <- phi
    OMEGA[, , mkeep1] <- omega
    #if(rp%%(keep*10)==0){
    #  mkeep2 <- rp/(keep*5)
    #  W_SEG[mkeep2, ] <- word_z
    #  A_SEG[mkeep2, ] <- aux_z
    #}
    
    #サンプリング結果を確認
    print(rp)
    print(round(c(sum(lognew), mean(alpha)), 3))
    print(round(lambda, 4))
    print(round(cbind(theta[1:10, ], theta0[1:10, ]), 3))
    #print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
    #print(round(cbind(omega[, 1:10], omegat[, 1:10]), 3))
    
    #予測分布を計算
    logit <- Data2 %*% cbind(oldtheta, 0)
    Pr <- exp(logit) / rowSums(exp(logit))
    print(mean(y[index_test, ] %*% 1:select==apply(Pr, 1, which.max)))
  }
}


round(cbind(oldtheta, b0), 3)
round(Pr, 3)
round(cbind(Pr, y[index_test, ]), 3)
