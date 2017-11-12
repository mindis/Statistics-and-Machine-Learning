#####結合確率的潜在意味解析#####
library(MASS)
library(lda)
library(RMeCab)
detach("package:bayesm", unload=TRUE)
library(extraDistr)
library(matrixStats)
library(monomvn)
library(lars)
library(glmnet)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(58079)

####データの発生####
#set.seed(423943)
#文書データの設定
k <- 15   #トピック数
d <- 8000   #文書数
v <- 450   #語彙数
w0 <- rpois(d, rgamma(d, 12.5, 0.175))   #1文書あたりの単語数
w <- ifelse(w0 < 20, 20, w0)
a <- 150   #補助変数数
x0 <- rpois(d, rgamma(d, 12.5, 1.25))   #1補助文書あたりの単語数
x <- ifelse(x0 < 2, 2, x0)
select <- 8   #応答変数数


#パラメータの設定
alpha0 <- rep(0.2, k)   #文書のディレクリ事前分布のパラメータ
alpha1 <- rep(0.15, v)   #単語のディレクリ事前分布のパラメータ
alpha2 <- rep(0.1, a)   #補助データのディクレリ事前分布のパラメータ

#ディレクリ乱数の発生
theta <- extraDistr::rdirichlet(d, alpha0)   #文書のトピック分布をディレクリ乱数から発生
phi <- extraDistr::rdirichlet(k, alpha1)   #単語のトピック分布をディレクリ乱数から発生
omega <- extraDistr::rdirichlet(k, alpha2)   #補助データのトピック分布をディクレリ乱数から発生


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
  
  #単語の変量効果のパラメータ
  sparse2 <- matrix(rbinom((select-1)*v, 1, 0.3), nrow=v, ncol=select-1)   #パラメータのスパース行列
  cov0 <- diag(runif(select-1, 0.025, 0.25))
  a0 <- mvrnorm(v, rep(0, select-1), cov0) * sparse2

  ##文書ごとに確率と応答変数を発生
  for(i in 1:d){
    logit <- c(c(1, log(colSums(Z1[[i]])+1), log(colSums(Z2[[i]])+1)) %*% b0, 0)
    Pr[i, ] <- exp(logit) / sum(exp(logit))
    y[i, ] <- t(rmultinom(1, 1, Pr[i, ]))
  }
  
  t1 <- sum(apply(Pr, 1, which.max)==y %*% 1:select)/d
  print(round(c(t1, min(colSums(y))), 3))
  
  if(t1 > 0.85 &  min(colSums(y)) > 250) break
}


####EMアルゴリズムでトピックモデルを推定####
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
for(i in 1:length(unique(ID1_d))) {doc1_list[[i]] <- subset(1:length(ID1_d), ID1_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- subset(1:length(wd), wd==i)}
for(i in 1:length(unique(ID2_d))) {doc2_list[[i]] <- subset(1:length(ID2_d), ID2_d==i)}
for(i in 1:length(unique(ad))) {aux_list[[i]] <- subset(1:length(ad), ad==i)}
gc(); gc()


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

####EMアルゴリズムの初期値を設定する####
##初期値をランダムに設定
#phiの初期値
freq_v <- matrix(colSums(WX), nrow=k, ncol=v, byrow=T)   #単語の出現率
rand_v <- matrix(trunc(rnorm(k*v, 0, (colSums(WX)/2))), nrow=k, ncol=v, byrow=T)   #ランダム化
phi_r <- abs(freq_v + rand_v) / rowSums(abs(freq_v + rand_v))   #トピックごとの出現率をランダムに初期化

#thetaの初期値
theta_r <- rdirichlet(d, runif(k, 0.2, 4))   #ディレクリ分布から初期値を設定

#omegaの初期値
freq_v <- matrix(colSums(AX)/sum(AX), nrow=k, ncol=a, byrow=T)
rand_v <- matrix(trunc(rnorm(k*a, 0, (colSums(AX)/2))), nrow=k, ncol=a, byrow=T)   #ランダム化
omega_r <- abs(freq_v + rand_v) / rowSums(abs(freq_v + rand_v))   #トピックごとの出現率をランダムに初期化

###パラメータの更新
##単語レベルの負担率とパラメータを初期化
#単語レベルの負担率の更新
word_fr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, w=w, k=k)
Bw <- word_fr$Br   #負担率
r1 <- word_fr$r   #混合率

#thetaの更新
wsum <- (data.frame(id=ID1_d, Br=Bw) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_all(funs(sum)))[, 2:(k+1)]
theta_r <- wsum / matrix(w, nrow=d, ncol=k)   #パラメータを計算

##phiの更新
vf <- (data.frame(id=wd, Br=Bw) %>%
         dplyr::group_by(id) %>%
         dplyr::summarize_all(funs(sum)))[, 2:(k+1)]
phi_r <- t(vf) / matrix(colSums(vf), nrow=k, ncol=v)


##補助情報レベルの負担率とパラメータを初期化
#補助情報レベルの負担率の更新
aux_fr <- burden_fr(theta=theta_r, phi=omega_r, wd=ad, w=x, k=k)
Ba <- aux_fr$Br   #負担率
r2 <- aux_fr$r   #混合率

##omegaの更新
af <- (data.frame(id=ad, Br=Ba) %>%
         dplyr::group_by(id) %>%
         dplyr::summarize_all(funs(sum)))[, 2:(k+1)]
omega_r <- t(af) / matrix(colSums(af), nrow=k, ncol=a)


#対数尤度の計算
LLw <- sum(log(rowSums(word_fr$Bur)))   #単語レベルの対数尤度
LLa <- sum(log(rowSums(aux_fr$Bur)))   #補助情報レベルの対数尤度
LL <- LLw + LLa


####EMアルゴリズムでパラメータを更新####
#更新ステータス
iter <- 1
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 0.1
LL1 <- LL   #対数尤度の初期値
LLs <- c()

##EMアルゴリズムが収束するまで反復させる
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  
  ##単語レベルのパラメータを最尤推定
  #単語レベルの負担率の更新
  word_fr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, w=w, k=k)
  Bw <- word_fr$Br   #負担率
  r1 <- word_fr$r   #混合率
  
  #thetaの更新
  wsum <- (data.frame(id=ID1_d, Br=Bw) %>%
             dplyr::group_by(id) %>%
             dplyr::summarize_all(funs(sum)))[, 2:(k+1)]
  theta_r <- wsum / matrix(w, nrow=d, ncol=k)   #パラメータを計算
  
  ##phiの更新
  vf <- (data.frame(id=wd, Br=Bw) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_all(funs(sum)))[, 2:(k+1)]
  phi_r <- t(vf) / matrix(colSums(vf), nrow=k, ncol=v)
  
  
  ##補助情報レベルのパラメータを最尤推定
  #補助情報レベルの負担率の更新
  aux_fr <- burden_fr(theta=theta_r, phi=omega_r, wd=ad, w=x, k=k)
  Ba <- aux_fr$Br   #負担率
  r2 <- aux_fr$r   #混合率
  
  ##omegaの更新
  af <- (data.frame(id=ad, Br=Ba) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_all(funs(sum)))[, 2:(k+1)]
  omega_r <- t(af) / matrix(colSums(af), nrow=k, ncol=a)
  
  
  ##観測データの対数尤度の計算
  LLw <- sum(log(rowSums(word_fr$Bur)))   #単語レベルの対数尤度
  LLa <- sum(log(rowSums(aux_fr$Bur)))   #補助情報レベルの対数尤度
  LL <- LLw + LLa
  
  ##アルゴリズムの更新
  iter <- iter+1
  dl <- LL1 - LL
  LL1 <- LL
  LLs <- c(LLs, LL)
  print(LL)
}


####推定結果と統計量####
plot(1:length(LLs), LLs, type="l", xlab="iter", ylab="LL", main="対数尤度の変化", lwd=2)

(PHI <- data.frame(round(t(phi_r), 3), t=round(t(phi), 3)))   #phiの真の値と推定結果の比較
(OMEGA <- data.frame(round(t(omega_r), 3), t=round(t(omega), 3)))   #omegaの真の値と推定結果の比較
(THETA <- data.frame(w, round(theta_r, 3), t=round(theta, 3)))   #thetaの真の値と推定結果の比較
r   #混合率の推定結果

round(colSums(THETA[, 2:(k+1)]) / sum(THETA[, 2:(k+1)]), 3)   #推定された文書中の各トピックの比率
round(colSums(THETA[, (k+1):(2*k)]) / sum(THETA[, (k+1):(2*k)]), 3)   #真の文書中の各トピックの比率

#AICとBIC
tp <- dim(theta_r)[1]*dim(theta_r)[2] 
pp <- dim(phi_r)[1]*dim(phi_r)[2] + dim(omega_r)[1]*dim(omega_r)[2]

(AIC <- -2*LL + 2*(tp+pp)) 
(BIC <- -2*LL + log(nrow(WX))*(tp+pp))

##結果をグラフ化
#thetaのプロット(50番目の文書まで)
barplot(theta_r[1:100, 1], ylim=c(0, 1), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(theta_r[, 1]), col=10, lty=5)
barplot(theta_r[1:100, 2], ylim=c(0, 1), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(theta_r[, 2]), col=10, lty=5)
barplot(theta_r[1:100, 3], ylim=c(0, 1), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(theta_r[, 2]), col=10, lty=5)
barplot(theta_r[1:100, 4], ylim=c(0, 1), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(theta_r[, 2]), col=10, lty=5)
barplot(theta_r[1:100, 5], ylim=c(0, 1), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(theta_r[, 2]), col=10, lty=5)

#phiのプロット(50番目の単語まで)
barplot(phi_r[1, 1:50], ylim=c(0, 0.05), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(phi_r[1, ]), col=10, lty=5)
barplot(phi_r[2, 1:50], ylim=c(0, 0.05), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(phi_r[2, ]), col=10, lty=5)
barplot(phi_r[3, 1:50], ylim=c(0, 0.05), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(phi_r[3, ]), col=10, lty=5)
barplot(phi_r[4, 1:50], ylim=c(0, 0.05), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(phi_r[4, ]), col=10, lty=5)
barplot(phi_r[5, 1:50], ylim=c(0, 0.05), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(phi_r[5, ]), col=10, lty=5)

#omegaのプロット
barplot(omega_r[1, ], col=c(1:ncol(omega_r)), density=50)
abline(h=mean(omega_r[1, ]), col=10, lty=5)
barplot(omega_r[2, ], col=c(1:ncol(omega_r)), density=50)
abline(h=mean(omega_r[2, ]), col=10, lty=5)
barplot(omega_r[3, ], col=c(1:ncol(omega_r)), density=50)
abline(h=mean(omega_r[3, ]), col=10, lty=5)
barplot(omega_r[4, ], col=c(1:ncol(omega_r)), density=50)
abline(h=mean(omega_r[4, ]), col=10, lty=5)
barplot(omega_r[5, ], col=c(1:ncol(omega_r)), density=50)
abline(h=mean(omega_r[5, ]), col=10, lty=5)


####階層ベイズ多項ロジットモデルで分類モデルを作成####
####データの設定####
##推定されたトピックから応答変数を作成
#文書トピックの出現を集計
W0 <- as.matrix((data.frame(id=ID1_d, Br=Bw) %>%
                   dplyr::group_by(id) %>%
                   dplyr::summarize_all(funs(sum))))[, 2:(k+1)]

#タイトルトピックの出現を集計
A0 <- as.matrix((data.frame(id=ID2_d, Ba=Ba) %>%
                   dplyr::group_by(id) %>%
                   dplyr::summarize_all(funs(sum))))[, 2:(k+1)]

#データを結合して説明変数とする
Data1 <- cbind(1, log(W0 + 1), log(A0 + 1))
Data <- matrix(0, nrow=d, ncol=3*k+1)
round(cbind(theta_r, theta), 3)

for(i in 1:d){
  Data[i, ] <- c(1, colSums(Z1[[i]]), log(colSums(Z1[[i]])+1), log(colSums(Z2[[i]])+1))
}
round(cbind(W0, Data[, 2:21]), 0)

####マルコフ連鎖モンテカルロ法の設定####
##多項ロジットモデルの対数尤度
fr <- function(beta, y, x, hh, select){
  
  #ロジットと確率の計算
  logit <- t(x %*% t(beta))
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #対数尤度を設定
  LLi <- rowSums(y*log(Pr)) 
  return(LLi)
}

##学習データとテストデータに分割
index_test <- sample(1:nrow(Data), 1000)
n1 <- d-length(index_test)
n2 <- length(index_test)
Data_train <- Data[-index_test, ]
Data_test <- Data[index_test, ]
y_train <- y[-index_test, ]
y_test <- y[index_test, ]
y_vec <- y_test %*% 1:select

##アルゴリズムの設定
R <- 20000
keep <- 4
sbeta <- 1.5
iter <- 0

##初期値の設定
par <- 2*k+1
beta0 <- scale(colSums(y_train))
oldbeta <- mvrnorm(n1, beta0[-select], diag(0.2, select-1))
oldtheta <- matrix(0, nrow=par, ncol=select-1)

oldcov <- diag(0.1, select-1)
inv_cov <- solve(oldcov)
mu <- Data_train %*% oldtheta

##サンプリング結果の保存用配列
THETA <- array(0, dim=c(par, select-1, R/keep))
COV <- array(0, dim=c(select-1, select-1, R/keep))

##アルゴリズム推定用配列
lognew <- rep(0, n1)
logold <- rep(0, n1)
logpnew <- rep(0, n1)
logpold <- rep(0, n1)
lambda <- c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001)
x_diag <- diag(select)[, -select]
er_new <- matrix(0, nrow=n1, select-1)
er_old <- matrix(0, nrow=n1, select-1)


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##多項ロジットモデルのパラメータをサンプリング
  #新しいパラメータをサンプリング
  betad <- oldbeta
  betan <- betad + mvrnorm(n1, rep(0, select-1), diag(0.015, select-1))

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
    if(lambda[j] > 0.5){
      lambda[j] <- 0.001
    }
    res <- blasso(X=Data_train[, -1], y=oldbeta[, j], beta=oldtheta[-1, j], lambda2=lambda[j], s2=diag(oldcov)[j], 
                  normalize=TRUE, T=2)
    oldtheta[, j] <- c(res$mu[2], res$beta[2, ])
    lambda[j] <- res$lambda2[2]
  }
 
  ##階層モデルの分散共分散行列を推定
  mu <- Data_train %*% oldtheta
  oldcov <- var(oldbeta - mu)
  inv_cov <- solve(oldcov)
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- oldtheta
    COV[, , mkeep] <- oldcov
    print(rp)
    print(round(lambda, 3))
    print(round(mean(alpha), 3))
    print(sum(lognew))
    print(round(cbind(oldtheta[1:15, ], b0[1:15, ]), 3))

    ##予測分布を推定
    logit <- Data_test %*% oldtheta
    Pr <- exp(logit) / rowSums(exp(logit))
    print(mean(y_vec==apply(Pr, 1, which.max)))
  }
}

colSums(y)
round(cbind(oldtheta, b0), 3)



      