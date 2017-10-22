#####不釣り合いデータの確率的潜在意味解析(トピックモデル)#####
library(MASS)
library(lda)
library(RMeCab)
detach("package:bayesm", unload=TRUE)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####データの発生####
#set.seed(423943)
#データの設定
k <- 8   #トピック数
d <- 5000   #文書数
v <- 300   #語彙数
w <- rpois(d, 250)   #1文書あたりの単語数

#パラメータの設定
alpha0 <- round(runif(k, 0.1, 1.25), 3)   #文書のディレクリ事前分布のパラメータ
alpha1 <- rep(0.25, v)   #単語のディレクリ事前分布のパラメータ

#ディレクリ乱数の発生
theta <- rdirichlet(d, alpha0)   #文書のトピック分布をディレクリ乱数から発生
phi <- rdirichlet(k, alpha1)   #単語のトピック分布をディレクリ乱数から発生

#多項分布の乱数からデータを発生
WX <- matrix(0, nrow=d, ncol=v)
Z <- list()
for(i in 1:d){
  z <- t(rmultinom(w[i], 1, theta[i, ]))   #文書のトピック分布を発生
  zn <- z %*% c(1:k)   #0,1を数値に置き換える
  zdn <- cbind(zn, z)   #apply関数で使えるように行列にしておく
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi[x[1], ])))   #文書のトピックから単語を生成
  wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
  WX[i, ] <- wdn  
  Z[[i]] <- zdn[, 1]
  print(i)
}

####EMアルゴリズムでトピックモデルを推定####
####トピックモデルのためのデータと関数の準備####
##それぞれの文書中の単語の出現をベクトルに並べる
##データ推定用IDを作成
ID_list <- list()
wd_list <- list()

#求人ごとに求人IDおよび単語IDを作成
for(i in 1:nrow(WX)){
  print(i)
  ID_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * c(1:v) 
  num2 <- subset(num1, num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
}

#リストをベクトルに変換
ID_d <- unlist(ID_list)
wd <- unlist(wd_list)

##インデックスを作成
doc_list <- list()
word_list <- list()
for(i in 1:length(unique(ID_d))) {doc_list[[i]] <- subset(1:length(ID_d), ID_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- subset(1:length(wd), wd==i)}
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

###パラメータの更新
##負担率の計算
bfr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, w=w, k=k)
Br <- bfr$Br   #負担率
r <- bfr$r   #混合率

##thetaの更新
tsum <- (data.frame(id=ID_d, Br=Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
theta_r <- tsum / matrix(w, nrow=d, ncol=k)   #パラメータを計算

##phiの更新
vf <- (data.frame(id=wd, Br=Br) %>%
         dplyr::group_by(id) %>%
         dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
phi_r <- t(vf) / matrix(colSums(vf), nrow=k, ncol=v)

#対数尤度の計算
(LLS <- sum(log(rowSums(bfr$Bur))))


####EMアルゴリズムでパラメータを更新####
#更新ステータス
iter <- 1
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 0.1
LLo <- LLS   #対数尤度の初期値
LLw <- LLS

out <- cbind(ID_d, wd, round(bfr$Br, 3))

###パラメータの更新
##負担率の計算
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  bfr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, w=w, k=k)
  Br <- bfr$Br   #負担率
  r <- bfr$r   #混合率
  
  ##thetaの更新
  tsum <- (data.frame(id=ID_d, Br=Br) %>%
             dplyr::group_by(id) %>%
             dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
  theta_r <- tsum / matrix(w, nrow=d, ncol=k)   #パラメータを計算

  ##phiの更新
  vf <- (data.frame(id=wd, Br=Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
  phi_r <- t(vf) / matrix(colSums(vf), nrow=k, ncol=v)
  
  #対数尤度の計算
  LLS <- sum(log(rowSums(bfr$Bur)))
  
  iter <- iter+1
  dl <- LLS-LLo
  LLo <- LLS
  LLw <- c(LLw, LLo)
  print(LLo)
}

####推定結果と統計量####
plot(1:length(LLw), LLw, type="l", xlab="iter", ylab="LL", main="対数尤度の変化", lwd=2)

(PHI <- data.frame(round(t(phi_r), 3), t=round(t(phi), 3)))   #phiの真の値と推定結果の比較
(THETA <- data.frame(w, round(theta_r, 3), t=round(theta, 3)))   #thetaの真の値と推定結果の比較
r   #混合率の推定結果

round(colSums(THETA[, 1:k]) / sum(THETA[, 1:k]), 3)   #推定された文書中の各トピックの比率
round(colSums(THETA[, (k+1):(2*k)]) / sum(THETA[, (k+1):(2*k)]), 3)   #真の文書中の各トピックの比率

#AICとBIC
tp <- dim(theta_r)[1]*dim(theta_r)[2]
pp <- dim(phi_r)[1]*dim(phi_r)[2]

(AIC <- -2*LLo + 2*(tp+pp)) 
(BIC <- -2*LLo + log(nrow(WX))*(tp+pp))

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