#####入れ子型マルチレベルモデル#####
library(MASS)
library(nlme)
library(glmm)
library(reshape2)
library(plyr)
library(lattice)
library(ggplot2)

####データの発生####
#set.seed(94327)
n <- 1000   #評価者数
g1 <- 50   #対象アニメ数
g2 <- round(runif(g1, 2.3, 5.8))   #登場キャラ数
g2s <- sum(g2)

####対象アニメを視聴しているかどうかを発生####
##説明変数の発生
cont <- 3; bin <- 4; multi <- 4
X.cont <- matrix(rnorm(n*cont), nrow=n, ncol=cont)
X.bin <- matrix(0, nrow=n, ncol=bin)
X.multi <- matrix(0, nrow=n, ncol=multi)

#二値説明変数を設定
for(i in 1:bin){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(n, 1, p)
}

#多値説明変数を設定
p <- runif(multi)
X.multi <- t(rmultinom(n, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))] #冗長な変数は削除

#データを結合
X <- cbind(1, X.cont, X.bin, X.multi)

##アニメの割当を発生
#パラメータの設定
alpha0 <- runif(g1, -1.5, 0.9)
alpha1 <- matrix(runif(g1*cont, 0, 0.7), nrow=cont, ncol=g1)
alpha2 <- matrix(runif(g1*(bin+multi-1), -0.8, 0.6), nrow=bin+multi-1, ncol=g1)
alpha <- rbind(alpha0, alpha1, alpha2)

#ロジットと確率の計算
logit <- X %*% alpha
Pr <- exp(logit)/(1+exp(logit))

#二項分布から割当を発生
R <- apply(Pr, 2, function(x) rbinom(n, 1, x))
colMeans(R); mean(R)


####デザイン行列を定義####
#デザイン行列の格納用配列
Z1 <- matrix(0, nrow=n*g2s, ncol=n)
Z2 <- matrix(0, nrow=n*g2s, ncol=g1)
Z3 <- matrix(0, nrow=n*g2s, ncol=g2s)

#インデックスを作成
index_g21 <- c(1, cumsum(g2))
index_g21[2:length(index_g21)] <- index_g21[2:length(index_g21)] + 1
index_g22 <- cumsum(g2)

for(i in 1:n){
  print(i)
  #個人別の格納用配列
  z2 <- matrix(0, nrow=g2s, ncol=g1)
  z3 <- matrix(0, nrow=g2s, ncol=g2s)
  
  r <- ((i-1)*g2s+1):((i-1)*g2s+g2s)
  Z1[r, i] <- 1
  
  for(j in 1:g1){
    if(R[i, j]==1){
      z2[index_g21[j]:index_g22[j], j] <- 1
      z3[index_g21[j]:index_g22[j], index_g21[j]:index_g22[j]] <- diag(g2[j])
    }
  }
  Z2[r, ] <- z2
  Z3[r, ] <- z3
}

Z <- cbind(Z1, Z2, Z3)   #データを結合

#評価していないアニメは欠損させる
index_zeros <- subset(1:nrow(Z2), rowSums(Z2)==0)
Z <- Z[-index_zeros, ]

##IDの設定
#ユーザーIDを設定
freq <- colSums(Z[, 1:n])
u.id <- rep(1:n, freq)

#評価回数を設定
t.id <- c()
for(i in 1:n) {t.id <- c(t.id, 1:freq[i])}

#アニメIDを設定
a.id <- rep(0, nrow(Z))
anime <- Z[, (n+1):(n+1+g1-1)]

for(i in 1:ncol(anime)){
  index <- subset(1:nrow(anime), anime[, i] > 0)
  a.id[index] <- i
}

#キャラIDを設定
c.id <- rep(0, nrow(Z))
chara <- Z[, (n+1+g1):(ncol(Z))]

for(i in 1:ncol(chara)){
  index <- subset(1:nrow(chara), chara[, i] > 0)
  c.id[index] <- i
}

#IDを結合
ID <- data.frame(no=1:nrow(Z), t=t.id, u.id=u.id, a.id=a.id, c.id=c.id)
table(ID$a.id); table(ID$c.id)


##固定効果の説明変数をパネル形式に変更
XM <- list()
for(i in 1:n) {XM[[i]] <- matrix(X[i, ], nrow=freq[i], ncol=ncol(X), byrow=T)}
X.panel <- do.call(rbind, XM)


####応答変数(評価データ)の発生####
##パラメータの設定
#固定効果のパラメータ
b.fix <- c(runif(1, 0.4, 0.7), runif(cont, 0, 0.6), runif(bin, -0.5, 0.7), runif(multi-1, -0.6, 0.9))   

#変量効果のパラメータ
random1 <- 1.0; random2 <- 1.5; random3 <- 1.25
b.g1 <- rnorm(n, 0, random1)
b.g2 <- rnorm(g1, 0, random2)
b.g3 <- rnorm(g2s, 0, random3)
b.random <- c(b.g1, b.g2, b.g3)

#個体内分散の設定
Cov <- 0.8

##応答変数を発生
mu <- X.panel %*% b.fix  + Z %*% b.random   #平均構造
y <- mu + rnorm(length(mu), 0, Cov)   #平均構造 + 誤差


####マルコフ連鎖モンテカルロ法でマルチレベルモデルの推定####
#アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4

##事前分布の設定
#固定効果の事前分布
beta_prior <- rep(0, ncol(X))   #固定効果の回帰係数の事前分布の平均
sigma_prior <- 0.01*diag(ncol(X))   #固定効果の回帰係数の事前分布の分散
tau_prior1 <- 0.01   #逆ガンマ分布の形状パラメータ
tau_prior2 <- 0.01   #逆ガンマ分布のスケールパラメータ

#変量効果の事前分布
alpha_random <- 0   #変量効果の事前分布の平均
tau_random1 <- 1   #逆ガンマ分布の形状パラメータ
tau_random2 <- 0.01   #逆ガンマ分布のスケールパラメータ


