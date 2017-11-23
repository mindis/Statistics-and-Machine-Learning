#####ベイジアン潜在因子モデル#####
library(MASS)
library(lda)
library(RMeCab)
detach("package:bayesm", unload=TRUE)
library(extraDistr)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

set.seed(654978)

####データの発生####
#データを設定
hh <- 2000   #サンプル数
colums <- 200   #変数数
k <- 10   #潜在変数数

#バイナリ行列を発生
Z <- matrix(0, nrow=hh, ncol=k)
for(j in 1:k){
  p <- runif(1, 0.25, 0.5)
  Z[, j] <- rbinom(hh, 1, p)
}

#因子行列を発生
X <- matrix(rnorm(hh*k), nrow=k, ncol=colums)

#応答変数を発生
Cov <- diag(0.25, colums)
y = Z %*% X + mvrnorm(hh, rep(0, colums), Cov)


####マルコフ連鎖モンテカルロ法で潜在特徴モデルを推定####
##アルゴリズムの設定
R <- 10000
keep <- 2
iter <- 0

##事前分布の設定
#正規分布の分散の事前分布のパラメータ
s1 <- 0.01
v1 <- 0.01
s2 <- 0.01
v1 <- 0.01


##初期値の設定



