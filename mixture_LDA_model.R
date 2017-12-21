#####多項LDAモデル#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(2578)

####データの発生####
##データの設定
k1 <- 10
k2 <- 8
d <- 2000   #文書数
v1 <- 300   #内容に関係のある語彙数
v2 <- 100   #内容に関係のない語彙数
v <- v1 + v2   #語彙数 
s <- rpois(d, 11.5)   #文章数
s[s < 5] <- ceiling(runif(sum(s < 5), 5, 10))
a <- sum(s)   #総文章数
w <- rpois(a, 9)   #文章あたりの単語数
w[w < 5] <- ceiling(runif(sum(w < 5), 5, 10))
f <- sum(w)   #総単語数

#文書IDの設定
u_id <- rep(1:d, s)
t_id <- c()
for(i in 1:d){t_id <- c(t_id, 1:s[i])}

##パラメータを設定
#ディレクリ分布のパラメータ
alpha0 <- rep(0.3, k1)
alpha1 <- c(rep(0.4, v1), rep(0.005, v2))
alpha2 <- c(rep(0.1, v1), rep(5.0, v2))

#ディレクリ分布よりパラメータを生成
thetat <- theta <- extraDistr::rdirichlet(d, alpha0)
phit <- phi <- extraDistr::rdirichlet(k1, alpha1)
gammat <- gamma <- extraDistr::rdirichlet(1, alpha2)


##文章ごとに単語を生成する
WX <- matrix(0, nrow=s, ncol=v)
x <- rep(0, f)
Z <- list()
index_v1 <- 1:v1
index_v2 <- (v1+1):v 

for(i in 1:d){
  i <- 1
  x <- rmnom(s[i], 1, theta[i, ])
}

rbind(theta[1, ], extraDistr::rdirichlet(1, colSums(x)*10+0.5))

