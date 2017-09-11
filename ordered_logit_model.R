#####順序ロジットモデル(比例オッズモデル)#####
library(MASS)
library(ordinal)
library(mlogit)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)
set.seed(74819)

####データの発生####
N <- 5000   #サンプル数
R <- 5   #順序数
k <- 12   #変数数

####説明変数の発生####
#連続変数の発生
cont <- 4
X.cont <- matrix(rnorm(N*cont, 0, 1), nrow=N, ncol=cont)

#二値変数の発生
bin <- 4
X.bin <- matrix(0, nrow=N, ncol=bin)
for(i in 1:bin2){
  X.bin[, i] <- rbinom(N, 1, runif(1, 0.35, 0.65))
}

#多値変数の発生
multi <- 4
r <- runif(multi, 0.3, 1.8)
p <- r/sum(r)
X.multi <- t(rmultinom(N, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))]

#データの結合
X <- cbind(X.cont, X.bin, X.multi)


####応答変数の発生####
##パラメータの設定
b0 <- c(1.8, 1.0, 0.5, -0.5)
b1 <- c(runif(cont, 0, 0.8), runif(bin+multi-1, -0.8, 0.9))
b0

logit <- X %*% b1
exp(1.8 + logit)
exp(1.8 - logit) 
exp(1.0 - logit)
exp(0.5 - logit)
exp(-0.5 - logit)
exp(-0.5 + logit)


