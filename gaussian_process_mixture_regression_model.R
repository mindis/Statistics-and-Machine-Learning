#####ガウス過程混合回帰モデル#####
library(GPfit)
library(kernlab)
library(MASS)
library(Matrix)
library(matrixStats)
library(extraDistr)
library(reshape2)
library(dplyr)

#set.seed(85347)

####データの発生####
##データの発生
week <- 7   #周期成分数
n <- week*104    #サンプル数
Data <- cbind(1, 1:n/10, matrix(diag(week), nrow=n, ncol=week, byrow=T)[, -1])   #データフレーム

##説明変数をグラム行列に変換
Kern1 <- kernelMatrix(rbfdot(sigma = 0.0001), Data[1:(n/2), ])
Kern2 <- kernelMatrix(rbfdot(sigma = 0.0001), Data[1:(n/2), ])

y1 <- mvrnorm(1, rep(35, n/2), Kern1)
y2 <- mvrnorm(1, rep(35, n/2), Kern2)
plot(1:n, c(y1, y2), type="l")

Kern1



