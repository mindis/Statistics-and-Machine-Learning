#####不完全観測データの生存時間解析#####
library(MASS)
library(matrixStats)
library(survival)
library(FAdist)
library(actuar)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(5362879)

####データの発生####
##データの設定
hh <- 3000
seg <- 2
dt <- 36   #観測期間
seg_id <- rep(1:seg, rep(hh/seg, seg))   #セグメントを設定

##説明変数の発生
type <- rbinom(hh, 1, 0.5)   #治療群の設定
level <- t(rmultinom(hh, 1, c(0.3, 0.25, 0.2, 0.15, 0.1))) %*% 1:5   #疾患レベル
Bun <- runif(hh, 0.5, 2)   #血中尿素窒素
Ca <- runif(hh, 0.1, 1.5)   #血清カルシウム
Hb <- runif(hh, 0.4, 1.5)   #血中ヘモグロビン
Prot <- rbinom(hh, 1, 0.4)   #ベンスジョーンズ蛋白の有無
sex <- rbinom(hh, 1, 0.6)   #性別
age <- t(rmultinom(hh, 1, c(0.2, 0.3, 0.3, 0.2)))   #年代
colnames(age) <- c("40代", "50代", "60代", "70代以上")
age <- age[, -1]

#データの結合
X <- data.frame(切片=1, type, level, Bun, Ca, Hb, Prot, sex, age)
XM <- as.matrix(X)

####応答変数と打ち切り変数の発生####
##パラメータの設定
#形状パラメータの設定
alpha0 <- c(runif(1, 5.5, 9), runif(1, 5, 8))

#スケールパラメータの設定
beta00 <- c(runif(1, 2.3, 2.7), runif(1, 2.8, 3.2))
beta01 <- t(matrix(c(runif(2, 0.1, 0.25), runif(2, 0.04, 0.1), runif(2, -0.25, -0.05), runif(2, 0.05, 0.3), runif(2, 0.1, 0.25), 
                     runif(2, -0.15, 0.15), runif(2, -0.2, 0.2), -0.1, -0.12, -0.16, -0.15, -0.2, -0.23), nrow=2, ncol=ncol(X)-1))
beta0 <- rbind(beta00, beta01)   #パラメータの結合  

##対数ロジスティック分布およびワイブル分布から生存時間を発生
#線形結合
scale1 <- exp(XM[seg_id==1, ] %*% beta0[, 1])   
scale2 <- exp(XM[seg_id==2, ] %*% beta0[, 2])

#応答変数を発生
y01 <- rllogis(length(scale1), shape=alpha[1], scale=scale[, 1])   #対数ロジスティック分布
y02 <- rweibull(length(scale2), shape=alpha[2], scale=scale[, 2])   #ワイブル乱数
y0 <- c(y01, y02)

#可視化
hist(y01, col="grey", breaks=25, main="対数ロジスティック分布からの応答変数", xlab="生存時間")
hist(y02, col="grey", breaks=25, main="ワイブル分布からの応答変数", xlab="生存時間")

##打ち切りの設定
#右側打ち切りの設定
y1 <- ifelse(y0 > dt, dt, y0)

#不完全観測の打ち切りを設定
section <- 3
y1_lower <- floor(y1/section) * section
y1_upper <- ceiling(y1/section) * section
round(cbind(y1, y1_lower, y1_upper), 3)   #データの確認


####EMアルゴリズムで不完全観測データの生存時間解析モデルを推定####
##観測データの対数尤度を定義


L1 <- dllogis(y2, 1, shape=alpha[1], scale=scale[, 1]) 
L2 <- dweibull(y2, alpha[2], scale[, 2])
round(rate1 <- L1/(L1 + L2), 3)
round(rate2 <- L2/(L1 + L2), 3)
mean(rate1)
mean(rate2)
