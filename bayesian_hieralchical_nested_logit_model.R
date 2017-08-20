#####階層ベイズネステッドロジットモデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(extraDistr)
library(gtools)
library(mlogit)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
hh <- 500
member <- 9   #メンバー数
c_num <- 8   #衣装パターン
hhpt <- hh*member*c_num   #全変数数

##個体内説明変数の発生
#メンバーの説明変数の設定
Mus <- matrix(as.numeric(table(1:hhpt, rep(rep(1:member, rep(c_num, member)), hh))), nrow=hhpt, ncol=member)
colnames(Mus) <- c("hono", "koto", "umi", "rin", "hana", "maki", "nico", "nozo", "eri")

#衣装の説明変数の設定
ct <- c(1, rep(0, c_num))
cloth <- matrix(ct, nrow=hh*member*(c_num+1), ncol=c_num, byrow=T)
CLOTH <- subset(cloth, rowSums(cloth) > 0)[, -c_num]
colnames(CLOTH) <- c("A", "B", "C", "D", "F", "G", "H")

#カード種別
k <- 3   #種類数
card <- t(rmultinom(hhpt, 1, c(2/c_num, 2/c_num, (c_num-4)/c_num)))
colnames(card) <- c("UR", "SSR", "SR")
CARD <- card[, -k]

#プロモーション接触数
Prom <- rpois(hhpt, 5)

#データの結合
X <- data.frame(Mus, CLOTH, CARD, Prom)
XM <- as.matrix(X)


##個体間説明変数の発生
