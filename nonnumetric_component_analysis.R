#####非計量主成分分析#####
library(homals)
library(MASS)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)


####データの読み込み####
data(galo)
Data <- galo

####データクレンジング####
##カテゴリカル変数をダミー変数に変換
#adviceをダミー行列に変換
advice_table <- table(1:nrow(Data), Data$advice)
advice_name <- colnames(advice_table)
advice <- matrix(as.numeric(advice_table), nrow=nrow(Data), ncol=ncol(advice_table))
colnames(advice) <- advice_name

#SESをダミー行列に変換
SES_table <- table(1:nrow(Data), Data$SES)
SES_name <- colnames(SES_table)
SES <- matrix(as.numeric(SES_table), nrow=nrow(Data), ncol=ncol(SES_table))
colnames(SES) <- SES_name

#性別を変換
sex <- ifelse(Data$gender=="M", 1, 0)   

##データの結合
X = cbind(sex, IQ=Data$IQ, advice, SES, school=Data$School)


####交互最小二乗法で非計量主成分分析を推定####
##パラメータの初期値を設定

