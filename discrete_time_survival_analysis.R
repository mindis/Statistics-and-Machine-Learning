#####離散時間生存モデル#####
library(MASS)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(324)
##データの設定
t <- 24*3   #半月単位の3年
n <- 500   #サンプル数
N <- t*n   #総サンプル数
col <- 15   #変数数

#id、time、整理番号を設定
id <- rep(1:n, rep(t, n))
time <- rep(1:t, n)
no <- 1:N

##パラメータの設定
#トレンド成分
T <- 10000
for(t1 in 1:T){
tb <- 0.4   #初期値
  trend <- numeric()
  s <- seq(0.7, 0.2, length=t)
  for(i in 1:t){
    r <- rnorm(5, tb, 0.04)
    sort <- sort(r)
    bi <- rbinom(1, 1, s[i])
    bb <- ifelse(bi == 1, sort[4], sort[2])
    tb <- bb
    trend <- c(trend, bb)
  }
  if(max(trend) > 0.55 && min(trend) < 0.1) break
  print(t1)
}  
plot(trend, type="l", lwd=1, xlab="half_month")

##回帰成分の設定
cont <- 8   #連続変数
cnt <- 2   #カウントデータ数
bin <- 4   #二値変数

#連続変数の発生
X_cont <- matrix(runif(cont*N, 0, 1), N, cont)   

#カウントデータの発生
X_cnt <- matrix(rpois(cnt*N, 1.5), N, cnt)   

#二値変数の発生
X_bin <- matrix(0, N, bin)
pr_b <- runif(bin, 0.3, 0.7)
for(i in 1:bin){
  bin <- rbinom(N, 1, pr_b[i])
  X_bin[, i] <- bin
}

#前回の購買からの経過時間
X_a <- rep(0, N)

X <- data.frame(cont=X_cont, cnt=X_cnt, bin=X_bin, hist=X_a)

##パラメータの設定
beta <- c(runif((col-1), -0.44, 0.14), 0.04)

##前回の購買からの期間の説明変数を記録しながら、シミュレーションデータを発生させる
y <- matrix(NA, N, 1)   #反応変数を格納
for(i in 1:n){
  for(j in 1:t){
    r <- (i-1)*t+j
    #購買データの発生
    xb <- trend[j] + as.matrix(X[r, ]) %*% beta   #線形関数
    y[r, ] <- rbinom(1, 1, exp(xb) / (1 + exp(xb)))   #購買を発生させる
    
    #前回からの購買月数を記録
    if(id[r]==n && time[r]==t) break   #最終行はbreak
    if(y[r, 1]==0 && sum(y[is.na(y)==0 & id==i, ])==0 && time[r]!=t) X$hist[r+1] <- 0 
    if(y[r, 1]==1 && time[r]!=t) X$hist[r+1] <- 0.5 
    if(y[r, 1]==0 && sum(y[is.na(y)==0 & id==i, ])!=0 && time[r]!=t) 
      {h <- subset(1:t, y[id==i, 1]==1); X$hist[r+1] <- ((j+1)-h[length(h)])/2} 
    if(time[r]==t) X$hist[r+1] <- 0
  }
    print(i)
}

##発生させたデータの確認と集計
round(YX <- data.frame(id, time, y, X), 2)   #データの確認
mean(y)   #購買率
summary(YX[, 3:length(YX)])   #データの要約
by(YX[, 3:length(YX)], YX$id, summary)   #個人ごとの要約
by(YX[, 3:length(YX)], YX$time, summary)   #時間ごとの要約
by(YX$y, YX$time, mean)   #時間ごとの要約

####離散時間生存モデルで推定####


