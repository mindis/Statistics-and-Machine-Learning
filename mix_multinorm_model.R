#####混合正規多変量回帰モデルの推定#####
library(MASS)
library(plyr)
library(lavaan)
####データの発生####
#set.seed(4543)
k <- 4   #混合数
col <- 5   #変数数
n <- 4000   #サンプル数

##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

#混合分布ごとの相関行列を作成
corM1 <- corrM(col=5, lower=-0.8, upper=0.8)
corM2 <- corrM(col=5, lower=-0.4, upper=0.5)
corM3 <- corrM(col=5, lower=-0.3, upper=0.6)
corM4 <- corrM(col=5, lower=-0.5, upper=0.4)

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, mu, sig){
  m <- abs(rnorm(col, mu, sig))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- m[i] * m[j]
    }
  }
  diag(c) <- m
  cc <- c * corM
  #固有値分解で強制的に正定値行列に修正する
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

#分散共分散行列を作成
Sigma1 <- covmatrix(col=5, corM=corM1, mu=0, sig=6)
Sigma2 <- covmatrix(col=5, corM=corM2, mu=0, sig=3)
Sigma3 <- covmatrix(col=5, corM=corM3, mu=0, sig=8)
Sigma4 <- covmatrix(col=5, corM=corM4, mu=0, sig=4)

Sigma1[-2]
Sigma2[-2]
Sigma3[-2]
Sigma4[-2]

##多変量正規分布の平均を回帰モデルで表現する
#入力するデータの発生
#性別
sex1 <- rbinom(1000, 1, 0.7)
sex2 <- rbinom(1000, 1, 0.5)
sex3 <- rbinom(1000, 1, 0.4)
sex4 <- rbinom(1000, 1, 0.3)
sex <- c(sex1, sex2, sex3, sex4)

#年代
age1 <- t(rmultinom(1000, 1, c(0.2, 0.3, 0.2, 0.2, 0.1)))
age2 <- t(rmultinom(1000, 1, c(0.1, 0.4, 0.3, 0.1, 0.1)))
age3 <- t(rmultinom(1000, 1, c(0.1, 0.2, 0.2, 0.3, 0.2)))
age4 <- t(rmultinom(1000, 1, c(0.4, 0.3, 0.1, 0.05, 0.05)))
age <- rbind(age1, age2, age3, age4)

#職業
job1 <- t(rmultinom(1000, 1, c(0.5, 0.2, 0.2, 0.1)))
job2 <- t(rmultinom(1000, 1, c(0.6, 0.4, 0.05, 0.05)))
job3 <- t(rmultinom(1000, 1, c(0.4, 0.4, 0.1, 0.1)))
job4 <- t(rmultinom(1000, 1, c(0.4, 0.3, 0.1, 0.2)))
job <- rbind(job1, job2, job3, job4)

#累積来店回数
cnt1 <- rpois(1000, 9)
cnt2 <- rpois(1000, 6)
cnt3 <- rpois(1000, 12)
cnt4 <- rpois(1000, 4)
cnt <- c(cnt1, cnt2, cnt3, cnt4)

#セグメント番号
segment <- rep(1:4, c(rep(1000, 4)))

#データを結合
data <- as.data.frame(cbind(segment, sex, age, job, cnt))
head(data, 10); tail(data, 10)

#変数に名前をつける
names(data)[3:12] <- c("age20", "age30", "age40", "age50", "age60u", "work", "stud", "hwife", "others", "s_cnt")
head(data, 10); tail(data, 10)
summary(data)


#基準変数を削除する
data1 <- data[, -7]   #60代以上を削除
data_n <- data1[, -10]   #その他の職業を削除
data_de <- cbind(1, data_n[, -1])
names(data_de)[1] <- ("intercept")

#セグメントごとの回帰係数を決定
v <- dim(data_de)[2]

a1 <- matrix(rnorm(v*col, 1.3, 1.2), nrow=v, ncol=col) 
a2 <- matrix(rnorm(v*col, 2.5, 2.3), nrow=v, ncol=col) 
a3 <- matrix(rnorm(v*col, 0.6, 1), nrow=v, ncol=col) 
a4 <- matrix(rnorm(v*col, 0.4, 2), nrow=v, ncol=col) 

a1; a2; a3; a4   #変数を確認

#回帰係数をリスト構造にする
a <- list(a1, a2, a3, a4)

#平均構造を多変量回帰モデルで発生させる
y <- matrix(0, 0, 11)
for(i in 1:k){
  seg_x <- subset(data_de, data_n[, 1] == i)
  y_temp <- as.matrix(seg_x) %*% as.matrix(a[[i]])
  y_temp <- as.data.frame(y_temp)
  y <- rbind(y, y_temp)
}
round(y, 3)
names(y) <- c("y1", "y2", "y3", "y4", "y5")
y_seg <- cbind(segment, y)
by(y_seg[, 2:6], y_seg$segment, colMeans)   #セグメントごとの平均

#セグメントごとにyを抽出
y1 <- subset(y_seg[, 2:6], y_seg[, 1] == 1)
y2 <- subset(y_seg[, 2:6], y_seg[, 1] == 2)
y3 <- subset(y_seg[, 2:6], y_seg[, 1] == 3)
y4 <- subset(y_seg[, 2:6], y_seg[, 1] == 4)

n/col
dim(y_seg)

##多変量正規分布から回帰構造のある混合正規乱数を発生させる
k = 4
n = 4000
y <- matrix(0, nrow=n, ncol=col+1)
cov <- list(Sigma1$covariance, Sigma2$covariance,Sigma3$covariance, Sigma4$covariance)
for(k in 1:4){
  yy <- subset(y_seg[, 2:6], y_seg[, 1] == k)
  cc <- as.matrix(cov[[k]])
  for(i in 1:1000){
    r <- (k-1)*1000 + i
    y[r, ] <- c(k, mvrnorm(n=1, as.matrix(yy[i, ]), cc))
  }
}
#データを見る
round(y, 3)
yy <- as.data.frame(y)
by(y_seg[, 2:6], y_seg[, 1], colMeans)   #元のデータのセグメントごとの平均
by(yy[, 2:6], yy[, 1], colMeans)   #乱数発生させたデータのセグメントごとの平均
by(yy[, 2:6], yy[, 1], cor) 

#結果をプロット
boxplot(yy[, 2:6])   #セグメントを無視すると…

#セグメント別に箱ひげ図を描画
par(mfrow=c(2, 3))
for(i in 2:6){
  boxplot(yy[, i] ~ yy[, 1], data=yy)
}
par(mfrow=c(1, 1))

plot(yy[, 2:6], col=yy[, 1])   #散布図

#アイリスデータの箱ひげ図
#par(mfrow=c(2, 2))
#for (response in c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"))
#     boxplot(iris[, response] ~ Species, data=iris, ylab=response)
#par(mfrow=c(1, 1))

#####セグメントごとに多変量回帰モデルを当てはめる####
##回帰モデルを当てはめる他のデータの準備
#反応変数
Yn <- y[, -1]   
Yn[1:10, ]; Yn[1001:1010, ]; Yn[2001:2010, ]; Yn[3001:3010, ]

iris
#説明変数
head(data_de)   #切片がついたデザイン行列
head(segment); tail(segment)   #セグメント
segment <- as.data.frame(segment)
data_seg <- cbind(segment, data_de)   #切片とセグメントがついたデザイン行列

#すべてがついた行列
yx <- cbind(segment, Yn, data_de)
names(yx)[2:6] <- c("y1", "y2", "y3", "y4", "y5")
head(yx)


##セグメントごとの多変量回帰モデルを実行
beta_seg <- list()   #回帰係数を入れるリスト
S_seg <- list()   #分散共分散行列を入れるリスト
for(i in 1:4){
  yx_s <- subset(yx, yx$segment==i)
  #回帰係数を推定
  beta_seg[[i]] <- solve(t(yx_s[, 7:16]) %*% as.matrix(yx_s[, 7:16])) %*% t(yx_s[, 7:16]) %*% as.matrix(yx_s[, 2:6])
  bs <- as.matrix(beta_seg[[i]])
  #分散共分散行列を推定
  S_seg[[i]] <- t(as.matrix(yx_s[, 2:6]) - as.matrix(yx_s[, 7:16])%*%bs) %*% 
                 (as.matrix(yx_s[, 2:6]) - as.matrix(yx_s[, 7:16])%*%bs) / nrow(yx_s)
}

beta_seg   #推定された回帰係数行列
S_seg   #推定された分散共分散行列

#適合度を確認
xs <- subset(yx, yx$segment==1)
ys <- as.matrix(xs[, 7:16]) %*% as.matrix(beta_seg[[1]])   #推定された反応変数
yt <- xs[, 2:6]   #元の反応変数
cbind(yt[, 1], ys[, 1])   
round(error <- yt - ys, 3)   #誤差

##関数を用いて推定
names(yx)
res2 <- by(yx, yx$segment, function(x) summary(lm(cbind(y1, y2, y3, y4, y5) ~ 
                                                        sex + age20 + age30 + age40 + age50 + 
                                                        work + stud + hwife + s_cnt, data=x)))

#セグメントごとの推定結果
res2$`1`
res2$`2`
res2$`3`
res2$`4`

####混合正規多変量回帰モデルの推定####