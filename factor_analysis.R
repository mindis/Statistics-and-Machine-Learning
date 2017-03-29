#####因子分析#####
library(polycor)
library(psych)
library(MASS)
library(reshape2)
library(plyr)
####測定方程式からデータを発生####
#set.seed(489)
n <- 1000   #サンプル数
k <- 7   #変数数
ff <- 3   #因子数
mu <- rnorm(7, 5, 5)   #平均ベクトル
A <- matrix(runif(k*ff, 0, 0.6), k, ff)   #因子負荷行列 
F <- matrix(rnorm(n*ff, 0, 1), n, ff)   #共通因子行列

##測定方程式を設定して、観測変数を作成
MU <- matrix(mu, n, k, byrow=T)
X <- F %*% t(A) + matrix(rnorm(n*k, 0, 0.4), n, k, byrow=T)   #観測変数
A

#発生させたデータの要約
summary(X)
colMeans(X)
cor(X)
plot(as.data.frame(X), col=4)


####多変量正規分布からのデータ発生####
val <- 7   #説明変数の数
n <- 1000   #サンプル数

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

#群ごとの相関行列を作成(群ですべて同じ)
corM <- corrM(col=val, lower=-0.3, upper=0.8)
eigen(corM)

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
    diag(c) <- m   #対角行列を元の分散に戻す
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

#分散共分散行列を作成(群ですべて同じ)
Sigma <- covmatrix(col=val, corM=corM, lower=1, upper=3)

#群ごとの変数の平均を作成
mu <- c(rnorm(val, 3, 2))

##多変量正規分布からの乱数を発生させる
val; n
X <- mvrnorm(n=n, mu, Sigma$covariance)

##データを要約
round(colMeans(X), 2)   #変数ごとの平均
round(var(X), 2)   #分散共分散行列
round(cor(X), 2)   #相関行列
summary(X)   #要約

##散布図行列の作成
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#散布図行列
pairs(as.data.frame(X[, 1:7]), panel=panel.smooth, bg="lightblue", diag.panel=panel.hist,
      upper.panel=panel.cor)



##最尤推定による因子負荷行列の推定
##EMアルゴリズム
#初期値の設定
R <- var(X)
Aold <- matrix(rep(0.5, k*ff), k, ff)
Aold[1, 2:3] <- 0
Aold[2, 3] <- 0
D2old <- diag(R - Aold %*% t(Aold))
D2old <- diag(D2old)
I <- diag(1, ff)

#EMアルゴリズムの設定
max.iter <- 1000   #最大繰り返し数
iter <- 1   #繰り返しカウンター
tol <- 10^(-1)   #推定値の変化の許容度
S.zz <- R   #S.zzの条件付き期待値はデータの分散共分散行列

##EMアルゴリズムによる推定
while(iter < max.iter){
  #Eステップ
  Sigma <- Aold %*% t(Aold) + D2old
  delta <- t(Aold) %*% solve(Sigma)
  S.zf <- S.zz %*% t(delta)
  S.ff <- delta %*% S.zz %*% t(delta) + (I - delta %*% Aold)
  #Mステップ
  Anew <- S.zf %*% solve(S.ff)
  D2new <- diag(S.zz - S.zf %*% solve(S.ff) %*% t(S.zf))
  #推定値の変化
  diff <- max(abs(Anew - Aold), abs(D2new - diag(D2old)))
  if(diff < tol) break;
  Aold <- Anew
  D2old <- diag(D2new)
  iter <- iter+1
  print(diff)
}

#結果を確認
delta <- round(1-D2new, 2)
Anew
promax(Anew)

##関数を用いて推定
fa(X, nfactors=3, rotate="promax", fm="ml", cor="cov")
promax(Anew)
help(fa)
