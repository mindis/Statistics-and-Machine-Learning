#####ベイジアン階層有限混合多項プロビットモデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(gtools)
library(MNP)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
library(lattice)

####多変量正規乱数を発生させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1

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

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
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


####データの発生####
#set.seed(8437)
##データの設定
hh <- 1000   #プレイヤー数
pt <- rpois(hh, 3)   #選択機会数
pt <- ifelse(pt==0, 1, pt)   #選択機会数が0なら1に置き換え
hhpt <- sum(pt)   #総サンプル数
member <- 10   #選択可能メンバー数
st <- 10   #基準メンバー
g <- 3   #セグメント数

##IDの設定
id <- rep(1:hh, pt)
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}

ID <- data.frame(no=1:hhpt, id=id, t=t)   #IDデータの結合


####説明変数の発生####
##離散選択モデルの説明変数の発生
#条件付きの説明変数の発生
X1.cont <- matrix(rnorm(hhpt*member*2, 0, 1), nrow=hhpt, ncol=member*2)

X1.bin <- matrix(0, nrow=hhpt, ncol=member*2)
for(i in 1:(member*2)){
  X1.bin[, i]  <- rbinom(hhpt, 1, runif(1, 0.35, 0.6))
}

#基準メンバーとの相対説明変数
X1.cont_r <- cbind(X1.cont[, 1:(member-1)] - X1.cont[, member], X1.cont[, (member+1):(2*member-1)] - X1.cont[, 2*member])
X1.bin_r <- cbind(X1.bin[, 1:(member-1)] - X1.bin[, member], X1.bin[, (member+1):(2*member-1)] - X1.bin[, 2*member])

#多項型の説明変数の発生
X2.cont <- c()
X2.bin <- c()

for(i in 1:hh){
  bin <- rbinom(1, 1, runif(1, 0.35, 0.7))
  X2.cont <- c(X2.cont, rep(rnorm(1, 0, 1), pt[i]))
  X2.bin <- c(X2.bin, rep(bin, pt[i]))
}

##ベクトル型説明変数にデータフォーマットを変更
#切片の設定
p <- c(1, rep(0, member-1))
pop <- matrix(p, nrow=hhpt*member, ncol=member-1, byrow=T)
POP <- pop[rowSums(pop) > 0, ]

#条件付き説明変数の設定
X1.cont_v <- cbind(as.numeric(t(X1.cont_r[, 1:(member-1)])), as.numeric(t(X1.cont_r[, member:(2*(member-1))])))
X1.bin_v <- cbind(as.numeric(t(X1.bin_r[, 1:(member-1)])), as.numeric(t(X1.bin_r[, member:(2*(member-1))])))

#多項型説明変数の設定
X2.v <- matrix(0, nrow=hhpt*(member-1), ncol=(member-1)*2)
for(i in 1:hhpt){
  r <- ((i-1)*(member-1)+1):((i-1)*(member-1)+member-1)
  X2.v[r, ] <- cbind(diag(X2.cont[i], member-1), diag(X2.bin[i], member-1))
}

#データの結合
X <- data.frame(mu=POP, c=X1.cont_v, b=X1.bin_v, m=X2.v)
XM <- as.matrix(X)
round(XM, 2)


##階層モデルの説明変数の発生
#連続変数の発生
cont <- 3
Z.cont <- matrix(rnorm(hh*cont, 0, 1), nrow=hh, ncol=cont)

#二値変数の発生
bin <- 2
Z.bin <- matrix(0, nrow=hh, ncol=bin)
for(i in 1:bin){
  Z.bin[, i] <- rbinom(hh, 1, runif(1, 0.4, 0.6))
}

#多値変数の発生
multi <- 4
p <- runif(multi, 0.25, 2)
Z.multi <- t(rmultinom(hh, 1, p))
Z.multi <- Z.multi[, -which.min(colSums(Z.multi))]

#データの結合
Zx <- data.frame(c=Z.cont, b=Z.bin, m=Z.multi)


##ベクトル型説明変数にデータフォーマットを変更
#切片の設定
p <- c(1, rep(0, g))
int <- matrix(p, nrow=hh*(g+1), ncol=g, byrow=T)
INT <- int[rowSums(int) > 0, -g]

#説明変数をベクトル型に変更
Zi.v <- matrix(0, nrow=nrow(Zx)*g, ncol=(cont+bin+multi-1)*2)

for(i in 1:hh){
  index <- ((i-1)*g+1):((i-1)*g+g)
  
  diag.x <- c()
  for(j in 1:(cont+bin)){
    diag.x <- cbind(diag.x, diag(Zx[i, j], g)[, -g])
  }
  
  diag.m <- matrix(0, nrow=g, ncol=(multi-1)*2)
  for(j in 1:(g-1)){
    r <- ((j-1)*g+1):((j-1)*g+g)
    diag.m[j, r] <- as.numeric(Zx[i, (cont+bin+1):ncol(Zx)])
  }
 Zi.v[index, ] <- cbind(diag.x, diag.m)
}

#データの結合
Zx.v <- cbind(INT, Zi.v)
round(Zx.v, 3)   #データの確認


####応答変数の発生####
##多項ロジットモデルよりセグメント割当を発生
#パラメータの設定
for(i in 1:1000){
  theta.z <- c(runif(g-1, -0.75, 0.75), runif(cont*(g-1), 0, 1), runif(bin*(g-1), -1.2, 1.2), runif((multi-1)*(g-1), -1.3, 1.3))
  
  #ロジットと確率の計算
  logit <- matrix(Zx.v %*% theta.z, nrow=hh, ncol=g, byrow=T)   #ロジットの設定
  Pr.z <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=g)   #確率の計算
  
  #多項分布よりセグメントを生成
  Z <- t(apply(Pr.z, 1, function(x) rmultinom(1, 1, x)))
  if(min(colMeans(Z)) > 0.25) {break}
}
Zi <- Z %*% 1:g
colSums(Z); colMeans(Z)

##IDごとにセグメントを割当
zi <- c()
z <- c()

for(i in 1:hh){
  zi <- c(zi, rep(Zi[i], pt[i]))
  z <- rbind(z, matrix(Z[i, ], nrow=pt[i], ncol=g, byrow=T))
}

ID <- data.frame(ID, z=zi)   #IDに結合

##IDをベクトル形式に変更
#セグメント割り当てを変更
zv1 <- c(); zv2 <- c(); idv <- c(); tv <- c()

for(i in 1:(member-1)){
  zv1 <- cbind(zv1, z)
  zv2 <- cbind(zv2, ID$z)
  idv <- cbind(idv, ID$id)
  tv <- cbind(tv, ID$t)
}

#ベクトル形式のデータに変更
z.vec <- matrix(as.numeric(t(zv1)), nrow=hhpt*(member-1), ncol=g, byrow=T)
zi.v <- as.numeric(t(zv2))
id.v <- as.numeric(t(idv))
time.v <- as.numeric(t(tv))

#IDを変更
ID.v <- data.frame(no=1:length(id.v), id=id.v, t=time.v, z=zi.v)
cbind(ID.v, z.vec)   #データを確認

##多項プロビットモデルより好きなメンバーを発生
#回帰パラメータの設定
beta0.z <- matrix(runif((member-1)*g, 0.2, 2.0), nrow=g, ncol=member-1, byrow=T)
beta1.z <- matrix(runif(g*2, 0, 0.9), nrow=g, ncol=2, byrow=T)
beta2.z <- matrix(runif(g*2, -0.9, 1.1), nrow=g, ncol=2, byrow=T)
beta3.z <- matrix(runif(g*(member-1), 0, 0.8), nrow=g, ncol=member-1, byrow=T)
beta4.z <- matrix(runif(g*(member-1), -0.9, 1.0), nrow=g, ncol=member-1, byrow=T)
beta.z <- t(cbind(beta0.z, beta1.z, beta2.z, beta3.z, beta4.z))   #データを結合
rownames(beta.z) <- 1:nrow(beta.z)

#分散供分散行列を設定
Cov <- corrM(member-1, -0.6, 0.75)   #分散共分散パラメータはセグメントで共通


#効用関数の設定
XM %*% beta.z 

