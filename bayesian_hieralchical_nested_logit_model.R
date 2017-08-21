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

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
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
hh <- 500
member <- 9   #メンバー数
c_num <- 8   #衣装パターン
hhpt <- hh*member*c_num   #全変数数


####説明変数の発生####
##個体内説明変数の発生
#メンバーの説明変数の設定
Mus <- matrix(as.numeric(table(1:hhpt, rep(rep(1:member, rep(c_num, member)), hh))), nrow=hhpt, ncol=member)
colnames(Mus) <- c("hono", "koto", "umi", "rin", "hana", "maki", "nico", "eri", "nozo")
Mus <- Mus[, -ncol(Mus)]

#衣装の説明変数の設定
ct <- c(1, rep(0, c_num))
cloth <- matrix(ct, nrow=hh*member*(c_num+1), ncol=c_num, byrow=T)
CLOTH <- subset(cloth, rowSums(cloth) > 0)[, -c_num]
colnames(CLOTH) <- c("A", "B", "C", "D", "F", "G", "H")

#カード種別
type <- 3   #種類数
card <- t(rmultinom(hhpt, 1, c(2/c_num, 2/c_num, (c_num-4)/c_num)))
colnames(card) <- c("UR", "SSR", "SR")
CARD <- card[, -type]

#プロモーション接触数
Prom <- scale(rpois(hhpt, 5))

#データの結合
X <- data.frame(Mus, CLOTH, CARD, Prom)
XM1 <- as.matrix(X)
XM2 <- XM1[, -ncol(XM1)]

#パラメータ数
k1 <- 4 + ncol(XM1) + ncol(XM2) - (member-1) - (c_num-1)

##IDの設定
id <- rep(1:hh, rep(member*c_num, hh))
pt <- rep(1:(member*c_num), hh)
ID <- data.frame(no=1:hhpt, id=id, pt=pt)


##個体間説明変数の発生
#連続変数
cont.h <- 3
Z.cont <- matrix(runif(hh*cont.h, 0, 1), nrow=hh, ncol=cont.h) 

#二値変数
bin.h <- 3
Z.bin <- matrix(0, nrow=hh, ncol=bin.h)
for(i in 1:bin.h){
  p.bin <- runif(1, 0.2, 0.8)
  Z.bin[, i] <- rbinom(hh, 1, ph.bin)  
}

#多値変数
multi.h <- 4
p.multi <- runif(multi.h)
Z.multi <- t(rmultinom(hh, 1, ph.multi))
freq.min <- which.min(colSums(Xh.multi))
Z.multi <- Xh.multi[, -freq.min]

#データの結合
Z <- data.frame(cont=Xh.cont, bin=Xh.bin, multi=Xh.multi)
ZM <- as.matrix(Z)
k2 <- ncol(ZM)

####応答変数の発生####
##個体間パラメータの設定
#個体間分散共分散行列を設定
corM <- corrM(col=k1, lower=-0.55, upper=0.9, eigen_lower=0.025, eigen_upper=0.35)   #相関行列を作成
Sigma <- covmatrix(col=k1, corM=corM, lower=1, upper=1)   #分散共分散行列
Cov <- Sigma$covariance


#個体間回帰パラメータの設定
theta0 <- c(runif(1, -1.1, -0.55), runif(1, -1.4, -1.1), runif(member-1, -0.2, 0.85), runif(c_num-1, -0.5, 0.7), 
            runif(1, -0.6, -0.4), runif(1, -0.5, -0.3), runif(1, 0.05, 0.15), runif(1, -0.2, 0.5), 
            runif(1, -0.4, 0.3), runif(1, -0.8, -0.5), runif(1, -0.55, -0.3))

theta1 <- matrix(c(runif(2*k2, -0.4, 0.5), runif((member-1)*k2, -0.4, 0.55), runif((c_num-1)*k2, -0.4, 0.5), 
            runif(k2, -0.4, 0.4), runif(k2, -0.4, 0.4), runif(k2, -0.1, 0.15), runif(k2, -0.3, 0.4), runif(k2, -0.2, 0.3),
            runif(k2, -0.4, 0.4), runif(k2, -0.35, 0.35)), nrow=k2, ncol=k1, byrow=T)

#パラメータの結合
theta <- rbind(theta0, theta1)


##個体内回帰パラメータの設定
beta <- cbind(1, ZM) %*% theta + mvrnorm(hh, rep(0, k1), Cov)
beta1 <- beta[, c(1, 3:(ncol(beta)-2))]
beta2 <- beta[, c(2:(member-1 + c_num-1 + 2), ncol(beta)-1, ncol(beta))]

#ログサム変数のパラメータを0～1に収まるように変換しておく
beta1[, (ncol(beta1)-1):ncol(beta1)] <- exp(beta1[, (ncol(beta1)-1):ncol(beta1)])/
                                                  (1+exp(beta1[, (ncol(beta1)-1):ncol(beta1)]))


##応答変数の発生
#ロジットとログサム変数を定義
logit1.list <- list()
logit2.list <- list()
logsum.list <- list()

#全ユーザーの個人ごとのログサム変数とロジットを計算
for(i in 1:hh){
  print(i)
  #ログサム変数の定義
  logsum.list[[i]] <- log(1 + exp(cbind(1, XM2)[ID$id==i, 1:(ncol(XM2)-2)] %*% beta2[i, 1:(ncol(XM2)-2)]))
  
  #ロジットの定義
  logit1.list[[i]] <- cbind(1, XM1, matrix(logsum.list[[i]], nrow=hhpt, ncol=type-1)*CARD)[ID$id==i, ] %*% beta1[i, ]
  logit2.list[[i]] <- cbind(1, XM2)[ID$id==i, ] %*% beta2[i, ]
}

#リストを数値型に変更
logsum <- unlist(logsum.list)
logit1 <- unlist(logit1.list)
logit2 <- unlist(logit2.list)

#ネステッドロジットモデルににより確率を計算し、ベルヌーイ分布より応答変数を発生
#カードを持っているかどうか
Pr1 <- exp(logit1)/(1+exp(logit1))
y1 <- rbinom(hhpt, 1, Pr1)

#カードが覚醒しているかどうか
Pr2 <- exp(logit2)/(1+exp(logit2))
y2 <- rbinom(hhpt, 1, Pr2)
y2[y1==0] <- NA   #カードを持っている場合のみ覚醒有無を定義する


####マルコフ連鎖モンテカルロ法で階層ベイズネステッドロジットモデルを推定####
##ネステッドロジットモデルの対数尤度関数
loglike <- function(x, y1, y2, XM1, XM2, CARD, type, member, c_num, hhpt){

  #パラメータの設定
  beta1 <- x[, c(1, 3:(length(x)-2))]
  beta2 <- x[, c(2:(member-1 + c_num-1 + 2), length(x)-1, length(x))]
  
  #ログサム変数のパラメータを0～1に収まるように変換しておく
  beta1[, (length(beta1)-1):length(beta1)] <- exp(beta1[, (ncol(beta1)-1):ncol(beta1)])/
                                                          (1+exp(beta1[, (ncol(beta1)-1):ncol(beta1)]))
  
  #ロジットとログサム変数を定義
  logit2 <- cbind(1, XM2) %*% beta2   #覚醒有無のロジット
  logsum <- log(1 + exp(cbind(1, XM2)[, 1:(ncol(XM2)-2)] %*% beta2[, 1:(ncol(XM2)-2)]))   #ログサム変数
  logit1 <- cbind(1, XM1, matrix(logsum, nrow=hhpt, ncol=type-1)*CARD) %*% beta1   #カード所有有無のロジット
  
  #対数尤度を定義する
  #カード所有有無の対数尤度
  Pr.l <- exp(logit1) / (1 + exp(logit1))
  LLs.l <- y1*log(Pr.l) + (1-y1)*log(1-Pr.l)  
  LL.l <- sum(LLs.l)
  
  #覚醒有無の対数尤度
  Pr.b <- exp(logit2[y1==1]) / (1 + exp(logit2[y1==1]))
  LLs.b <- y2[y1==1]*log(Pr.b) + (1-y2[y1==1])*log(1-Pr.b)  
  LL.b <- sum(LLs.b)
  
  #対数尤度を合計
  LL <- LL.l + LL.b
  return(LL)
}


