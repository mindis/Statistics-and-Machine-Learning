#####相関構造のある階層ベイズRFモデル#####
library(MASS)
library(nlme)
library(glmm)
library(survival)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(431278)

####データの発生####
##データの設定
hh <- 1000   #サンプル数
pt <- rpois(hh, 10.0)
pt <- ifelse(pt==0, 1, pt)
hhpt <- sum(pt)
dt <- 100

##IDの設定
id <- rep(1:hh, pt)
time <- c()
for(i in 1:hh){time <- c(time, 1:pt[i])}
ID <- data.frame(no=1:hhpt, id=id, time=time)

####説明変数の発生####
##時点で変化しない説明変数
k1 <- 5
cont1 <- 2; bin1 <- 3
X1 <- matrix(0, nrow=hhpt, ncol=k1)
for(i in 1:hh){
  X1[ID$id==i, 1:cont1] <- matrix(rnorm(2, 0, 1), nrow=pt[i], ncol=2, byrow=T)
  for(j in 1:bin1){
    X1[ID$id==i, (cont1+j)] <- rep(rbinom(1, 1, runif(1, 0.4, 0.6)), pt[i])
  }
}

##時変の説明変数
Price <- runif(hhpt, 0.6, 1.0)   #旧作の平均値引率
Disc <- runif(hhpt, 0.1, 0.6)   #旧作の値引きゲームの割合
Prom <- rbinom(hhpt, 1, 0.4)   #プロモーション有無

#新作ゲーム本数とジャンルの発生
g <- 6
Genre <- matrix(0, nrow=hhpt, ncol=g)
pr_gen <- runif(g, 0.3, 1.0)

#サンプルごとに多項分布からゲームジャンルを発生させる
for(i in 1:hhpt){
  new <- rpois(1, 6.5)
  if(new==0){
    next
  } else {
    Genre[i, ] <- t(rmultinom(1, new, pr_gen))
  }
}

##データを結合
X <- data.frame(1, X1, Price, Disc, Prom, Genre)
colnames(X) <- c("切片", "cont1", "cont2", "bin1", "bin2", "bin3", "Price", "Disc", "Prom", "genre1", "genre2", "genre3", 
                 "genre4", "genre5", "genre6")
XM <- as.matrix(X)


####応答変数の発生####
for(i in 1:1000){
  ##パラメータの設定
  #生存モデルのパラメータの設定
  alpha0 <- runif(1, 0.8, 1.8)
  theta0 <- c(runif(1, 0, 5.0), runif(cont1, 0, 0.4), runif(bin1, -0.4, 0.8), runif(1, 0.4, 1.0), runif(1, -0.8, -0.3),
              runif(1, -0.7, -0.2), runif(g, -0.4, 0))
  
  #頻度モデルのパラメータの設定
  beta0 <- c(runif(1, 0, 0.3), runif(cont1, 0, 0.2), runif(bin1, -0.2, 0.1), runif(1, -0.4, -0.2), runif(1, 0, 0.2),
             runif(1, 0, 0.2), runif(g, 0, 0.15)) 
  
  ##生存モデルと頻度モデルのスケールパラメータを多変量正規分布から発生
  #分散共分散パラメータを設定
  Cov0 <- matrix(c(0.5, -0.35, -0.35, 0.7), nrow=2, ncol=2)
  
  #平均構造を設定
  scale0 <- XM %*% theta0 
  lambda0 <- XM %*% beta0
  
  #多変量対数正規分布からパラメータを発生
  Mu <- exp(t(apply(cbind(scale0, lambda0), 1, function(x) mvrnorm(1, x, Cov0))))
  scale_t <- Mu[, 1]
  lambda_t <- Mu[, 2]
  
  #ワイブル分布とポアソン分布より購買間隔と購買数を発生
  y1 <- rweibull(hhpt, alpha0, scale_t)
  y2 <- rpois(hhpt, lambda_t)
  
  print(round(c(min(y1), max(y1), min(y2), max(y2)), 3))
  if(min(y1) > 0.025 & max(y2) < 50) break
}
Y <- cbind(y1, y2)

#発生させた応答変数を可視化
hist(y1[y1 < 100], xlab="購買間隔", main="ゲーム店への訪問間隔", col="grey")
hist(y2, xlab="購買数", main="ゲームの購買数", col="grey")

####打ち切りの設定####
##ユーザーごとにT = 100まで観測
#変数の格納用リスト
ID.list <- list()
y1.list <- list()
y2.list <- list()
X.list <- list()
z.list <- list()

#個人ごとに打ち切り変数を設定
for(i in 1:hh){
  print(i)
  y1_ind <- Y[ID$id==i, 1]
  y2_ind <- Y[ID$id==i, 2]
  z <- rep(0, length(y1_ind))
  c_sum <- cumsum(y1_ind)
  
  #累積時間が100以上のイベントは打ち切り
  index1 <- subset(1:length(c_sum), c_sum <= 100)
  
  if(max(c_sum) <= 100){
    index2 <- index1
  } else {
    index2 <- c(index1, length(index1)+1)
  }
  
  #応答変数の打ち切りを設定
  if(max(c_sum) > dt & length(index1) > 0){
    print(1)
    y_vec <- c(y1_ind[index1], dt-c_sum[length(index1)])
    z[length(y_vec)] <- 1
  } else if(max(c_sum) > dt & length(index1)==0) {
    print(2)
    y_vec <- dt
    z <- 1
  } else {
    print(3)
    y_vec <- y1_ind[index2]
  }

  
  #打ち切られた変数を格納
  y1.list[[i]] <- y_vec[index2]
  y2.list[[i]] <- y2_ind[index2]
  ID.list[[i]] <- ID[ID$id==i, ][index2, ]
  X.list[[i]] <- X[ID$id==i, ][index2, ]
  z.list[[i]] <- z[index2]
}

#リストをベクトルあるいは行列化
y1 <- unlist(y1.list)
y2 <- unlist(y2.list) 
ID <- do.call(rbind, ID.list)
X <- do.call(rbind, X.list)
XM <- as.matrix(X)
z <- unlist(z.list)
hhpt <- nrow(ID)

#データの確認と可視化
round(cbind(ID, y1, y2, z, X), 3)
hist(y1, xlab="購買間隔", main="ゲーム店への訪問間隔", col="grey")
hist(y2, xlab="購買数", main="ゲームの購買数", col="grey")


####マルコフ連鎖モンテカルロ法で階層ベイズRFモデルを推定####
##RFモデルの対数尤度関数の設定
llike <- function(alpha, scale, lambda, y1, y2, X, z){
  #パラメータの設定
  scale1 <- exp(scale)
  lambda1 <- exp(lambda)
  
  #対数尤度を計算
  LL1 <- sum(z*(log(scale1)+log(alpha)+(alpha-1)*log(y1)) - scale1*y1^alpha)   #ワイブルモデルの対数尤度
  LL2 <- y2*log(lambda1)-lambda1 - lfactorial(y2)
  LL <- sum(LL1) + sum(LL2)
  return(LL)
}

##アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
len <- nrow(X)
par <- ncol(X)

##事前分布の設定
#形状パラメータの事前分布
alpha_mu <- 0
alpha_sigma <- 2.5

#階層モデルの事前分布
Deltabar <- matrix(0, nrow=ncol(X), ncol=2)
Adelta <- 0.01 * diag(2)
nu <- 5
V <- nu * diag(2)

##サンプリング結果の保存用配列
THETA <- matrix(0, nrow=R/keep, ncol=2)
ALPHA <- rep(0, R/keep)
BETA <- matrix(0, nrow=R/keep, ncol=ncol(X)*2)
SIGMA <- matrix(0, nrow=R/keep, ncol=2^2)

##棄却率と対数尤度の保存用配列
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##初期値の設定
#階層モデルの回帰パラメータの設定
theta_old <- c(runif(1, 2.5, 5.5), runif(cont1, 0, 0.5), runif(bin1, -0.6, 1.0), runif(1, 0, 1.0), runif(1, -1.0, 0),
              runif(1, -1.0, 0), runif(g, -0.7, 0))

#頻度モデルのパラメータの設定
beta_old <- c(runif(1, 0, 0.5), runif(cont1, 0, 0.5), runif(bin1, -0.3, 0.3), runif(1, -0.5, 0), runif(1, 0, 0.4),
           runif(1, 0, 0.4), runif(g, 0, 0.3)) 

oldbeta <- cbind(theta_old, beta_old)   #パラメータを結合

#階層モデルの分散共分散行列の設定
oldcov <- matrix(c(0.7, -0.3, -0.3, 0.7), nrow=2, ncol=2)   

#ワイブルモデルの形状パラメータの設定
oldalpha <- runif(1, 0.9, 1.5)

#ワイブルモデルおよびポアソンモデルの個人別パラメータを設定
oldtheta <- XM %*% oldbeta + mvrnorm(hhpt, rep(0, 2), oldcov)
oldscale <- oldtheta[, 1]
oldlambda <- oldtheta[, 2]


####マルコフ連鎖モンテカルロ法で階層ベイズRFモデルを推定####



