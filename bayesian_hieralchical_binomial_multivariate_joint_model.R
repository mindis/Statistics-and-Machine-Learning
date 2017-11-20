#####階層ベイズ二項-多変量同時プロビットモデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

#set.seed(3108)

####多変量正規分布を発生させる関数####
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
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
##データの設定
hh <- 1500   #ユーザー数
pt0 <- rpois(hh, 15.5)   #観測期間
pt <- ifelse(pt0 < 1, 1, pt0) 
hhpt <- sum(pt)   #総観測数
select <- 10   #観測ブランド数

##IDの設定
id <- rep(1:hh, pt)
time <- c()
for(i in 1:hh){time <- c(time, rep(1:pt[i]))}
ID <- data.frame(no=1:hhpt, id, time)

####説明変数の発生####
##来店したかどうかの説明変数
#チラシ有無
Flyer1 <- rbinom(hhpt, 1, 0.4)   
Flyer2 <- rbinom(hhpt, 1, 0.3)

#チラシ掲載商品数
Product1 <- rpois(hhpt, rgamma(hhpt, 25, 0.3)) * Flyer1
Product2 <- rpois(hhpt, rgamma(hhpt, 30, 0.25)) * Flyer2
Product1[Product1!=0] <- log(Product1[Product1!=0])
Product2[Product2!=0] <- log(Product2[Product2!=0])
freq_min <- min(Product1[Product1!=0])
Product1[Product1!=0] <- Product1[Product1!=0] - freq_min
Product2[Product2!=0] <- Product2[Product2!=0] - freq_min


#天気(降水量)
w0 <- exp(rnorm(hhpt, 0.35, 0.8))
Weather <- ifelse(w0 < 1, runif(1, 1, 2), w0) * rbinom(hhpt, 1, 0.3)
Weather[Weather!=0] <- log(Weather[Weather!=0]+0.1)

#休日かどうか
Holiday <- rbinom(hhpt, 1, 125/365)

##プロモータースコアの発生
#店舗1のスコア
Score01 <- rpois(hh, 3.4)
Score01[Score01 < 1] <- 1
Score01[Score01 > 5] <- 5
Score1 <- Score01[ID$id] - 3

#店舗2のスコア
Score02 <- rpois(hh, 3.0)
Score02[Score02 < 1] <- 1
Score02[Score02 > 5] <- 5
Score2 <- Score02[ID$id] - 3


#データの結合
X1 <- data.frame(切片=1, flyer=Flyer1, product=Product1, weather=Weather, holiday=Holiday, score=Score1, roy=0)
XM1 <- as.matrix(X1)


##ブランドを購買したかどうかの説明変数
##説明変数の発生
#通常価格の発生
PRICE <- matrix(runif(hhpt*select, 0.4, 1), nrow=hhpt, ncol=select)   

#ディスカウント率の発生
DISC <- matrix(runif(hhpt*select, 0, 0.5), nrow=hhpt, ncol=select)

#特別陳列の発生
DISP <- matrix(0, nrow=hhpt, ncol=select)
for(i in 1:select){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hhpt, 1, r)
}


#特別キャンペーンの発生
CAMP <- matrix(0, nrow=hhpt, ncol=select)
for(i in 1:select){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hhpt, 1, r)
}

#個別ブランドのスコア
Score_z0 <- matrix(0, nrow=hh, ncol=select)
for(i in 1:select){
  x <- runif(1, 2.3, 3.7)
  Score_z0[, i] <- rpois(hh, x)
}
Score_z0[Score_z0 < 1] <- 1
Score_z0[Score_z0 > 5] <- 5
Score_z <- Score_z0[ID$id, ] - 3


##データのベクトル変換
#IDのベクトル化
id_vec <- rep(1:hh, select*pt)
time_vec <- as.numeric(t(matrix(time, nrow=hhpt, ncol=select)))
select_vec <- rep(1:select, hhpt)
select_no <- rep(1:hhpt, rep(select, hhpt)) 
ID_vec0 <- data.frame(no=1:(hhpt*select), num=select_no, id=id_vec, time=time_vec, select=select_vec)

#説明変数のベクトル化
#回帰係数が全ブランドで共通の説明変数
DISP.vec <- as.numeric(t(DISP))
CAMP.vec <- as.numeric(t(CAMP))
SCORE.vec <- as.numeric(t(Score_z))

#回帰係数がブランドで異なる説明変数
BP.vec <- matrix(diag(select), nrow=hhpt*select, ncol=select, byrow=T)
PRICE.vec <- matrix(0, nrow=hhpt*select, ncol=select)
DISC.vec <- matrix(0, nrow=hhpt*select, ncol=select)

for(i in 1:hhpt){
  print(i)
  r <- which(ID_vec0$num==i)
  PRICE.vec[r, ] <- diag(PRICE[i, ])
  DISC.vec[r, ] <- diag(DISC[i, ])
}

#データを結合
X02<- data.frame(bp=BP.vec, price=PRICE.vec, disc=DISC.vec, disp=DISP.vec, camp=CAMP.vec, score=SCORE.vec)
XM02 <- as.matrix(X02)


####階層モデルの説明変数の発生####
##デモグラフィック変数の発生
#連続変数の発生
cont <- 2
x.cont <- matrix(rnorm(cont*hh, 0, 1), nrow=hh, ncol=cont)

#二値変数の発生
bin <- 3
x.bin <- matrix(0, nrow=hh, ncol=bin)
for(i in 1:bin){
  p <- runif(1, 0.3, 0.5)
  x.bin[, i] <- rbinom(hh, 1, p) 
}

#多値変数の発生
multi <- 3
p <- c(0.4, 0.3, 0.3)
x.multi <- t(rmultinom(hh, 1, p))[, -multi]


#データの結合
Z <- data.frame(切片=1, cont=x.cont, bin=x.bin, multi=x.multi)
ZX <- as.matrix(Z)


####応答変数の発生####
##パラメータの設定
#変量効果の分散共分散行列の設定
tau1 <- diag(runif(ncol(XM1), 0.05, 0.25))
tau2 <- diag(runif(ncol(XM02), 0.05, 0.25))

#パラメータ数を設定
par0 <- ncol(ZX)
par1 <- ncol(XM1)
par2 <- ncol(XM02)

##ブランド購買有無の応答変数を発生
for(i in 1:1000){
  print(i)
  
  #ブランド購買の分散共分散行列の発生
  Cor0 <- corrM(select, -0.6, 0.9, 0.01, 0.1)
  
  #階層モデルのパラメータを設定
  delta0 <- cbind(matrix(runif(par0*select, -0.55, 0.55), nrow=par0, ncol=select),
                  matrix(runif(par0*select, -0.6, 0.2), nrow=par0, ncol=select),
                  matrix(runif(par0*select, -0.2, 0.5), nrow=par0, ncol=select),
                  runif(par0, -0.25, 0.5), runif(par0, -0.25, 0.5), runif(par0, -0.3, 0.3))
  
  #ブランド購買有無のパラメータを階層モデルから生成
  beta0 <- ZX %*% delta0 + mvrnorm(hh, rep(0, par2), tau2)
  
  #ブランド購買有無の応答変数を発生
  mu2 <- matrix(rowSums(XM02 * beta0[ID_vec0$id, ]), nrow=hhpt, ncol=select, byrow=T)   #平均構造
  U2 <- t(apply(mu2, 1, function(x) mvrnorm(1, x, Cor0)))   #潜在効用
  y02 <- ifelse(U2 > 0, 1, 0)
  
  print(round(c(min(colMeans(y02)), max(colMeans(y02))), 3))
  if(min(colMeans(y02)) > 0.2 & max(colMeans(y02)) < 0.7) break
}

##来店有無の応答変数を発生
X1$roy <- scale(rowSums(mu2))
XM1[, ncol(XM1)] <- scale(rowSums(mu2))

for(i in 1:1000){
  print(i)
  
  #階層モデルのパラメータを設定
  gamma0 <- cbind(runif(par0, -0.6, 0.7), runif(par0, -0.3, 0.6), runif(par0, -0.3, 0.6), runif(par0, -0.7, 0.2),
                  runif(par0, -0.3, 0.6), runif(par0, -0.5, 0.5), runif(par0, -0.1, 0.25))
  
  #来店有無のパラメータを階層モデルから生成
  alpha0 <- ZX %*% gamma0 + rmvnorm(hh, rep(0, par1), tau1)
  
  #来店有無の応答変数を発生
  mu1 <- rowSums(XM1 * alpha0[ID$id, ])   #潜在効用
  U1 <- mu1 + rnorm(hhpt, 0, 1)   #潜在効用
  y1 <- ifelse(U1 > 0, 1, 0)
  
  if(mean(y1) > 0.4 & mean(y1) < 0.6) break
}

##来店がないサンプルは欠損させる
#インデックスを作成
index_y <- which(y1==0)

index_list <- list()
for(i in 1:length(index_y)){
  index_list[[i]] <- which(ID_vec0$num==index_y[i])
}
index_zeros <- unlist(index_list)

#ブランド購買有無を欠損させる
ID_vec <- ID_vec0[-index_zeros, ]
n <- length(unique(ID_vec$id))
XM2 <- XM02[-index_zeros, ]
PRICE_z <- PRICE[-index_y, ]
DISC_z <- DISC[-index_y, ]
DISP_z <- DISP[-index_y, ]
CAMP_z <- CAMP[-index_y, ]
Score_z <- Score1[-index_y]
y2 <- y02[-index_y, ]
ZX1 <- ZX
ZX2 <- ZX[unique(ID_vec$id), ]


####マルコフ連鎖モンテカルロ法で階層ベイズ二項-多変量同時プロビットモデルを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##多変量正規分布の条件付き期待値と条件付き分散を計算する関数
cdMVN <- function(mean, Cov, dependent, U){
  
  #分散共分散行列のブロック行列を定義
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent, drop=FALSE]
  Cov21 <- Cov[-dependent, dependent, drop=FALSE]
  Cov22 <- Cov[-dependent, -dependent]
  
  #条件付き分散と条件付き平均を計算
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mean[, dependent] + t(CDinv %*% t(U[, -dependent] - mean[, -dependent]))   #条件付き平均を計算
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #条件付き分散を計算
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

#プロビットモデルの対数尤度の定義
probit_LL <- function(x, Y, X){
  #パラメータの設定
  b0 <- x[1]
  b1 <- x[2:(ncol(X)+1)]
  
  #効用関数の定義
  U <- b0 + as.matrix(X) %*% b1
  
  #対数尤度を計算
  Pr <- pnorm(U)   #確率の計算
  LLi <- Y*log(Pr) + (1-Y)*log(1-Pr)
  LL <- sum(LLi)
  return(LL)
}


##アルゴリズムの設定
R <- 20000 
keep <- 2
disp <- 10
iter <- 0
sbeta <- 1.5

##事前分布の設定
#多変量プロビットモデルの事前分布
nu <- select   #逆ウィシャーと分布の自由度
V <- nu1 * diag(rep(1, select))   #逆ウィシャート分布のパラメータ

#階層モデルの事前分布
#来店有無の階層モデルのパラメータ
Deltabar <- matrix(rep(0, par0*par1), nrow=ncol(ZX), ncol=par1)   #階層モデルの回帰係数の事前分布の分散
ADelta <- 0.01 * diag(rep(1, par0))   #階層モデルの回帰係数の事前分布の分散
nu1 <- par1   #逆ウィシャート分布の自由度
V1 <- nu1 * diag(rep(1, par1)) #逆ウィシャート分布のパラメータ

#ブランド購買の階層モデルのパラメータ
Deltabar <- matrix(rep(0, par0*par2), nrow=ncol(ZX), ncol=par2)   #階層モデルの回帰係数の事前分布の分散
ADelta <- 0.01 * diag(rep(1, par0))   #階層モデルの回帰係数の事前分布の分散
nu2 <- par2   #逆ウィシャート分布の自由度
V2 <- nu2 * diag(rep(1, par2)) #逆ウィシャート分布のパラメータ


##サンプリング結果の保存用配列
ALPHA <- array(0, dim=c(hh, par1, R/keep))
BETA <- array(0, dim=c(n, par2, R/keep))
COR <- array(0, dim=c(select, select, R/keep))
THETA1 <- matrix(0, nrow=R/keep, ncol=par0*par1)
THETA2 <- matrix(0, nrow=R/keep, ncol=par0*par2)
COV <- matrix(0, nrow=R/keep, ncol=par1+par2)
gc(); gc()


##初期値の設定
#来店有無のプロビットモデルを推定
#初期パラメータとデータを設定
x <- rep(0, par1)
X <- XM1[, -1]

#準ニュートン法で対数尤度最大化
res <- optim(x, probit_LL, Y=y1, X=X, method="BFGS", hessian=FALSE, 
             control=list(fnscale=-1))
oldalpha <- matrix(res$par, nrow=hh, ncol=par1, byrow=T) + mvrnorm(hh, rep(0, par1), diag(0.2, par1))


#ブランド購買ごとに独立にプロビットモデルを推定
beta.init0 <- list()

for(j in 1:select){
  print(j)
  
  #初期パラメータとデータを設定
  X <- cbind(PRICE_z[, j], DISC_z[, j], DISP_z[, j], CAMP_z[, j], Score_z)
  x <- rep(0, ncol(X) + 1)

  #準ニュートン法で対数尤度最大化
  res <- optim(x, probit_LL, Y=y2[, j], X=X, method="BFGS", hessian=FALSE, 
               control=list(fnscale=-1))
  beta.init0[[j]] <- res$par
}

#個人別のプロビットモデルのパラメータを設定
beta.init <- do.call(rbind, beta.init0) 
betad <- c(as.numeric(t(beta.init[, 1:3])), colMeans(beta.init[, 4:ncol(beta.init)])) 
oldbeta <- matrix(betad, nrow=n, ncol=par2, byrow=T) + mvrnorm(n, rep(0, par2), diag(0.15, par2))

#分散共分散行列の初期値
oldcor <- diag(select)


##階層モデルのパラメータを設定
#来店有無の階層モデルのパラメータ
oldtheta1 <- solve(t(ZX1) %*% ZX1) %*% t(ZX1) %*% oldalpha   #回帰係数
oldcov1 <- var(oldalpha - ZX1 %*% oldtheta1)   #分散共分散行列

#ブランド購買の階層モデルのパラメータ
oldtheta2 <- solve(t(ZX2) %*% ZX2) %*% t(ZX2) %*% oldbeta   #回帰係数
oldcov2 <- var(oldbeta - ZX2 %*% oldtheta2)   #分散共分散行列


##切断領域を定義
a1 <- ifelse(y1==0, -100, 0)
b1 <- ifelse(y1==1, 100, 0)
a2 <- ifelse(y2==0, -100, 0)
b2 <- ifelse(y2==1, 100, 0)


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
##来店有無の潜在効用を発生
util_mu1 <- rowSums(XM1 * oldalpha[ID$id, ])   #潜在効用の平均
util1 <- rtnorm(util_mu1, 1, a1, b1)

##ブランド購買の潜在効用を発生




