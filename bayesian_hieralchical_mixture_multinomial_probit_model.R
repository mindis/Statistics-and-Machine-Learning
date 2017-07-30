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
hh <- 500   #プレイヤー数
pt <- rpois(hh, 5)   #選択機会数
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
id_r <- matrix(1:(hhpt*(member-1)), nrow=hhpt, ncol=member-1, byrow=T)


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
X <- data.frame(mu=POP, c=X1.cont_v, b=X1.bin_v, m=X2.v[, 1:(member-1)])
XM <- as.matrix(X)
round(XM, 2)


##階層モデルの説明変数の発生
#連続変数の発生
cont <- 3
Z.cont <- matrix(rnorm(hh*cont, 0, 1), nrow=hh, ncol=cont)

#二値変数の発生
bin <- 3
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
Zx <- data.frame(c=Z.cont, b=Z.bin)


##ベクトル型説明変数にデータフォーマットを変更
#切片の設定
p <- c(1, rep(0, g))
int <- matrix(p, nrow=hh*(g+1), ncol=g, byrow=T)
INT <- int[rowSums(int) > 0, -g]

#説明変数をベクトル型に変更
Zi.v <- matrix(0, nrow=nrow(Zx)*g, ncol=(cont+bin)*2)

for(i in 1:hh){
  index <- ((i-1)*g+1):((i-1)*g+g)
  
  diag.x <- c()
  for(j in 1:(cont+bin)){
    diag.x <- cbind(diag.x, diag(Zx[i, j], g)[, -g])
  }
  
  #  diag.m <- matrix(0, nrow=g, ncol=(multi-1)*2)
  #  for(j in 1:(g-1)){
  #    r <- ((j-1)*g+1):((j-1)*g+g)
  #    diag.m[j, r] <- as.numeric(Zx[i, (cont+bin+1):ncol(Zx)])
  # }
  Zi.v[index, ] <- cbind(diag.x)
}

#データの結合
Zx.v <- cbind(INT, Zi.v)
round(Zx.v, 3)   #データの確認


####応答変数の発生####
##多項ロジットモデルよりセグメント割当を発生
#パラメータの設定
for(i in 1:1000){
  theta.z <- c(runif(g-1, -0.75, 0.75), runif(cont*(g-1), 0, 1), runif(bin*(g-1), -1.2, 1.2))

  #ロジットと確率の計算
  logit <- matrix(Zx.v %*% theta.z, nrow=hh, ncol=g, byrow=T)   #ロジットの設定
  Pr.z <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=g)   #確率の計算
  
  #多項分布よりセグメントを生成
  Z <- t(apply(Pr.z, 1, function(x) rmultinom(1, 1, x)))
  if(min(colMeans(Z)) > 0.25) {break}
}
Zi <- Z %*% 1:g
r_rate <- colSums(Z)/sum(Z)   #混合率
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
beta0.z <- rbind(c(3.5, 2.4, 2.7, 0.7, 1.0, 0.8, -0.5, -0.2, 1.4),
                 c(1.6, 1.0, 0.8, 2.8, 4.0, 2.8, 1.2, -0.7, 0.6),
                 c(0.7, -0.3, -0.4, 1.2, 1.5,0.6, 2.6, 3.0, 3.3))
beta1.z <- matrix(runif(g*2, 0, 1.2), nrow=g, ncol=2, byrow=T)
beta2.z <- matrix(runif(g*2, -1.2, 1.2), nrow=g, ncol=2, byrow=T)
beta3.z <- matrix(runif(g*(member-1), 0, 1.2), nrow=g, ncol=member-1, byrow=T)
beta4.z <- matrix(runif(g*(member-1), -1.2, 1.3), nrow=g, ncol=member-1, byrow=T)
beta.z <- t(cbind(beta0.z, beta1.z, beta2.z, beta3.z))   #データを結合
rownames(beta.z) <- 1:nrow(beta.z)

#分散供分散行列を設定
Cov <- corrM(member-1, -0.55, 0.95)   #分散共分散パラメータはセグメントで共通


#効用関数の設定
U.mean <- matrix(rowSums(XM %*% beta.z * z.vec), nrow=hhpt, ncol=member-1, byrow=T)   #効用関数の平均構造
U <- U.mean + mvrnorm(hhpt, rep(0, member-1), Cov)

#対数尤度の真の値
inverse.Cov <- solve(Cov)
det.Cov <- det(Cov)
LLt <- apply(cbind(U, U.mean), 1, function(x) dmv(x[1:(member-1)], x[member:(2*(member-1))], inverse.Cov, det.Cov))
LL_true <- sum(log(LLt))


#効用最大化原理に基づき選択メンバー決定
y <- apply(U, 1, function(x) ifelse(max(x) < 0, member, which.max(x)))

#選択メンバーを0、1行列に変更
Y <- matrix(0, nrow=hhpt, ncol=member)
for(i in 1:hhpt){
  Y[i, y[i]] <- 1
}


table(y)   #選択メンバーの集計
round(cbind(ID, y, U, U.mean), 2)   #選択メンバーと効用を比較

####マルコフ連鎖モンテカルロ法による階層混合多項プロビットモデルの推定####
####MCMC推定のための推定準備####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##多変量正規分布の尤度関数
dmv <- function(x, mean.vec, inv.cov, det.cov){
  LLo <- det.cov^(-1/2) * exp(-1/2*(x - mean.vec) %*% inv.cov %*% (x - mean.vec))
  return(LLo)
}

##多項ロジットモデルの尤度関数
LL_logit <- function(theta, Z, Zx, hh, g){
  #ロジットの計算
  logit <- matrix(Zx %*% theta, nrow=hh, ncol=g, byrow=T)
  
  #確率と対数尤度の和を計算
  P <- exp(logit)/matrix(rowSums(exp(logit)), nrow=hh, ncol=g)
  LLi <- rowSums(Z * log(P))
  LL <- sum(LLi)
  val <- list(LL=LL, P=P)
  return(val)
}

LL_opt <- function(theta, Z, Zx, hh, g){
  #ロジットの計算
  logit <- matrix(Zx %*% theta, nrow=hh, ncol=g, byrow=T)
  
  #確率と対数尤度の和を計算
  P <- exp(logit)/matrix(rowSums(exp(logit)), nrow=hh, ncol=g)
  LLi <- rowSums(Z * log(P))
  LL <- sum(LLi)
  return(LL)
}

##多変量正規分布の条件付き期待値と分散を計算する関数
cdMVN <- function(mu, Cov, dependent, U){
  
  #分散共分散行列のブロック行列を定義
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent]
  Cov21 <- Cov[-dependent, dependent]
  Cov22 <- Cov[-dependent, -dependent]
  
  #条件付き分散と条件付き平均を計算
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mu[, dependent] + t(CDinv %*% t(U[, -dependent] - mu[, -dependent]))   #条件付き平均を計算
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #条件付き分散を計算
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##MCMCアルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- array(0, dim=c(R/keep))   #対数尤度の保存用


##説明変数を多次元配列化
X.array <- array(0, dim=c(member-1, ncol(X), hhpt))
selop <- matrix(1:nrow(XM), nrow=hhpt, ncol=member-1, byrow=T)

for(i in 1:hhpt){
  X.array[, , i] <- XM[selop[i, ], ]
}

YX.array <- array(0, dim=c(member-1, ncol(X)+1, hhpt))

#推定プロセスの格納配列
UM <- matrix(0, nrow=hhpt, ncol=member-1)
util.M <- matrix(0, nrow=hhpt, ncol=member-1) 


##事前分布の設定
nu <- member   #逆ウィシャート分布の自由度
V <- solve((1/10)*diag(member-1))    #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(X))  #回帰係数の平均の事前分布
Adelta <- solve(100 * diag(rep(1, ncol(X)))) #回帰係数の事前分布の分散
zeta <- rep(0, ncol(Zx.v))   #階層モデルの回帰係数の平均の事前分布
Azeta <- solve(100 * diag(rep(1, ncol(Zx.v))))    #階層モデルの回帰係数の事前分布の分散

##サンプリング結果の保存用配列
Util <- array(0, dim=c(hhpt, member-1, R/keep))
BETA1 <- matrix(0, nrow=R/keep, ncol=ncol(XM))
BETA2 <- matrix(0, nrow=R/keep, ncol=ncol(XM))
BETA3 <- matrix(0, nrow=R/keep, ncol=ncol(XM))
SIGMA <- matrix(0, nrow=R/keep, ncol=nrow(Cov)^2)
THETA <- matrix(0, nrow=R/keep, ncol=ncol(Zx.v))
Z.VEC <- matrix(0, nrow=R/keep, ncol=hh)
Prob <- matrix(0, nrow=R/keep, ncol=hh)


##初期値の設定
##選択モデルの初期値の設定
#選択モデルの回帰係数の初期値
beta00.z <- matrix(colSums(Y[, -member])/(hhpt/10), nrow=g, ncol=member-1, byrow=T) +
  matrix(runif((member-1)*g, -1, 1), nrow=g, ncol=member-1)
beta01.z <- matrix(runif(g*2, 0, 1.4), nrow=g, ncol=2, byrow=T)
beta02.z <- matrix(runif(g*2, -1.2, 1.3), nrow=g, ncol=2, byrow=T)
beta03.z <- matrix(runif(g*(member-1), 0, 1.2), nrow=g, ncol=member-1, byrow=T)
beta04.z <- matrix(runif(g*(member-1), -1.2, 1.2), nrow=g, ncol=member-1, byrow=T)
oldbeta <- t(cbind(beta00.z, beta01.z, beta02.z, beta03.z))   #データを結合

#分散共分散行列の初期値
oldcov <- corrM(member-1, 0, 0)
invcov <- solve(oldcov)

##潜在変数モデルの初期値の設定
#thetaの初期値
oldtheta <- rep(0, ncol(Zx.v))

#潜在変数モデルの事前確率の計算
Pr <- rep(1/g, g)

#セグメントを生成
Z1 <- t(rmultinom(hh, 1, Pr))
z1 <- Z1 %*% 1:g
ZZ1 <- Z1

##効用の初期値
#zをベクトル形式に変更
Zi1 <- matrix(0, nrow=nrow(Y), ncol=g)
Z.vec <- matrix(0, nrow=nrow(XM), ncol=g)
pt.zeros <- c(0, pt)

for(i in 1:hh){
  r1 <- (sum(pt.zeros[1:i])*(member-1)+1):(sum(pt.zeros[1:i])*(member-1)+pt[i]*(member-1))
  r2 <- (sum(pt.zeros[1:i])+1):(sum(pt.zeros[1:i])+pt[i])
  
  Z.vec[r1, ] <- matrix(Z1[i, ], nrow=pt[i]*(member-1), ncol=g, byrow=T)
  Zi1[r2, ] <- matrix(Z1[i, ], nrow=pt[i], ncol=g, byrow=T)
}

#効用を計算
old.utilm <- matrix(rowSums(XM %*% oldbeta * Z.vec), nrow=hhpt, ncol=member-1, byrow=T)   #効用の平均構造
old.util <- old.utilm + mvrnorm(hhpt, rep(0, member-1), oldcov)   #誤差を加える


####マルコフ連鎖モンテカルロ法で階層有限混合多項プロビットモデルを推定####
for(rp in 1:R){
  
  ##選択結果と整合的な潜在効用を発生させる
  #条件付き期待値と条件付き分散を計算
  S <- rep(0, member-1)
  
  for(j in 1:(member-1)){
    MVR <- cdMVN(mu=old.utilm, Cov=oldcov, dependent=j, U=old.util)   #条件付き分布を計算
    UM[, j] <- MVR$CDmu   #条件付き期待値を取り出す
    S[j] <- sqrt(MVR$CDvar) #条件付き分散を取り出す
    
    #潜在効用を発生させる
    #切断領域を定義
    max.u <- apply(cbind(old.util[, -j], 0), 1, max)
    max.u <- ifelse(y==member, 0, max.u)
    
    #切断正規分布より潜在変数を発生
    old.util[, j] <- ifelse(y==j, rtnorm(mu=UM[, j], sigma=S[j], a=max.u, b=100), 
                            rtnorm(mu=UM[, j], sigma=S[j], a=-100, b=max.u))
    old.util[, j] <- ifelse(is.infinite(old.util[, j]), ifelse(y==j, max.u + runif(1), max.u - runif(1)), old.util[, j])
  }
  util.v <- as.numeric(t(old.util))
  
  ##betaの分布とパラメータの計算
  #z.vecとX.vec
  YX.bind <- cbind(util.v, XM)
  for(i in 1:hhpt){
    YX.array[, , i] <- YX.bind[id_r[i, ], ]
  }
  
  ##回帰モデルの回帰パラメータをセグメント別にギブスサンプリング
  #セグメントのインデックスを作成
  zi1 <- Zi1 %*% 1:g
  index.zi1 <- list()
  
  for(i in 1:g){
    index.zi1[[i]] <- subset(1:length(zi1), zi1==i)
  }
  
  #セグメントごとに　betaの平均と分散を計算
  invcov <- solve(oldcov)
  B <- matrix(0, nrow=ncol(XM), ncol=g)
  er <- matrix(0, nrow=nrow(XM), ncol=g)
  util.all <- matrix(0, nrow=nrow(XM), ncol=g)
  inv_XVX <- list()
  
  for(i in 1:g){
    xvx.vec <- rowSums(apply(X.array[, , index.zi1[[i]]], 3, function(x) t(x) %*% invcov %*% x))
    XVX <- matrix(xvx.vec, nrow=ncol(XM), ncol=ncol(XM), byrow=T)
    XVY <- rowSums(apply(YX.array[, , index.zi1[[i]]], 3, function(x) t(x[, -1]) %*% invcov %*% x[, 1]))
    
    #betaの分布の分散共分散行列のパラメータ
    inv_XVX[[i]] <- solve(XVX + Adelta)
    
    #betaの分布の平均パラメータ
    B[, i] <- inv_XVX[[i]] %*% (XVY + Adelta %*% Deltabar)   #betaの平均
    
    #多変量正規分布から回帰係数をサンプリング
    oldbeta[, i] <- mvrnorm(1, B[, i], inv_XVX[[i]])
  }
  
  #誤差を計算
  util.all <- XM %*% oldbeta
  er <- matrix(util.v, nrow=length(util.v), ncol=g) - util.all
  
  ##Covの分布のパラメータの計算とmcmcサンプリング(全体で共通)
  #逆ウィシャート分布のパラメータを計算
  R.error <- matrix(rowSums((er * Z.vec)), nrow=hhpt, ncol=member-1, byrow=T)
  Sn <- nu + hhpt
  IW.R <- V + matrix(rowSums(apply(R.error, 1, function(x) x %*% t(x))), nrow=member-1, ncol=member-1, byrow=T)
  
  #逆ウィシャート分布からCovをサンプリング
  oldcov <- rwishart(Sn, solve(IW.R))$IW
  
  ##多変量正規分布の尤度関数をセグメント別に計算
  det_Cov_hat <- det(oldcov)
  inv_Cov_hat <- solve(oldcov)
  LLi <- matrix(0, nrow=hhpt, ncol=g)
  
  for(i in 1:g){
    util.seg <- matrix(util.all[, i], nrow=hhpt, ncol=member-1, byrow=T)
    util.bind <- cbind(old.util, util.seg)
    
    #多変量正規分布の尤度を計算
    LLi[, i] <- apply(util.bind, 1, function(x) dmv(x[1:(member-1)], x[member:(2*(member-1))], inv_Cov_hat, det_Cov_hat))
  }
  
  #ID別に尤度の積を計算
  LLind <- as.matrix(data.frame(id=ID$id, L=LLi) %>%
                       dplyr::group_by(id) %>%
                       dplyr::summarise_each(funs(prod), L.1, L.2, L.3))
  
  logl <- sum(log(rowSums(LLind[, -1] * Z1)))   #対数尤度の和
  
  
  ##潜在セグメントを生成
  #潜在変数zの割当確率
  d1 <- exp(logit) * LLind[, -1]
  Pr_Z <- d1 / matrix(rowSums(d1), nrow=hh, ncol=g)
  
  #セグメントを多項分布より発生
  Z1 <- t(apply(Pr_Z, 1, function(x) rmultinom(1, 1, x)))
  z1 <- Z1 %*% 1:g
  
  #zをベクトル形式に変更
  Zi1 <- matrix(0, nrow=nrow(Y), ncol=g)
  Z.vec <- matrix(0, nrow=nrow(XM), ncol=g)
  pt.zeros <- c(0, pt)
  
  for(i in 1:hh){
    r1 <- (sum(pt.zeros[1:i])*(member-1)+1):(sum(pt.zeros[1:i])*(member-1)+pt[i]*(member-1))
    r2 <- (sum(pt.zeros[1:i])+1):(sum(pt.zeros[1:i])+pt[i])
    
    Z.vec[r1, ] <- matrix(Z1[i, ], nrow=pt[i]*(member-1), ncol=g, byrow=T)
    Zi1[r2, ] <- matrix(Z1[i, ], nrow=pt[i], ncol=g, byrow=T)
  }
  
  
  ##z1のセグメント割当から多項ロジットモデルのthetaをサンプリング
  #MHサンプリングで固定効果betaのサンプリング
  if(rp==1 | rp==100 | rp==500 | rp%%1000==0){
    print("最適化")
    res <- optim(oldtheta, LL_opt, gr=NULL, Z=Z1, Zx=Zx.v, hh=hh, g=g, method="BFGS", hessian=TRUE, 
                 control=list(fnscale=-1))
    oldtheta <- res$par
    rw <- -diag(solve(res$hessian))
    
  } else {
    
    newtheta <- oldtheta + rnorm(length(oldtheta), 0, rw)   #ランダムウォーク
    
    #ランダムウォークサンプリング
    #対数尤度と対数事前分布を計算
    lognew <- LL_logit(theta=newtheta, Z=Z1, Zx=Zx.v, hh, g)$LL
    logold <- LL_logit(theta=oldtheta, Z=Z1, Zx=Zx.v, hh, g)$LL
    logpnew <- lndMvn(newtheta, zeta, Azeta)
    logpold <- lndMvn(oldtheta, zeta, Azeta)
    
    #MHサンプリング
    alpha <- min(1, exp(lognew + logpnew - logold - logpold))
    if(alpha == "NAN") alpha <- -1
    
    #一様乱数を発生
    u <- runif(1)
    
    #u < alphaなら新しい固定効果betaを採択
    if(u < alpha){
      oldtheta <- newtheta
      
      #そうでないなら固定効果betaを更新しない
    } else {
      oldtheta <- oldtheta
    }
  }
  
  ##効用関数をとロジットの更新
  oldcov <- cov2cor(oldcov)
    
  old.utilm <- matrix(rowSums(util.all * Z.vec), nrow=hhpt, ncol=member-1, byrow=T)
  logit <- matrix(Zx.v %*% oldtheta, nrow=hh, ncol=g, byrow=T)
  
  ##サンプリング結果の保存とパラメータの確認
  if(rp%%keep==0){
    print(rp)
    mkeep <- rp/keep
    
    #サンプリング結果を保存
    BETA1[mkeep, ] <- oldbeta[, 1]
    BETA2[mkeep, ] <- oldbeta[, 2]
    BETA3[mkeep, ] <- oldbeta[, 3]
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    THETA[mkeep, ] <- oldtheta
    Z.VEC[mkeep, ] <- z1
    Util[, , mkeep] <- old.util 
    
    #パラメータを確認
    print(round(rbind(t(oldbeta), t(beta.z)), 2))
    print(round(cbind(cov2cor(oldcov), Cov), 2))
    print(round(rbind(oldtheta, theta.z), 2))
    print(round(c(colSums(Z1)/sum(Z1), r_rate), 2))
    print(logl)
    print(alpha)
  }
}

####推定結果と要約####
##サンプリング結果をプロット
matplot(BETA1[, 1:4], type="l", ylab="beta1の回帰係数", xlab="サンプリング数")
matplot(BETA1[, 5:9], type="l", ylab="beta1の回帰係数", xlab="サンプリング数")
matplot(BETA1[, 10:13], type="l", ylab="beta1の回帰係数", xlab="サンプリング数")
matplot(BETA1[, 14:17], type="l", ylab="beta1の回帰係数", xlab="サンプリング数")
matplot(BETA1[, 18:22], type="l", ylab="beta1の回帰係数", xlab="サンプリング数")
matplot(BETA2[, 1:4], type="l", ylab="beta2の回帰係数", xlab="サンプリング数")
matplot(BETA2[, 5:9], type="l", ylab="beta2の回帰係数", xlab="サンプリング数")
matplot(BETA2[, 10:13], type="l", ylab="beta2の回帰係数", xlab="サンプリング数")
matplot(BETA2[, 14:17], type="l", ylab="beta2の回帰係数", xlab="サンプリング数")
matplot(BETA2[, 18:22], type="l", ylab="beta2の回帰係数", xlab="サンプリング数")
matplot(BETA3[, 1:4], type="l", ylab="beta3の回帰係数", xlab="サンプリング数")
matplot(BETA3[, 5:9], type="l", ylab="beta3の回帰係数", xlab="サンプリング数")
matplot(BETA3[, 10:13], type="l", ylab="beta3の回帰係数", xlab="サンプリング数")
matplot(BETA3[, 14:17], type="l", ylab="beta3の回帰係数", xlab="サンプリング数")
matplot(BETA3[, 18:22], type="l", ylab="beta3の回帰係数", xlab="サンプリング数")

matplot(SIGMA[, 1:4], type="l", ylab="beta3の回帰係数", xlab="サンプリング数")
matplot(SIGMA[, 5:9], type="l", ylab="beta3の回帰係数", xlab="サンプリング数")

cbind(Z, Z1, ZZ1)

round(table(y[(Zi1 %*% 1:g)==1])/sum(table(y[(Zi1 %*% 1:g)==1])), 3)
round(table(y[(Zi1 %*% 1:g)==2])/sum(table(y[(Zi1 %*% 1:g)==2])), 3)
round(table(y[(Zi1 %*% 1:g)==3])/sum(table(y[(Zi1 %*% 1:g)==3])), 3)
