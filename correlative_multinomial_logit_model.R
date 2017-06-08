#####相関のある多項ロジットモデル#####
library(MASS)
library(mlogit)
library(MCMCpack)
library(bayesm)
library(caret)
library(reshape2)
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
hh <- 1500   #サンプル数
choise <- 5   #選択可能数
st <- 5   #基準ブランド
k <- 5   #説明変数の数

##説明変数の発生
#通常価格の発生
PRICE <- matrix(runif(hh*choise, 0.7, 1), nrow=hh, ncol=choise)   

#ディスカウント率の発生
DISC <- matrix(runif(hh*choise, 0, 0.3), nrow=hh, ncol=choise)

#特別陳列の発生
DISP <- matrix(0, nrow=hh, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hh, 1, r)
}

#特別キャンペーンの発生
CAMP <- matrix(0, nrow=hh, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hh, 1, r)
}

#カテゴリーロイヤルティ
ROYL <- matrix(runif(hh), nrow=hh, ncol=1)

##分散共分散行列の設定
corM <- corrM(col=choise, lower=-0.7, upper=0.8)   #相関行列を作成
Sigma <- covmatrix(col=choise, corM=corM, lower=1, upper=1.5)   #分散共分散行列
Cov <- Sigma$covariance

##パラメータの設定
beta1 <- -6.5   #価格のパラメータ
beta2 <- 6.3   #割引率のパラメータ
beta3 <- 2.0   #特別陳列のパラメータ
beta4 <- 1.8   #キャンペーンのパラメータ
beta5 <- c(1.1, 0.6, -0.5, 0.3)   #カテゴリーロイヤルティのパラメータ
beta0 <- c(0.5, 1.1, 1.4, 2.2)   #ブランド1～4の相対ベース販売力
betat <- c(beta0, beta1, beta2, beta3, beta4, beta5)

##効用を発生させ、選択されたブランドを決定
#多変量正規分布からロジットを発生
logit.l <- matrix(0, nrow=hh, ncol=st)
for(i in 1:(st-1)){
  logit.l[, i] <- beta0[i] + beta1*PRICE[, i]  + beta2*DISC[, i] + beta3*DISP[, i] + 
                beta4*CAMP[, i] + beta5[i]*ROYL
}
#基準変数のロジットを計算
logit.l[, st] <- beta1*PRICE[, st]  + beta2*DISC[, st] + beta3*DISP[, st] + beta4*CAMP[, st]

##多変量正規乱数を発生
logit <- logit.l + mvrnorm(n=hh, rep(0, choise), Cov)

##発生させたロジットから選択ブランドを決定
#ブランド選択確率を計算
Pr <- exp(logit)/rowSums(exp(logit))
colMeans(Pr); apply(Pr, 2, summary)

#選択ブランドを発生
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
colMeans(Y); apply(Y, 2, table)

round(cbind(Y %*% 1:choise, Pr), 3)

####マルコフ連鎖モンテカルロ法で相関構造のある混合型ロジスティック回帰モデルを推定####
##回帰モデルを推定するために説明変数をベクトル形式に変更設定
#IDの設定
ID <- rep(1:hh, rep(choise, hh))

#切片の設定
p <- c(1, rep(0, choise))
bp <- matrix(p, nrow=hh*(choise+1), ncol=choise, byrow=T)
BP <- subset(bp, rowSums(bp) > 0)
BP <- BP[, -st] 

#カテゴリロイヤルティの設定
ROYL.v <- matrix(0, nrow=hh*choise, ncol=choise)
for(i in 1:hh){
  ROYL.v[ID==i, ] <- diag(c(rep(ROYL[i, ], choise-1), 0))
}
ROYL.v <- ROYL.v[, -st]

#説明変数の設定
PRICE.v <- as.numeric(t(PRICE))
DISC.v <- as.numeric(t(DISC))
DISP.v <- as.numeric(t(DISP))
CAMP.v <- as.numeric(t(CAMP))

round(X <- data.frame(BP=BP, PRICE=PRICE.v, DISC=DISC.v, DISP=DISP.v, CAMP=CAMP.v, ROYL=ROYL.v), 2)   #データの結合
XM <- as.matrix(X)

#IDの設定
brand <- rep(1:choise, hh)
id <- rep(1:hh, rep(choise, hh))
ID <- data.frame(id, brand)

##MCMCアルゴリズムの設定
R <- 15000
sbeta <- 1.5
keep <- 2
llike <- c()   #対数尤度の保存用

##事前分布の設定
nu <- ncol(X)-2   #逆ウィシャート分布の自由度
V <- nu*diag(choise)   #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(X))  #回帰係数の平均の事前分布
Adelta <- 100 * diag(rep(1, ncol(X)))   #回帰係数の事前分布の分散

##サンプリング結果の保存用
Util <- array(0, dim=c(hh, choise, R/keep))
BETA <- matrix(0, nrow=R/keep, ncol=choise)
SIGMA <- array(0, dim=c(choise, choise, R/keep))

##初期値の設定
#回帰係数の初期値
oldbeta <- c(runif(choise-1, 0, 3), -3.0, 3.0, runif(2, 0, 2), runif(choise-1, -2, 3))  

#分散共分散行列の初期値
corM.f <- corrM(col=choise, lower=-0.6, upper=0.6)   #相関行列を作成
Sigma.f <- covmatrix(col=choise, corM=corM.f, lower=1, upper=5)   #分散共分散行列
oldcov <- Sigma.f$covariance


#効用の平均構造の初期値
b <- oldbeta[(choise+1):length(oldbeta)]

round(cbind(Y, Pr, exp(old.util)/rowSums(exp(old.util)), 
exp(matrix(u, nrow=hh, ncol=choise, byrow=T))/rowSums(exp(matrix(u, nrow=hh, ncol=choise, byrow=T)))), 2)


####マルコフ連鎖モンテカルロ法で推定####
##データ拡大法で多変量正規分布から潜在効用を発生させる
u <- XM %*% oldbeta   #ベクトル形式の効用の平均構造
U <- matrix(u, nrow=hh, ncol=choise, byrow=T)   #行列形式の効用の平均構造
util.M <- t(apply(U, 1, function(x) mvrnorm(1, x, oldcov)))   #多変量正規乱数で潜在効用を発生
U.old <- as.numeric(t(util.M))
I <- diag(hh)

for(rp in 1:50000){
##回帰モデルのギブスサンプリングでbetaとsigmaを推定
  u.vec <- as.numeric(t(util.M))   #潜在効用をベクトルに変更
  
  ##betaのギブスサンプリング
  oldcovi <- solve(oldcov)
  SIGMA.B <- kronecker(I, oldcovi)
  
  #回帰係数の平均構造
  B <- solve(t(XM) %*% SIGMA.B %*% XM) %*% t(XM) %*% SIGMA.B %*% u.vec   #回帰係数の最小二乗推定量
  XVX <- t(XM) %*% SIGMA.B %*% XM
  BETA.M <- solve(XVX + solve(Adelta)) %*% (XVX %*% B + solve(Adelta) %*% Deltabar)
  
  #回帰係数の分散共分散行列
  BETA.SIG <- solve(XVX + solve(Adelta))
  
  #多変量正規分布から回帰係数をサンプリング
  oldbeta <- mvrnorm(1, as.numeric(BETA.M), BETA.SIG)
  
  ##sigmaのギブスサンプリング
  #逆ウィシャート分布の自由度を計算
  Sn <- nu + hh

  #逆ウィシャート分布のパラメータを計算
  #二乗誤差を計算して和を取る
  Vi <- solve(V) 
  EE <- matrix(0, nrow=choise, ncol=choise)
  redi <- u.vec - XM %*% oldbeta
  
  for(i in 1:hh){
    r <- (i-1)*(choise)
    error <- redi[(r+1):(r+choise)] %*% t(redi[(r+1):(r+choise)])
    EE <- EE + error
  }
  R <- solve(EE + Vi)
  
  #逆ウィシャート乱数を発生
  oldcov <- rwishart(Sn, R)$IW
 
  u <- XM %*% oldbeta   #ベクトル形式の効用の平均構造
  U <- matrix(u, nrow=hh, ncol=choise, byrow=T)   #行列形式の効用の平均構造
  util.new <- t(apply(U, 1, function(x) mvrnorm(1, x, oldcov)))   #多変量正規乱数で潜在効用を発生
  util.old <- util.M
  
  dnew <- exp(util.new[, 1]) + exp(util.new[, 2]) + exp(util.new[, 3]) + 
          exp(util.new[, 4]) + exp(util.new[, 5])
  dold <- exp(util.old[, 1]) + exp(util.old[, 2]) + exp(util.old[, 3]) + 
          exp(util.old[, 4]) + exp(util.old[, 5])
  
  LLind.new <- Y[, 1]*util.new[, 1] + Y[, 2]*util.new[, 2] + Y[, 3]*util.new[, 3] + 
               Y[, 4]*util.new[, 4] + Y[, 5]*util.new[, 5] - log(dnew)
  LLind.old <- Y[, 1]*util.old[, 1] + Y[, 2]*util.old[, 2] + Y[, 3]*util.old[, 3] + 
               Y[, 4]*util.old[, 4] + Y[, 5]*util.old[, 5] - log(dold)                                                                        

  rand <- matrix(runif(hh), nrow=hh, ncol=choise)
  LLind.diff <- exp(LLind.new - LLind.old)
  alpha <- matrix(ifelse(LLind.diff > 1, 1, LLind.diff), nrow=hh, ncol=choise)

  util.M <- ifelse(alpha > rand, util.M <- util.new, util.M <- util.old)
  logl <- ifelse(alpha > rand, logl <- LLind.new, logl <- LLind.old)
  LL <- sum(logl)
  
  print(LL)
  print(rp)
  print(rbind(round(oldbeta, 3), round(betat, 3)))
  print(cbind(round(cov2cor(oldcov), 3), round(cov2cor(Cov), 3)))
}
