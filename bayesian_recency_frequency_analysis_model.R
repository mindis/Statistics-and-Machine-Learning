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
pt <- rpois(hh, 15.0)
pt <- ifelse(pt==0, 1, pt)
hhpt <- sum(pt)
dt <- 150

##IDの設定
id <- rep(1:hh, pt)
time <- c()
for(i in 1:hh){time <- c(time, 1:pt[i])}
ID <- data.frame(no=1:hhpt, id=id, time=time)

####説明変数の発生####
##階層モデルの説明変数の発生
k1 <- 5
cont1 <- 2; bin1 <- 3
X1 <- matrix(0, nrow=hh, ncol=k1)
for(i in 1:hh){
  X1[i, 1:cont1] <- rnorm(2, 0, 1)
  for(j in 1:bin1){
    X1[i, (cont1+j)] <- rbinom(1, 1, runif(1, 0.4, 0.6))
  }
}

##時変の説明変数
Price <- scale(rbeta(hhpt, 7.5, 1.5))   #旧作の平均値引率
Disc <- scale(rbeta(hhpt, 2.0, 6.5))   #旧作の値引きゲームの割合
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
ZX <- X1
colnames(ZX) <- c("cont1", "cont2", "bin1", "bin2", "bin3")

XM <- cbind(1, Price, Disc, Prom, Genre)
colnames(XM) <- c("切片", "Price", "Disc", "Prom", "genre1", "genre2", "genre3", "genre4", "genre5", "genre6")

##変量効果のデザイン行列の設定
Z <- matrix(0, nrow=hhpt, ncol=hh)
for(i in 1:hh){
  Z[ID$id==i, i] <- 1 
}

####応答変数の発生####
for(i in 1:1000){
  ##パラメータの設定
  #階層モデルのパラメータの設定
  theta01 <- c(runif(cont1, 0, 0.4), runif(bin1, -0.4, 0.8))
  theta02 <- c(runif(cont1, 0, 0.4), runif(bin1, -0.3, 0.2))
  theta0 <- cbind(theta01, theta02)
  
  ##変量効果のパラメータを多変量正規分布から発生
  #分散共分散パラメータを設定
  Cov0 <- matrix(c(0.6, -0.35, -0.35, 0.5), nrow=2, ncol=2)
  
  #変量効果のパラメータを発生
  random <- mvrnorm(hh, rep(0, 2), Cov0)
  theta <- ZX %*% theta0 + random
 
  ##固定効果のパラメータを発生
  #生存モデルのパラメータの設定
  alpha0 <- runif(1, 0.8, 1.8)   #形状パラメータ
  beta0 <- c(runif(1, 0, 5.0), runif(1, 0.1, 0.4), runif(1, -0.4, -0.1), runif(1, -0.7, -0.2), runif(g, -0.4, 0))
  
  #頻度モデルのパラメータの設定
  gamma0 <- c(runif(1, 0, 0.5), runif(1, -0.25, 0), runif(1, 0, 0.25), runif(1, 0, 0.35), runif(g, 0, 0.2)) 
  
  ##生存モデルと頻度モデルの応答変数を発生
  #平均構造を設定
  scale0 <- exp(XM %*% beta0 + Z %*% theta[, 1])
  lambda0 <- exp(XM %*% gamma0 + Z %*% theta[, 2])
  
  #ワイブル分布とポアソン分布より購買間隔と購買数を発生
  y1 <- rweibull(hhpt, alpha0, scale0)
  y2 <- rpois(hhpt, lambda0)
  
  print(round(c(min(y1), max(y1), min(y2), max(y2)), 3))
  if(min(y1) > 0.01 & max(y2) < 50) break
}

Y <- cbind(y1, y2)

#発生させた応答変数を可視化
hist(y1[y1 < 100], xlab="購買間隔", main="ゲーム店への訪問間隔", col="grey")
hist(y2[y2 < 20], xlab="購買数", main="ゲームの購買数", col="grey")

####打ち切りの設定####
##ユーザーごとにT = 150まで観測
#変数の格納用リスト
ID.list <- list()
y1.list <- list()
y2.list <- list()
X.list <- list()
Z.list <- list()
z.list <- list()

#個人ごとに打ち切り変数を設定
for(i in 1:hh){
  print(i)
  y1_ind <- Y[ID$id==i, 1]
  y2_ind <- Y[ID$id==i, 2]
  z <- rep(0, length(y1_ind))
  c_sum <- cumsum(y1_ind)
  
  #累積時間がdt以上のイベントは打ち切り
  index1 <- subset(1:length(c_sum), c_sum <= dt)
  
  if(max(c_sum) <= dt){
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
  if(class(XM[ID$id==i, ])=="numeric"){
    X.list[[i]] <- XM[ID$id==i, ]
    Z.list[[i]] <- Z[ID$id==i, ]
  } else {
    X.list[[i]] <- XM[ID$id==i, ][index2, ]
    Z.list[[i]] <- Z[ID$id==i, ][index2, ]
  }
  z.list[[i]] <- z[index2]
}

#リストをベクトルあるいは行列化
y1 <- unlist(y1.list)
y2 <- unlist(y2.list) 
ID <- do.call(rbind, ID.list)
XM <- do.call(rbind, X.list)
Z <- do.call(rbind, Z.list)
z <- 1-unlist(z.list)
hhpt <- nrow(ID)

#データの確認と可視化
round(cbind(ID, y1, y2, z, XM), 3)
hist(y1, xlab="購買間隔", main="ゲーム店への訪問間隔", col="grey")
hist(y2, xlab="購買数", main="ゲームの購買数", col="grey")

cbind(ID, y1, z)

####マルコフ連鎖モンテカルロ法で階層ベイズRFモデルを推定####
##RFモデルの対数尤度関数の設定
loglike <- function(alpha, beta, gamma, theta, y1, y2, z, X, Z){
  
  #平均構造の設定
  scale <- exp(X %*% beta + Z %*% theta[, 1])
  lambda <- exp(X %*% gamma + Z %*% theta[, 2])
  
  #対数尤度を計算
  LL1 <- z*(log(scale)+log(alpha)+(alpha-1)*log(y1)) - scale*y1^alpha   #ワイブルモデルの対数尤度
  LL2 <- y2*log(lambda)-lambda - lfactorial(y2)   #ポアソンモデルの対数尤度
  LLi <- LL1 + LL2
  LL <- sum(LL1 + LL2)
  LL_val <- list(LLi=LLi, LL=LL)
  return(LL_val)
}

##対数尤度関数の最大化用関数の設定
rf_model <- function(theta, y1, y2, z, X, index_par1, index_par2){
  
  #パラメータの設定
  alpha <- theta[1]
  beta <- theta[index_par1]
  gamma <- theta[index_par2]
  
  #平均構造の設定
  scale <- exp(X %*% beta)
  lambda <- exp(X %*% gamma)
  
  #対数尤度を計算
  LL1 <- z*(log(scale)+log(alpha)+(alpha-1)*log(y1)) - scale*y1^alpha   #ワイブルモデルの対数尤度
  LL2 <- y2*log(lambda)-lambda - lfactorial(y2)   #ポアソンモデルの対数尤度
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

#固定効果のパラメータの事前分布
betas <- rep(0, ncol(XM)*2)   #回帰係数の平均の事前分布
sigma <- diag(rep(0.01), ncol(XM)*2)   #回帰係数の事前分布の分散

#階層モデルの事前分布
Deltabar <- matrix(0, nrow=ncol(ZX), ncol=2)
Adelta <- 0.01 * diag(ncol(ZX))
nu <- 2 + ncol(ZX)
V <- nu * diag(2)

##サンプリング結果の保存用配列
ALPHA <- rep(0, R/keep)
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XM)*2)
THETA0 <- matrix(0, nrow=R/keep, ncol=ncol(ZX)*2)
THETA1 <- matrix(0, nrow=R/keep, ncol=hhpt)
THETA2 <- matrix(0, nrow=R/keep, ncol=hhpt)
SIGMA <- matrix(0, nrow=R/keep, ncol=2^2)


##棄却率と対数尤度の保存用配列
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##初期値の設定
#階層モデルのパラメータの設定
oldtheta01 <- -c(runif(cont1, 0, 0.5), runif(bin1, -0.8, 0.8))
oldtheta02 <- c(runif(cont1, 0, 0.6), runif(bin1, -0.3, 0.3))
oldtheta0 <- cbind(oldtheta01, oldtheta02)

##変量効果のパラメータを多変量正規分布から発生
#分散共分散パラメータを設定
oldcov <- matrix(c(0.5, -0.2, -0.2, 0.5), nrow=2, ncol=2)
cov_inv <- solve(oldcov)

#変量効果のパラメータを発生
oldrandom <- mvrnorm(hh, rep(0, 2), Cov0)
theta_mu <- ZX %*% oldtheta0
oldtheta <- ZX %*% oldtheta0 + oldrandom

##固定効果のパラメータを発生
##パラメータの初期値とランダムウォークの分散を設定するためにRFモデルの対数尤度を最大化
#初期値の設定
oldalpha <- runif(1, 0.8, 1.8)   #形状パラメータ
oldbeta <- -c(runif(1, 1, 5.5), runif(1, 0.4, 1.0), runif(1, -0.8, -0.3), runif(1, -0.7, -0.2), runif(g, -0.4, 0))
oldgamma <- c(runif(1, 0, 0.6), runif(1, -0.5, -0.2), runif(1, 0, 0.3), runif(1, 0, 0.35), runif(g, 0, 0.2)) 
par <- c(oldalpha, oldbeta, oldgamma)

#パラメータのインデックス
index1 <- 2:(1+ncol(XM))
index2 <- (2+ncol(XM)):length(par)

#準ニュートン法で対数尤度を最大化
res <- try(optim(par, rf_model, y1=y1, y2=y2, z=z, X=XM, index_par1=index1, index_par2=index2, method="BFGS", 
                 hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)

#推定結果
par <- res$par

#ランダムウォークの分散を設定
rw <- diag(solve(-res$hessian))
rw1 <- rw[1]
rw2 <- diag(rw[2:length(par)])

#生存モデルのパラメータの初期値を設定
oldalpha <- res$par[1]   #形状パラメータの初期値
oldbeta <- par[index1]

#頻度モデルのパラメータの設定
oldgamma <- par[index2]

#パラメータを結合
olddelta <- c(oldbeta, oldgamma)
index_par1 <- 1:length(oldbeta)
index_par2 <- (length(oldbeta)+1):length(olddelta)


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##MH法で固定効果パラメータをサンプリング
  deltad <- olddelta 
  deltan <- deltad + 0.15 * mvrnorm(1, rep(0, length(deltad)), rw2)
  
  #対数尤度と対数事前分布を計算
  lognew1 <- loglike(oldalpha, deltan[index_par1], deltan[index_par2], oldtheta, y1, y2, z, XM, Z)$LL
  logold1 <- loglike(oldalpha, deltad[index_par1], deltad[index_par2], oldtheta, y1, y2, z, XM, Z)$LL
  logpnew1 <- lndMvn(deltan, betas, sigma)
  logpold1 <- lndMvn(deltad, betas, sigma)
  
  #MHサンプリング
  alpha1 <- min(1, exp(lognew1 + logpnew1 - logold1 - logpold1))
  if(alpha1 == "NAN") alpha2 <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しい固定効果betaを採択
  if(u < alpha1){
    olddelta <- deltan
    logl <- lognew1
    
    #そうでないなら固定効果betaを更新しない
  } else {
    olddelta <- deltad
    logl <- logold1
  }  
  
  ##MH法で形状パラメータをサンプリング
  alphad <- abs(oldalpha)
  alphan <- abs(alphad + rnorm(1, 0, sqrt(rw1)))

  #対数尤度と対数事前分布を計算
  lognew2 <- loglike(alphan, olddelta[index_par1], olddelta[index_par2], oldtheta, y1, y2, z, XM, Z)$LL
  logold2 <- loglike(alphad, olddelta[index_par1], olddelta[index_par2], oldtheta, y1, y2, z, XM, Z)$LL
  logpnew2 <- -1/2 * alpha_sigma^(-1) * (log(alphan) - alpha_mu)^2
  logpold2 <- -1/2 * alpha_sigma^(-1) * (log(alphad) - alpha_mu)^2
  
  #MHサンプリング
  alpha2 <- min(1, exp(lognew2 + logpnew2 - logold2 - logpold2))
  if(alpha2 == "NAN") alpha2 <- -1

  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しい固定効果betaを採択
  if(u < alpha2){
    oldalpha <- alphan
    
    #そうでないなら固定効果betaを更新しない
  } else {
    oldalpha <- alphad
  }  
  
  
  ##MH法で変量効果をサンプリング
  oldthetad <- oldtheta
  oldthetan <- oldthetad + 0.5 * mvrnorm(hh, rep(0, 2), oldcov)
  
  #事前分布の誤差を計算
  er_new <- oldthetan - theta_mu
  er_old <- oldthetad - theta_mu
  
  #対数尤度と対数事前分布を計算
  lognew3 <- loglike(oldalpha, olddelta[index_par1], olddelta[index_par2], oldthetan, y1, y2, z, XM, Z)$LLi
  logold3 <- loglike(oldalpha, olddelta[index_par1], olddelta[index_par2], oldthetad, y1, y2, z, XM, Z)$LLi
  logpnew3 <- apply(er_new, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  logpold3 <- apply(er_old, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  
  #ID別に対数尤度の和を取る
  logl.ind <- as.matrix(data.frame(logl1=lognew3, logl2=logold3, id=ID$id) %>%
                          dplyr::group_by(id) %>%
                          dplyr::summarize(new=sum(logl1), old=sum(logl2)))[, 2:3]
  
  ##MHサンプリング
  #サンプリングを採択するかどうかを決定
  rand <- runif(hh)   #一様乱数から乱数を発生
  LLind.diff <- exp(logl.ind[, 1] + logpnew3 - logl.ind[, 2] - logpold3)   #採択率を計算
  alpha <- ifelse(LLind.diff > 1, 1, LLind.diff)
  alpha <- matrix(alpha, nrow=hh, ncol=2)
  
  #alphaに基づきbetaを採択
  oldtheta.r <- ifelse(alpha > rand, oldthetan, oldthetad)   #alphaがrandを上回っていたら採択
  adopt <- sum(oldtheta[, 1]!=oldtheta.r[, 1])/hh   #採択率
  oldtheta <- oldtheta.r   #パラメータを更新
  
  
  ##多変量回帰モデルによる階層モデルのサンプリング
  out <- rmultireg(Y=oldtheta, X=ZX, Bbar=Deltabar, A=Adelta, nu=nu, V=V)
  oldtheta0 <- out$B
  oldcov <- out$Sigma
  
  #階層モデルのパラメータを更新
  cov_inv <- solve(oldcov)
  theta_mu <- ZX %*% oldtheta0
  

  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    ALPHA[mkeep] <- oldalpha
    BETA[mkeep, ] <- olddelta
    THETA0[mkeep, ] <- as.numeric(oldtheta0)
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    colnames(oldtheta0) <- c("theta01", "theta02")
   
    print(rp)
    print(round(c(logl, res$value), 2))
    print(round(c(oldalpha, alpha0), 3))
    print(round(rbind(olddelta, delta_t=c(beta0, gamma0), ML_par=par[2:length(par)]), 3))
    print(round(rbind(theta=t(oldtheta0), t(theta0)), 3))
    print(round(cbind(oldcov, Cov0), 3))
    print(round(cbind(cov2cor(oldcov), cov2cor(Cov0)), 3))
    print(round(adopt, 3))
  }
}

matplot(BETA[, 6:10], type="l")
matplot(THETA0[, 1:5], type="l")
matplot(SIGMA, type="l")
XM
