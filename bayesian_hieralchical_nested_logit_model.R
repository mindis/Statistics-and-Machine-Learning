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
hh <- 1000
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
card <- t(rmultinom(hhpt, 1, c(2/member, 2/member, (member-4)/member)))
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
  Z.bin[, i] <- rbinom(hh, 1, p.bin)  
}

#多値変数
multi.h <- 4
p.multi <- runif(multi.h)
Z.multi <- t(rmultinom(hh, 1, p.multi))
freq.min <- which.min(colSums(Z.multi))
Z.multi <- Z.multi[, -freq.min]

#データの結合
Z <- data.frame(cont=Z.cont, bin=Z.bin, multi=Z.multi)
ZM <- as.matrix(Z)
k2 <- ncol(ZM)

####応答変数の発生####
##個体間パラメータの設定
#個体間分散共分散行列を設定
Cov1 <- corrM(col=k1-2, lower=-0.55, upper=0.9, eigen_lower=0.025, eigen_upper=0.35)   #相関行列を作成
Cov2 <- corrM(col=2, lower=0, upper=0, eigen_lower=0.025, eigen_upper=0.35)   #相関行列を作成


#個体間回帰パラメータの設定
theta10 <- c(runif(1, -1.1, -0.55), runif(1, -1.4, -1.1), runif(member-1, -0.2, 0.85), runif(c_num-1, -0.5, 0.7), 
             runif(1, -0.6, -0.4), runif(1, -0.5, -0.3), runif(1, 0.05, 0.15), runif(1, -0.8, -0.5), 
             runif(1, -0.55, -0.3))

theta11 <- matrix(c(runif(2*k2, -0.4, 0.5), runif((member-1)*k2, -0.4, 0.55), runif((c_num-1)*k2, -0.4, 0.5), 
                    runif(k2, -0.4, 0.4), runif(k2, -0.4, 0.4), runif(k2, -0.1, 0.15),
                    runif(k2, -0.4, 0.4), runif(k2, -0.35, 0.35)), nrow=k2, ncol=k1-2, byrow=T)

#ログサム変数のパラメータの設定
theta20 <- c(runif(1, -0.2, 0.5), runif(1, -0.4, 0.3))
theta21 <- matrix(c(runif(k2, -0.3, 0.4), runif(k2, -0.2, 0.3)), nrow=k2, ncol=2, byrow=T)


#パラメータの結合
theta1 <- rbind(theta10, theta11)
theta2 <- rbind(theta20, theta21) 


##個体内回帰パラメータの設定
beta <- cbind(1, ZM) %*% theta1 + mvrnorm(hh, rep(0, k1-2), Cov1/3)   #回帰係数のパラメータ
rho.par <- cbind(1, ZM) %*% theta2 + mvrnorm(hh, rep(0, 2), Cov2/3)   #ログサム変数のパラメータ

#回帰係数の設定
beta1 <- beta[, c(1, 3:(ncol(beta)-2))]
beta2 <- beta[, c(2:(member-1 + c_num-1 + 2), ncol(beta)-1, ncol(beta))]

#ログサム変数のパラメータを0～1に収まるように変換しておく
rho <- exp(rho.par)/(1+exp(rho.par)) 


##応答変数の発生
#ロジットとログサム変数を定義
logit1.list <- list()
logit2.list <- list()
logsum.list <- list()

#全ユーザーの個人ごとのログサム変数とロジットを計算
for(i in 1:hh){
  print(i)
  #ログサム変数の定義
  logsum.list[[i]] <- log(1 + exp(cbind(1, XM2[ID$id==i, 1:(ncol(XM2)-2)]) %*% beta2[i, 1:(ncol(XM2)-1)]))
  
  #ロジットの定義
  logit1.list[[i]] <- cbind(1, XM1[ID$id==i, ]) %*% beta1[i, ] + 
    (matrix(logsum.list[[i]], nrow=member*c_num, ncol=2)*CARD[ID$id==i, ]) %*% rho[i, ]
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
####マルコフ連鎖モンテカルロ法の推定のための準備####
##ネステッドロジットモデルの対数尤度関数
loglike <- function(x, y1, y2, XM1, XM2, CARD, type, member, c_num, hhpt){
  
  #パラメータの設定
  beta1 <- x[c(1, 3:(length(x)-4))]
  beta2 <- x[c(2:(member-1 + c_num-1 + 2), length(x)-3, length(x)-2)]
  rho <- exp(x[(length(x)-1):length(x)])/(1+exp(x[(length(x)-1):length(x)]))   #パラメータは0～1
  
  
  #ログサム変数の定義
  logsum <- log(1 + exp(cbind(1, XM2[, 1:(ncol(XM2)-2)]) %*% beta2[1:(ncol(XM2)-1)]))  #ログサム変数
  
  #ロジットの定義
  logit1 <- cbind(1, XM1) %*% beta1 + (matrix(logsum, nrow=hhpt, ncol=2)*CARD) %*% rho    #カード所有有無のロジット
  logit2 <- cbind(1, XM2) %*% beta2   #覚醒有無のロジット
  
  
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

##準ニュートン法で対数尤度を最大化する
#推定されたパラメータを初期値の参考にする　
x <- c(runif(k1-2, -0.5, 0.5), 1.0, 0.6)
res <- optim(x, loglike, y1=y1, y2=y2, XM1=XM1, XM2=XM2, CARD=CARD, type=type, member=member, c_num=c_num, hhpt=hhpt, 
             method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))


#推定結果
x1 <- res$par
H <- sqrt(-diag(solve(res$hessian)))


##MCMCアルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
ZMi <- cbind(1, ZM)

#データの設定
n <- member*c_num
index.x <- matrix(1:hhpt, nrow=hh, ncol=member*c_num, byrow=T)   #データのインデックス
index.y <- subset(1:hhpt, is.na(y2)==TRUE)
XM1_ind <- list()
XM2_ind <- list()
CARD_ind <- list()
y1_ind <- list()
y2_ind <- list()

for(i in 1:hh){
  XM1_ind[[i]] <- cbind(1, XM1[index.x[i, ], ])
  XM2_ind[[i]] <- cbind(1, XM2[index.x[i, ], ])
  CARD_ind[[i]] <- CARD[index.x[i, ], ]
  y1_ind[[i]] <- y1[index.x[i, ]]
  y2_ind[[i]] <- y2[index.x[i, ]]
}


#対数尤度の保存用配列
lognew1 <- matrix(0, nrow=hh, ncol=1)
logold1 <- matrix(0, nrow=hh, ncol=1)
logpnew1 <- matrix(0, nrow=hh, ncol=1)
logpold1 <- matrix(0, nrow=hh, ncol=1)
lognew2 <- matrix(0, nrow=hh, ncol=1)
logold2 <- matrix(0, nrow=hh, ncol=1)
logpnew2 <- matrix(0, nrow=hh, ncol=1)
logpold2 <- matrix(0, nrow=hh, ncol=1)


##事前分布の設定
Deltabar1 <- matrix(0, nrow=ncol(cbind(1, ZM)), ncol=k1-2)   #回帰パラメータの階層モデルの回帰係数の平均の事前分布
Adelta1 <- 0.01 * diag(rep(1, ncol(ZMi)))   #回帰パラメータの階層モデルの回帰係数の分散の事前分布
Deltabar2 <- matrix(0, nrow=ncol(cbind(1, ZM)), ncol=2)   #ログサム変数のパラメータの階層モデルの回帰係数の平均の事前分布
Adelta2 <- 0.01 * diag(rep(1, ncol(ZMi)))   #ログサム変数のパラメータの階層モデルの回帰係数の分散の事前分布

nu1 <- (ncol(ZM)+1)+k1-2   #逆ウィシャート分布の自由度
V1 <- nu1 * diag(k1-2)   #逆ウィシャート分布のパラメータ
nu2 <- (ncol(ZM)+1)+2   #逆ウィシャート分布の自由度
V2 <- nu2 * diag(2)   #逆ウィシャーと分布のパラメータ


##サンプリング結果の保存用配列
#パラメータの保存用配列
BETA <- array(0, dim=c(hh, k1-2, R/keep))
RHO <- array(0, dim=c(hh, 2, R/keep))
THETA1 <- matrix(0, nrow=R/keep, ncol=ncol(cbind(1, ZM))*(k1-2))
THETA2 <- matrix(0, nrow=R/keep, ncol=ncol(cbind(1, ZM))*2)
SIGMA1 <- matrix(0, nrow=R/keep, ncol=(k1-2)*(k1-2))
SIGMA2 <- matrix(0, nrow=R/keep, ncol=2*2)

#棄却率と対数尤度の保存用配列
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##初期値の設定
#回帰ペラメータの初期値
tau1 <- mvrnorm(hh, rep(0, k1-2), diag(0.3, k1-2))
oldbetas <- matrix(x1[1:(k1-2)], nrow=hh, ncol=k1-2, byrow=T) + tau1
oldDelta1 <- ginv(t(ZMi) %*% ZMi) %*% t(ZMi) %*% oldbetas
oldVbeta1 <- 1/hh * (t(oldbetas - ZMi %*% oldDelta1) %*% (oldbetas - ZMi %*% oldDelta1))
oldVbeta_inv1 <- solve(oldVbeta1)


#ログサム変数の初期値
tau2 <- mvrnorm(hh, rep(0, 2), diag(0.3, 2))
oldrho <- matrix(c(1.0, 0.7), nrow=hh, ncol=2, byrow=T) + tau2
oldDelta2 <- ginv(t(ZMi) %*% ZMi) %*% t(ZMi) %*% oldrho
oldVbeta2 <- 1/hh * (t(oldrho - ZMi %*% oldDelta2) %*% (oldrho - ZMi %*% oldDelta2))
oldVbeta_inv2 <- solve(oldVbeta2)

#パラメータのインデックスを作成
index.par1 <- c(1, 3:(ncol(XM1_ind[[1]])+1))
index.par2 <- c(2:(member-1 + c_num-1 + 2), ncol(XM1_ind[[1]])+2, ncol(XM2_ind[[1]])+3)
index.logsum <- c(2:(member-1 + c_num-1 + 2))


##ネステッドロジットモデルの対数尤度関数
LL_nest <- function(beta, rho, y1, y2, XM1, XM2, CARD, type, member, c_num, n, index.par1, index.par2, index.logsum){
  
  #パラメータの設定
  beta1 <- beta[index.par1]
  beta2 <- beta[index.par2]
  rho1 <- exp(rho)/(1+exp(rho))   #パラメータは0～1
  
  #ログサム変数の定義
  logsum <- log(1 + exp(XM2[, 1:(ncol(XM2)-2)] %*% beta[index.logsum]))  #ログサム変数
  
  #ロジットの定義
  logit1 <- XM1 %*% beta1 + (matrix(logsum, nrow=n, ncol=2)*CARD) %*% rho1    #カード所有有無のロジット
  logit2 <- XM2 %*% beta2   #覚醒有無のロジット
  
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

LLind <- matrix(0, nrow=hh, ncol=1)
for(i in 1:hh){
  LLind[i, ] <- LL_nest(beta=beta[i, ], rho=rho.par[i, ], y1=y1_ind[[i]], y2=y2_ind[[i]], XM1=XM1_ind[[i]], XM2=XM2_ind[[i]], 
                        CARD=CARD_ind[[i]], type=type, member=member, c_num=c_num, n=n, index.par1=index.par1, 
                        index.par2=index.par2, index.logsum=index.logsum)
}
LL_true <- sum(LLind)


####マルコフ連鎖モンテカルロ法で階層ベイズネステッドロジットモデルのパラメータをサンプリング####
for(rp in 1:R){
  
  ##MH法で個人別にbetaとrhoをサンプリング
  #ランダムウォークサンプリング
  rw1 <- mvrnorm(hh, rep(0, k1-2), diag(H[1:(k1-2)]))
  rw2 <- mvrnorm(hh, rep(0, 2), diag(0.1, 2))
  betad <- oldbetas.r <- oldbetas
  betan <- oldbetas + rw1
  rhod <- oldrho.r <- oldrho
  rhon <- oldrho + rw2
  
  #パラメータの事前分布と誤差を計算
  er_new1 <- betan - ZMi %*% oldDelta1
  er_old1 <- betad - ZMi %*% oldDelta1
  er_new2 <- rhon - ZMi %*% oldDelta2
  er_old2 <- rhod - ZMi %*% oldDelta2
  
  
  ##ID別に対数尤度と対数事前分布を計算
  for(i in 1:hh){
    #回帰パラメータの対数尤度と対数事前分布
    lognew1[i, ] <- LL_nest(beta=betan[i, ], rho=oldrho[i, ], y1=y1_ind[[i]], y2=y2_ind[[i]], XM1=XM1_ind[[i]], XM2=XM2_ind[[i]], 
                            CARD=CARD_ind[[i]], type=type, member=member, c_num=c_num, n=n, index.par1=index.par1, 
                            index.par2=index.par2, index.logsum=index.logsum)
    logold1[i, ] <- LL_nest(beta=betad[i, ], rho=oldrho[i, ], y1=y1_ind[[i]], y2=y2_ind[[i]], XM1=XM1_ind[[i]], XM2=XM2_ind[[i]], 
                            CARD=CARD_ind[[i]], type=type, member=member, c_num=c_num, n=n, index.par1=index.par1, 
                            index.par2=index.par2, index.logsum=index.logsum)
    logpnew1[i, ] <- -0.5 * (er_new1[i, ] %*% oldVbeta_inv1 %*% er_new1[i, ])
    logpold1[i, ] <- -0.5 * (er_old1[i, ] %*% oldVbeta_inv1 %*% er_old1[i, ])
    
    #ログサム変数のパラメータの対数尤度と対数事前分布
    lognew2[i, ] <- LL_nest(beta=oldbetas[i, ], rho=rhon[i, ], y1=y1_ind[[i]], y2=y2_ind[[i]], XM1=XM1_ind[[i]], XM2=XM2_ind[[i]], 
                            CARD=CARD_ind[[i]], type=type, member=member, c_num=c_num, n=n, index.par1=index.par1, 
                            index.par2=index.par2, index.logsum=index.logsum)
    logold2[i, ] <- LL_nest(beta=oldbetas[i, ], rho=rhod[i, ], y1=y1_ind[[i]], y2=y2_ind[[i]], XM1=XM1_ind[[i]], XM2=XM2_ind[[i]], 
                            CARD=CARD_ind[[i]], type=type, member=member, c_num=c_num, n=n, index.par1=index.par1, 
                            index.par2=index.par2, index.logsum=index.logsum)
    logpnew2[i, ] <- -0.5 * (er_new2[i, ] %*% oldVbeta_inv2 %*% er_new2[i, ])
    logpold2[i, ] <- -0.5 * (er_old2[i, ] %*% oldVbeta_inv2 %*% er_old2[i, ])
  }
  
  ##サンプリングを採択するかどうかを決定
  #一様乱数から乱数を発生
  rand1 <- runif(hh)  
  rand2 <- runif(hh)
  
  #採択率を計算
  LLind.diff1 <- exp(lognew1 + logpnew1 - logold1 - logpold1)
  alpha1 <- ifelse(LLind.diff1 > 1, 1, LLind.diff1)
  LLind.diff2 <- exp(lognew2 + logpnew2 - logold2 - logpold2) 
  alpha2 <- ifelse(LLind.diff2 > 1, 1, LLind.diff2)
  
  #alphaに基づきbetaとrhoを採択するかどうか決定
  index.adopt1 <- subset(1:hh, alpha1 > rand1)
  oldbetas.r[index.adopt1, ] <- betan[index.adopt1, ]
  index.adopt2 <- subset(1:hh, alpha2 > rand2)
  oldrho.r[index.adopt2, ] <- rhon[index.adopt2, ]
  
  #採択率と対数尤度を計算
  adopt1 <- sum(oldbetas[, 1] != oldbetas.r[, 1])/hh
  adopt2 <- sum(oldrho[, 1] != oldrho.r[, 1])/hh
  LLho <- sum(lognew1[index.adopt1, ]) + sum(logold1[-index.adopt1, ])   #対数尤度
  
  #パラメータを更新
  oldbetas <- oldbetas.r
  oldrho <- oldrho.r
  
  
  ##多変量回帰モデルによる階層モデルのサンプリング
  #回帰パラメータの多変量回帰
  out1 <- rmultireg(Y=oldbetas, X=ZMi, Bbar=Deltabar1, A=Adelta1, nu=nu1, V=V1)
  oldDelta1 <- out1$B
  oldVbeta1 <- out1$Sigma
  oldVbeta_ind1 <- solve(oldVbeta1)
  
  #ログサムパラメータの多変量回帰
  out2 <- rmultireg(Y=oldrho, X=ZMi, Bbar=Deltabar2, A=Adelta2, nu=nu2, V=V2)
  oldDelta2 <- out2$B
  oldVbeta2 <- out2$Sigma
  oldVbeta_ind2 <- solve(oldVbeta2)
  
  
  ##サンプリング結果を保存
  if(rp%%keep==0){
    cat("ζ*'ヮ')ζ <うっうー ちょっとまってね", paste(rp/R*100, "％"), fill=TRUE)
    mkeep <- rp/keep
    BETA[, , mkeep] <- oldbetas
    RHO[, , mkeep] <- oldrho
    THETA1[mkeep, ] <- as.numeric(oldDelta1)
    THETA2[mkeep, ] <- as.numeric(oldDelta2)
    SIGMA1[mkeep, ] <- as.numeric(oldVbeta1)
    SIGMA2[mkeep, ] <- as.numeric(oldVbeta2)
    print(round(c(adopt1, adopt2), 3))
    print(round(c(LLho, LL_true, res$value), 2))
    print(round(rbind(oldbetas[1:4, ], beta[1:4, ]), 2))
    print(round(rbind(colMeans(oldbetas), colMeans(beta)), 2))
    print(round(t(cbind(exp(oldrho[1:15, ])/(1+exp(oldrho[1:15, ])), rho[1:15, ])), 3))
    print(round(c(colMeans(exp(oldrho)/(1+exp(oldrho))), colMeans(rho)), 3))
    #print(round(t(cbind(oldDelta2, theta2)), 3))
    print(round(t(rbind(oldDelta1[, 1:5], theta1[, 1:5])), 2))
    print(round(t(rbind(oldDelta2, theta2)), 2))
    #print(round(cbind(cov2cor(oldVbeta1[1:5, 1:5]), Cov1[1:5, 1:5]), 3))
    print(rp)
  }
}
