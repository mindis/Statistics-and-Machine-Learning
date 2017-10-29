#####不完全観測データの生存時間解析#####
library(MASS)
library(matrixStats)
library(survival)
library(FAdist)
library(actuar)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(5362879)

####データの発生####
##データの設定
hh <- 10000
seg <- 2
dt <- 36   #観測期間
seg_id <- rep(1:seg, rep(hh/seg, seg))   #セグメントを設定

##説明変数の発生
type <- rbinom(hh, 1, 0.5)   #治療群の設定
level <- t(rmultinom(hh, 1, c(0.3, 0.25, 0.2, 0.15, 0.1))) %*% 1:5   #疾患レベル
Bun <- runif(hh, 0.5, 2)   #血中尿素窒素
Ca <- runif(hh, 0.1, 1.5)   #血清カルシウム
Hb <- runif(hh, 0.4, 1.5)   #血中ヘモグロビン
Prot <- rbinom(hh, 1, 0.4)   #ベンスジョーンズ蛋白の有無
sex <- rbinom(hh, 1, 0.6)   #性別
age <- t(rmultinom(hh, 1, c(0.2, 0.3, 0.3, 0.2)))   #年代
colnames(age) <- c("40代", "50代", "60代", "70代以上")
age <- age[, -1]

#データの結合
X <- data.frame(切片=1, type, level, Bun, Ca, Hb, Prot, sex, age)
XM <- as.matrix(X)

####応答変数と打ち切り変数の発生####
for(i in 1:1000){
  
  ##パラメータの設定
  #形状パラメータの設定
  alpha0 <- c(runif(1, 4, 6), runif(1, 5, 6))

  #スケールパラメータの設定
  beta00 <- c(runif(1, 2.0, 2.4), runif(1, 2.5, 2.95))
  beta01 <- t(matrix(c(0.03, 0.25, runif(2, 0.04, 0.1), runif(2, -0.25, -0.05), runif(2, 0.05, 0.25), runif(2, 0.1, 0.25), 
                       runif(2, -0.15, 0.15), runif(2, -0.2, 0.2), -0.1, -0.12, -0.16, -0.15, -0.2, -0.23), nrow=2, ncol=ncol(X)-1))
  beta0 <- rbind(beta00, beta01)   #パラメータの結合  
  
  ##対数ロジスティック分布およびワイブル分布から生存時間を発生
  #線形結合
  scale1 <- exp(XM[seg_id==1, ] %*% beta0[, 1])   
  scale2 <- exp(XM[seg_id==2, ] %*% beta0[, 2])

  #応答変数を発生
  y01 <- rllogis(length(scale1), shape=alpha0[1], scale=scale1)   #対数ロジスティック分布
  y02 <- rweibull(length(scale2), shape=alpha0[2], scale=scale2)   #ワイブル乱数
  
  if(abs(mean(y01)-mean(y02)) > 9) break
}

##可視化と集計
hist(y01, main="対数ロジスティック分布からの生存時間", xlab="生存時間", col="grey", breaks=20)
hist(y02, main="ワイブル分布からの生存時間", xlab="生存時間", col="grey", breaks=20)
hist(y01, breaks=seq(0, 100, 2.5), col="#FF00007F", xlim=c(0,max(y02)), main="混合分布からの生存時間", xlab="生存時間")
hist(y02, breaks=seq(0, 100, 2.5), col="#0000FF7F", add=T)
summary(y01)
summary(y02)

##打ち切りの設定
#右側打ち切りの設定
y0 <- c(y01, y02)
y1 <- ifelse(y0 > dt, dt, y0)
z <- ifelse(y0 > dt, 0, 1)

#不完全観測の打ち切りを設定
section <- 3
y1_lower <- floor(y1/section) * section
y1_upper <- ceiling(y1/section) * section
round(cbind(y1, y1_lower, y1_upper), 3)   #データの確認


####EMアルゴリズムで不完全観測データの生存時間解析モデルを推定####
##観測データの対数尤度を定義
obsll <- function(alpha, beta, y_upper, y_lower, X, z, r, hh){
  
  #線形結合
  lambda1 <- exp(X %*% beta[, 1])
  lambda2 <- exp(X %*% beta[, 2])
  
  #尤度と対数尤度の計算
  #対数ロジスティックモデルの尤度
  Li1 <- rep(0, hh)
  Li1[z==1] <- pllogis(y_upper[z==1], shape=alpha[1], scale=lambda1[z==1]) - 
                              pllogis(y_lower[z==1], shape=alpha[1], scale=lambda1[z==1])   #区間打ち切りの尤度
  Li1[z==0] <- 1 - pllogis(y_upper[z==0], shape=alpha0[1], scale=lambda1[z==0])   #右側打ち切りの尤度
  
  #ワイブルモデルの尤度
  Li2 <- rep(0, hh)
  Li2[z==1] <- pweibull(y_upper[z==1], alpha[2], lambda2[z==1]) - pweibull(y_lower[z==1], alpha[2], lambda2[z==1])   #区間打ち切りの尤度
  Li2[z==0] <- 1 - pweibull(y_upper[z==0], alpha[2], lambda2[z==0]); Li2 <- ifelse(Li2==0, 10^-100, Li2)   #右側打ち切りの尤度
  
  #尤度の結合
  Li <- cbind(Li1, Li2)
  Li/matrix(rowSums(Li), hh, 2)
  
  #潜在確率zの計算
  z0 <- matrix(r, nrow=hh, ncol=2, byrow=T) * Li
  z1 <- z0 / matrix(rowSums(z0), nrow=hh, ncol=2)   #潜在変数z
  
  #観測データの対数尤度
  LLho <- apply(matrix(r, nrow=hh, ncol=2, byrow=T) * Li, 1, sum)
  LLobz <- sum(log(LLho))
  rval <- list(LLobz=LLobz, z1=z1, Li=Li)
  return(rval)
}

##完全データの対数尤度
fr <- function(theta, y_upper, y_lower, z1, z, X, hh, index_alpha, index_beta1, index_beta2){

  #パラメータの設定
  alpha <- theta[index_alpha]
  beta <- cbind(theta[index_beta1], theta[index_beta2])
  
  #線形結合
  lambda1 <- exp(XM %*% beta[, 1])
  lambda2 <- exp(XM %*% beta[, 2])
  
  #対数尤度を計算
  #対数ロジスティックモデルの尤度
  Li1 <- rep(0, hh)
  Li1[z==1] <- pllogis(y_upper[z==1], shape=alpha[1], scale=lambda1[z==1]) - 
                        pllogis(y_lower[z==1], shape=alpha[1], scale=lambda1[z==1])   #区間打ち切りの尤度
  Li1[z==0] <- 1 - pllogis(y_upper[z==0], shape=alpha0[1], scale=lambda1[z==0])   #右側打ち切りの尤度
  
  #ワイブルモデルの尤度
  Li2 <- rep(0, hh)
  Li2[z==1] <- pweibull(y_upper[z==1], alpha[2], lambda2[z==1]) - pweibull(y_lower[z==1], alpha[2], lambda2[z==1])   #区間打ち切りの尤度
  Li2[z==0] <- 1 - pweibull(y_upper[z==0], alpha[2], lambda2[z==0]); Li2 <- ifelse(Li2==0, 10^-100, Li2)   #右側打ち切りの尤度
  
  #重み付き対数尤度の和
  LLi <- log(cbind(Li1, Li2))
  LL <- sum(z1 * LLi)   #潜在変数zの重み付き対数尤度の和
  return(LL)
}

##EMアルゴリズムの設定
iter <- 0
dl <- 100   #EMステップでの対数尤度の初期値の設定
tol <- 1

#インデックスの設定
index_alpha <- 1:2
index_beta1 <- 3:(2+ncol(X))
index_beta2 <- (3+ncol(X)):(ncol(X)*2+2)

##パラメータの初期値の設定
r <- rep(0.5, 2)   #混合率の初期値
z1 <- matrix(1, nrow=hh, ncol=2, byrow=T)

#準ニュートン法で初期パラメータを設定
alpha <- c(1, 1)
beta <- rep(0, ncol(X)*2)
theta <- c(alpha, beta)
theta[index_beta1[1]] <- 1
theta[index_beta2[1]] <- 1

#準ニュートン法で最尤推定
res <- try(optim(theta, fr, gr=NULL, y1_upper, y1_lower, z1, z, XM, hh, index_alpha, index_beta1, index_beta2, method="BFGS", 
                 hessian=FALSE, control=list(fnscale=-1, trace=TRUE)), silent=FALSE)

#推定されたパラメータを格納
alpha <- res$par[index_alpha]
beta <- matrix(res$par[-index_alpha], nrow=ncol(XM), 2)

##観測データの対数尤度と潜在変数zの初期値を設定
obzll <- obsll(alpha, beta, y1_upper, y1_lower, XM, z, r, hh)
z1 <- obzll$z1
LL1 <- obzll$LLobz


####EMアルゴリズムでパラメータを最尤推定####
while(abs(dl) >= tol){

  ##Nelder-Mead法で完全データを最尤推定(Mステップ)
  theta <- c(alpha, as.numeric(beta))
  res <- try(optim(theta, fr, gr=NULL, y1_upper, y1_lower, z1, z, XM, hh, index_alpha, index_beta1, index_beta2, method="Nelder-Mead", 
                   hessian=FALSE, control=list(fnscale=-1)), silent=TRUE)
  
  #パラメータを更新
  alpha <- res$par[index_alpha]
  beta <- matrix(res$par[-index_alpha], nrow=ncol(XM), ncol=2)
  r <- colSums(z1) / hh   #混合率の更新
  
  
  ##観測データの対数尤度を評価(Eステップ)
  #観測データの対数尤度と潜在変数zの更新
  obzll <- obsll(alpha, beta, y1_upper, y1_lower, XM, z, r, hh)
  z1 <- obzll$z1
  LL <- obzll$LLobz
  
  #アルゴリズムの収束判定
  iter <- iter + 1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

####推定結果の確認と適合度####
##推定されたパラメータと真のパラメータの比較
round(rbind(alpha, alpha0), 3)   #形状パラメータ
round(rbind(beta=as.numeric(beta), beta0=as.numeric(beta0)), 2)   #回帰パラメータ
round(rbind(r, r0=table(seg_id)/hh), 3)   #混合率
round(cbind(z1, seg=seg_id), 3)

##適合度
round(LL <- obzll$LLobz, 3)   #観測データの対数尤度
round(-2*(LL) + 2*(length(theta)+length(r)), 3)   #AIC

