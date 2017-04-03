#####ネステッドロジットモデル#####
library(mlogit)
library(nnet)
library(MASS)
library(plyr)
library(reshape2)
####データの発生####
#set.seed(58904)
##パラメータの設定
##ブランド1をベースブランドに設定
beta1 <- -4.2   #割引率のパラメータ
beta2 <- 3.1   #特別陳列のパラメータ
beta3 <- 2.3   #限定キャンペーンのパラメータ
beta4 <- 0.9   #ブランドロイヤルティのパラメータ
beta02 <- 2.0   #ブランド2のベース販売力
beta03 <- 1.0   #ブランド3のベース販売力  
beta04 <- 1.8   #ブランド4のベース販売力
beta05 <- 3.1   #ブランド5のベース販売力
rho1 <- 0.3   #クラスター1(ブランド1、2、3)の相関パラメータ
rho2 <- 0.5   #クラスター2(ブランド4、5)の相関パラメータ
lambda <- 0.6   #ブランドロイヤルティの繰越パラメータ
betaT <- c(beta1, beta2, beta3, beta4, beta02, beta03, beta04, beta05, rho1, rho2, lambda)   #真のパラメータ

hh <- 1000   #家計数 
pt <- round(runif(hh, 15, 35), 0)   #期間中の購買数は15～35回
hhpt <- sum(pt)   #総購買数

ID <- matrix(0, hhpt, 3)   #個人IDと購買回数
ID[, 1] <- c(1:hhpt)   #識別番号を入れる 
P <- matrix(0, hhpt, 5)   #購買確率
BUY <- matrix(0, hhpt, 5)   #購買ダミー
PRICE <- matrix(0, hhpt, 5)    #価格
DISP <- matrix(0, hhpt, 5)   #特別陳列
CAMP <- matrix(0, hhpt, 5)   #キャンペーン
ROY <- matrix(0, hhpt, 5)   #ブランドロイヤルティ

#ブランドロイヤルティの初期値
firstroy <- matrix(runif(hhpt*5), hhpt, 5)

##データを発生させる
#不釣り合いデータの繰り返し
for(i in 1:hh){
  for(j in 1:pt[i]){  
    r <- sum(pt[0:(i-1)])+j   
    #ID番号、購買回数を設定
    ID[r, 2] <- i
    ID[r, 3] <- j
    
    #ブランド1の販売価格、特別陳列、キャンペーンの有無の発生
    rn <- runif(3)
    #確率0.6で価格は0.85、確率0.2で価格は0.7、確率0.2で価格は0.6
    if(rn[1] < 0.6) SP <- 0.85 else
    {if(rn[1] < 0.8) SP <- 0.7 else SP <- 0.6}
    PRICE[r, 1] <- SP
    #確率0.3で特別陳列あり
    DISP[r, 1] <- (rn[2] > 0.7)
    #確率0.1でキャンペーンあり
    CAMP[r, 1] <- (rn[3] > 0.9)
    
    #ブランド2の販売価格、特別陳列、キャンペーンの有無の発生
    rn <- runif(3)
    #確率0.8で価格は1、確率0.15で価格は0.9、確率0.05で価格は0.65
    if(rn[1] < 0.6) SP <- 1 else
    {if(rn[1] < 0.9) SP <- 0.85 else SP <- 0.65}
    PRICE[r, 2] <- SP
    #確率0.4で特別陳列あり
    DISP[r, 2] <- (rn[2] > 0.6)
    #確率0.2でキャンペーンあり
    CAMP[r, 2] <- (rn[3] > 0.8)
    
    #ブランド3の販売価格、特別陳列、キャンペーンの有無の発生
    rn <- runif(3)
    #確率0.5で価格は0.9、確率0.3で価格は0.8、確率0.2で価格は0.6
    if(rn[1] < 0.7) SP <- 0.9 else
    {if(rn[1] < 0.85) SP <- 0.8 else SP <- 0.6}
    PRICE[r, 3] <- SP
    #確率0.3で特別陳列あり
    DISP[r, 3] <- (rn[2] > 0.7)
    #確率0.2でキャンペーンあり
    CAMP[r, 3] <- (rn[3] > 0.8)
    
    #ブランド4の販売価格、特別陳列、キャンペーンの有無の発生
    rn <- runif(3)
    #確率0.7で価格は1、確率0.15で価格は0.8、確率0.15で価格は0.7
    if(rn[1] < 0.7) SP <- 1 else
    {if(rn[1] < 0.85) SP <- 0.8 else SP <- 0.7}
    PRICE[r, 4] <- SP
    #確率0.15で特別陳列あり
    DISP[r, 4] <- (rn[2] > 0.85)
    #確率0.3でキャンペーンあり
    CAMP[r, 4] <- (rn[3] > 0.7)

    #ブランド5の販売価格、特別陳列、キャンペーンの有無の発生
    rn <- runif(3)
    #確率0.8で価格は1、確率0.1で価格は0.9、確率0.1で価格は0.85
    if(rn[1] < 0.8) SP <- 1 else
    {if(rn[1] < 0.9) SP <- 0.9 else SP <- 0.85}
    PRICE[r, 5] <- SP
    #確率0.2で特別陳列あり
    DISP[r, 5] <- (rn[2] > 0.8)
    #確率0.15でキャンペーンあり
    CAMP[r, 5] <- (rn[3] > 0.85)
    
    #ブランドロイヤルティ変数の作成
    if(j == 1) ROY[r, ] <- firstroy[r, ] else
    {ROY[r, ]<- lambda*ROY[r-1, ] + BUY[r-1, ]}
    
    ##選択確率の計算
    #効用の計算
    U1 <- beta1*PRICE[r, 1] + beta2*DISP[r, 1] + beta3*CAMP[r, 1] + beta4*ROY[r, 1]
    U2 <- beta1*PRICE[r, 2] + beta2*DISP[r, 2] + beta3*CAMP[r, 2] + beta4*ROY[r, 2] + beta02
    U3 <- beta1*PRICE[r, 3] + beta2*DISP[r, 3] + beta3*CAMP[r, 3] + beta4*ROY[r, 3] + beta03
    U4 <- beta1*PRICE[r, 4] + beta2*DISP[r, 4] + beta3*CAMP[r, 4] + beta4*ROY[r, 4] + beta04
    U5 <- beta1*PRICE[r, 5] + beta2*DISP[r, 5] + beta3*CAMP[r, 5] + beta4*ROY[r, 5] + beta05
    
    #ログサム変数の定義
    logsum1 <- log(exp(U1/rho1) + exp(U2/rho1) + exp(U3/rho1))    #クラスター1(ブランド1、2、3)のログサム変数
    logsum2 <- log(exp(U4/rho2) + exp(U5/rho2))   #クラスター2(ブランド4、5)のログサム変数
    
    #クラスターごとの選択確率
    CL1 <- exp(rho1*logsum1)/(exp(rho1*logsum1) + exp(rho2*logsum2))
    CL2 <- exp(rho2*logsum2)/(exp(rho1*logsum1) + exp(rho2*logsum2))
    
    #ブランドごとの選択確率
    P1 <- CL1 * exp(U1/rho1)/(exp(U1/rho1) + exp(U2/rho1) + exp(U3/rho1))
    P2 <- CL1 * exp(U2/rho1)/(exp(U1/rho1) + exp(U2/rho1) + exp(U3/rho1))
    P3 <- CL1 * exp(U3/rho1)/(exp(U1/rho1) + exp(U2/rho1) + exp(U3/rho1))
    P4 <- CL2 * exp(U4/rho2)/(exp(U4/rho2) + exp(U5/rho2))
    P5 <- CL2 * exp(U5/rho2)/(exp(U4/rho2) + exp(U5/rho2))
    
    Pr <- c(P1, P2, P3, P4, P5)   #選択確率の格納
    P[r, ] <- Pr
    ##選択確率より選択結果を発生させる
    BUY[r, ] <- t(rmultinom(1, 1, Pr))
  }
}
YX <- cbind(ID, BUY, PRICE, DISP, CAMP, ROY)   #データを結合
round(YX, 3)

##発生させたデータを要約集計
apply(BUY, 2, mean)   #購買率
apply(BUY, 2, table)   #購買回数
apply(PRICE, 2, mean)   #平均割引率
apply(DISP, 2, mean)   #特別陳列率
apply(CAMP, 2, mean)   #キャンペーン率
apply(ROY, 2, max)   #最大ブランドロイヤルティ
apply(ROY, 2, mean)   #平均ブランドロイヤルティ

####ネステッドロジットモデルのパラメータ推定####
##ブランドロイヤルティ変数を新しく定義
lambdaE <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)   #グリッドサーチでロイヤルティの繰越値を決めるために設定
ROYl <- list()
for(lam in 1:length(lambdaE)){
  ROYm <- matrix(0, nrow=hhpt, ncol=5)
  for(i in 1:hh){
    for(j in 1:pt[i]){
      r <- sum(pt[0:(i-1)])+j
      if(j==1) ROYm[r, ] <- firstroy[r, ] else
        ROYm[r, ] <- lambdaE[lam]*ROY[r-1, ] + BUY[r-1, ]
    }
  }
  ROYl[[lam]] <- ROYm
}  

##ネステッドロジットモデルの対数尤度の定義
frnested <- function(theta, BUY, PRICE, DISP, CAMP, ROY){
  b1 <- theta[1]   #価格弾力性
  b2 <- theta[2]   #特別陳列の係数
  b3 <- theta[3]   #キャンペーンの係数
  b4 <- theta[4]   #ブランドロイヤルティの係数
  b02 <- theta[5]   #ブランド2のベース販売力
  b03 <- theta[6]   #ブランド3のベース販売力
  b04 <- theta[7]   #ブランド4のベース販売力
  b05 <- theta[8]   #ブランド5のベース販売力
  rho1 <- theta[9]   #クラスター1の相関係数
  rho2 <- theta[10]   #クラスター2の相関係数
  
  #効用の定義
  U1 <- b1*PRICE[, 1] + b2*DISP[, 1] + b3*CAMP[, 1] + b4*ROY[, 1]
  U2 <- b1*PRICE[, 2] + b2*DISP[, 2] + b3*CAMP[, 2] + b4*ROY[, 2] + b02
  U3 <- b1*PRICE[, 3] + b2*DISP[, 3] + b3*CAMP[, 3] + b4*ROY[, 3] + b03 
  U4 <- b1*PRICE[, 4] + b2*DISP[, 4] + b3*CAMP[, 4] + b4*ROY[, 4] + b04
  U5 <- b1*PRICE[, 5] + b2*DISP[, 5] + b3*CAMP[, 5] + b4*ROY[, 5] + b05
  
  #ログサム変数の定義
  logsum1 <- log(exp(U1/rho1) + exp(U2/rho1) + exp(U3/rho1))    #クラスター1(ブランド1、2、3)のログサム変数
  logsum2 <- log(exp(U4/rho2) + exp(U5/rho2))   #クラスター2(ブランド4、5)のログサム変数
  
  #クラスターの選択確率
  CL1 <- exp(rho1*logsum1)/(exp(rho1*logsum1) + exp(rho2*logsum2))
  CL2 <- exp(rho2*logsum2)/(exp(rho1*logsum1) + exp(rho2*logsum2))
  
  #選択確率の計算
  P1 <- CL1 * exp(U1/rho1)/(exp(U1/rho1) + exp(U2/rho1) + exp(U3/rho1)) 
  P2 <- CL1 * exp(U2/rho1)/(exp(U1/rho1) + exp(U2/rho1) + exp(U3/rho1)) 
  P3 <- CL1 * exp(U3/rho1)/(exp(U1/rho1) + exp(U2/rho1) + exp(U3/rho1)) 
  P4 <- CL2 * exp(U4/rho2)/(exp(U4/rho2) + exp(U5/rho2)) 
  P5 <- CL2 * exp(U5/rho2)/(exp(U4/rho2) + exp(U5/rho2)) 
  
  #対数尤度の計算
  L <- BUY[, 1]*log(P1) + BUY[, 2]*log(P2) + BUY[, 3]*log(P3) + BUY[, 4]*log(P4) + BUY[, 5]*log(P5)
  LL <- sum(L)
  return(LL)
}

##ブランドロイヤルティのパラメータを動かしながら対数尤度を最大化
res <- list()
b0 <- c(-1, 1, 1, 0.5, 1.5, 0.6, 0.8, 1.5, 0.5, 0.5)   #パラメータの初期値
for(i in 1:length(lambdaE)){
res[[i]] <- optim(b0, frnested, gr=NULL, BUY, PRICE, DISP, CAMP, ROY=ROYl[[i]], 
                  method="BFGS", hessian=TRUE, control=list(fnscale=-1))
b0 <- res[[i]]$par
}

#対数尤度が最大のlambdaを選ぶ
value <- numeric()
for(lam in 1:length(lambdaE)){
  val <- res[[lam]]$value
  value <- c(value, val)
}
value   #推定された最大対数尤度
(maxres <- res[[which.max(value)]])   #対数尤度が最大の推定結果

##推定されたパラメータと統計量の推定値
op <- which.max(value)
(maxlambda <- lambdaE[which.max(value)])   #推定された繰越パラメータ
0.6   #真の繰越パラメータ
round(theta <- c(maxres$par, maxlambda), 2)   #推定されたパラメータ
betaT   #真のパラメーター

(tval <- theta[1:length(maxres$par)]/sqrt(-diag(solve(maxres$hessian))))   #t値
(AIC <- -2*maxres$value + 2*length(maxres$par))   #AIC
(BIC <- -2*maxres$value + log(nrow(BUY))*length(maxres$par))   #BIC

##推定された選択確率
#効用の計算
U1 <- theta[1]*PRICE[, 1] + theta[2]*DISP[, 1] + theta[3]*CAMP[, 1] + theta[4]*ROYl[[op]][, 1]
U2 <- theta[1]*PRICE[, 2] + theta[2]*DISP[, 2] + theta[3]*CAMP[, 2] + theta[4]*ROYl[[op]][, 2] + theta[5]
U3 <- theta[1]*PRICE[, 3] + theta[2]*DISP[, 3] + theta[3]*CAMP[, 3] + theta[4]*ROYl[[op]][, 3] + theta[6]
U4 <- theta[1]*PRICE[, 4] + theta[2]*DISP[, 4] + theta[3]*CAMP[, 4] + theta[4]*ROYl[[op]][, 4] + theta[7]
U5 <- theta[1]*PRICE[, 5] + theta[2]*DISP[, 5] + theta[3]*CAMP[, 5] + theta[4]*ROYl[[op]][, 5] + theta[8]

#ログサム変数の定義
logsum1 <- log(exp(U1/theta[9]) + exp(U2/theta[9]) + exp(U3/theta[9]))   #クラスター1(ブランド1、2、3)のログサム変数
logsum2 <- log(exp(U4/theta[10]) + exp(U5/theta[10]))   #クラスター2(ブランド4、5)のログサム変数

#クラスターごとの選択確率
CL1 <- exp(theta[9]*logsum1)/(exp(theta[9]*logsum1) + exp(theta[10]*logsum2))
CL2 <- exp(theta[10]*logsum2)/(exp(theta[9]*logsum1) + exp(theta[10]*logsum2))

#ブランドごとの選択確率
P1 <- CL1 * exp(U1/theta[9])/(exp(U1/theta[9]) + exp(U2/theta[9]) + exp(U3/theta[9]))
P2 <- CL1 * exp(U2/theta[9])/(exp(U1/theta[9]) + exp(U2/theta[9]) + exp(U3/theta[9]))
P3 <- CL1 * exp(U3/theta[9])/(exp(U1/theta[9]) + exp(U2/theta[9]) + exp(U3/theta[9]))
P4 <- CL2 * exp(U4/theta[10])/(exp(U4/theta[10]) + exp(U5/theta[10]))
P5 <- CL2 * exp(U5/theta[10])/(exp(U4/theta[10]) + exp(U5/theta[10]))

Pr <- data.frame(ID[, -1], P1, P2, P3, P4, P5, P)   #選択確率の格納
names(Pr) <- c("hh", "pt", "P1", "P2", "P3", "P4", "P5","Pr1", "Pr2", "Pr3", "Pr4", "Pr5")   #名前の変更
round(Pr, 2)   #データの確認

##要約集計
round(meanP <- apply(Pr[, 3:12], 2, mean), 2)   #平均選択確率
round(quantileP <- apply(Pr[, 3:ncol(Pr)], 2, quantile), 2)   #選択確率の四分位点
round(summaryP <- apply(Pr[, 3:ncol(Pr)], 2, summary), 2)   #選択確率の要約統計量