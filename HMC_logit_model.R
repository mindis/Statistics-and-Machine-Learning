#####ハミルトニアンモンテカルロ法によるベイジアン多項ロジットモデル#####
library(Matrix)
library(MASS)
library(bayesm)
library(R2WinBUGS)
library(LEAPFrOG)
library(matrixStats)
library(extraDistr)
library(rstan)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(238605)

####データの発生####
#データの設定
hh <- 5000   #サンプル数
select <- 10   #選択肢数
st <- 10   #基準変数
k1 <- 3   #無条件説明変数数
k2 <- 3   #条件付き説明変数数

#IDの設定
u_id <- rep(1:hh, rep(select, hh))
s_id <- rep(1:select, hh)

##説明変数の発生
#無条件説明変数
BR_vec <- matrix(diag(1, select), nrow=hh*select, ncol=select, byrow=T)
HIST_vec <- ROY_vec <- matrix(0, nrow=hh*select, ncol=select)
for(i in 1:hh){
  index <- which(u_id==i)
  ROY_vec[index, ] <- diag(rnorm(1, 0, 1.5), select)
  HIST_vec[index, ] <- diag(rbinom(1, 1, 0.5), select)
}

#条件付き説明変数
PRICE_vec <- runif(hh*select, 0, 1.5)
DISP_vec <-  rbinom(hh*select, 1, 0.4)
CAMP_vec <- rbinom(hh*select, 1, 0.3)

#データの結合
Data <- as.matrix(data.frame(br=BR_vec[, -st], roy=ROY_vec[, -st], hist=HIST_vec[, -st], price=PRICE_vec, 
                             disp=DISP_vec, camp=CAMP_vec))
sparse_data <- as(Data, "CsparseMatrix")


##ロジットモデルから応答変数を生成
#パラメータの生成
beta_br <- runif(select-1, -2.0, 2.0)
beta_roy <- runif(select-1, -1.5, 1.5)
beta_hist <- runif(select-1, -1.2, 1.2)
beta_price <- runif(1, 1.4, 2.2)
beta_disp <- runif(1, 0.6, 1.2)
beta_camp <- runif(1, 0.7, 1.3)
beta <- betat <- c(beta_br, beta_roy, beta_hist, beta_price, beta_disp, beta_camp)

#ロジットと確率を生成
logit <- matrix(sparse_data %*% beta, nrow=hh, ncol=select, byrow=T)
Pr <- exp(logit) / rowSums(exp(logit))

#多項分布から応答変数を生成
y <- rmnom(hh, 1, Pr)
colSums(y)

#####HMCでベイジアン多項ロジットモデルを推定####
##Leap Frog法を解く関数
leapfrog <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, XM, Y, par) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, XM, Y, par) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##多項ロジットモデルの対数尤度関数
loglike <- function(y, X, beta, hh, select){
  
  #ロジットと確率の計算
  logit <- matrix(X %*% beta, nrow=hh, ncol=select, byrow=T)
  Pr <- exp(logit) / rowSums(exp(logit))
  
  LLi <- rowSums(y * log(Pr))
  LL <- sum(LLi)
  val <- list(LLi=LLi, LL=LL)
  return(val)
}

##多項ロジットモデルの対数尤度の微分関数
dloglike

mean(y[y==1] * Pr[y==1])
Pr
(y - Pr) 
sparse_data
Data

colSums((as.numeric(t(y)) - as.numeric(t(Pr))) * sparse_data)

Pr <- matrix(runif(hh*select, 0, 1), nrow=hh, ncol=select)

round(cbind(as.numeric(t(y)) , as.numeric(t(Pr)), as.numeric(t(y)) - as.numeric(t(Pr)),
            (as.numeric(t(y)) - as.numeric(t(Pr))) * Data), 3)
beta
