#####教師ありレコメンデーション#####
library(MASS)
library(nlme)
library(ncvreg)
library(glmnet)
library(lars)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)


####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho <- ifelse(abs(rho) < 0.1, 0, rho)
  rho[upper.tri(rho)] <- 0
  
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  X.Sigma <- eigen(Sigma)
  Lambda <- diag(X.Sigma$values)
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda)
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
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
hh <- 10000
item <- 500

##アイテム共起行列を発生させる
k <- 4
Pr <- matrix(0, nrow=hh, ncol=item)
Z <- matrix(0, nrow=hh, ncol=item)
X <- matrix(rnorm(hh*k), nrow=hh, ncol=k)

for(i in 1:item){
  print(i)
  logit <- runif(1, -5.5, -1.0) + X %*% rnorm(k, 0, 0.3) + rnorm(hh, 0, 0.5)
  Pr[, i] <- exp(logit) / (1+exp(logit))
  Z[, i] <- rbinom(hh, 1, Pr[, i])
}

#発生させた共起行列の集計と可視化
colMeans(Z)
hist(colMeans(Z), col="grey", main="アイテムの出現率", xlab="")

##応答変数を発生させる
#応答確率が妥当な水準になるまで繰り返す
#パラメータの設定
for(i in 1:1000){
  print(i)
  beta00 <- runif(1, -2.5, -1.5)
  beta01 <- rnorm(item, 0.2, 0.4)
  beta11 <- ifelse(beta01 > 0, beta01-0.2, beta01)
  
  #ロジットと確率の計算
  logit <- beta00 + Z %*% beta11
  P <- exp(logit)/(1+exp(logit))
  y <- rbinom(hh, 1, P)
  if(mean(y) < 0.45 & mean(y) > 0.4) break
}

####教師ありレコメンデーションモデルを推定####
##L1正則化ロジスティック回帰モデルでCV率を回帰
#クロスバリデーションで最適なパラメータを設定
k <- 5   #CV分割数
p1 <- 5
p2 <- 10
alpha_par <- seq(0.4, 1, length=p1)
lambda_par <- seq(0.001, 0.02, length=p2)
pars <- matrix(0, nrow=p1*p2, ncol=2)
cv_right <- c()

#サンプルを分割させておく
split0 <- split(sample(1:length(y), length(y)), 1:k)
index <- c()
split <- c()

for(i in 1:k){
  split <- c(split, split0[[i]])
  index <- c(index, ifelse(split0[[i]] > 0, i, 0))
}

##クロスバリデーションでハイパーパラメータを推定
for(i in 1:length(alpha_par)){
  print(i)
  alpha <- alpha_par[i]   #L1とL2の重みパラメータを選択
  for(j in 1:length(lambda_par)){
    r <- (i-1)*p1+j
    lambda <- lambda_par[j]   #正則化パラメータを選択
    pars[r, ] <- c(alpha, lambda)   #ハイパーパラメータを格納
    
    ##5分割クロスバリデーションでCVエラー率を算出
    cv_vec <- c() 
    
    for(v in 1:k){
      #サンプルを分割
      z1 <- Z[split[index!=v], ]
      z2 <- Z[split[index==v], ]
      y1 <- y[split[index!=v]]
      y2 <- y[split[index==v]]
      
      #glmnetでL1正則化ロジスティック回帰モデルを推定
      res <- glmnet(x=z1, y=y1, lambda=lambda, family="binomial", alpha=alpha)
      alpha <- res$a0
      beta <- res$beta
      
      #CVエラーを算出
      logit <- alpha + z2 %*% beta
      Pr <- exp(logit) / (1+exp(logit))
      cv_vec <- c(cv_vec, y2==ifelse(Pr > 0.5, 1, 0))
    }
    cv_right <- c(cv_right, mean(cv_vec))
    print(cv_right[r])
  }
}

