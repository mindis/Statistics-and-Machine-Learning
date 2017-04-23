#####L1正則化ロジスティック回帰モデル#####
library(MASS)
library(ncvreg)
library(glmnet)
library(lars)
library(reshape2)
library(plyr)

####データの発生####
#set.seed(319)
p <- 100   #説明変数の個数
pm <- 10   #意味のある変数の個数
n <- 2000   #サンプル数
X <- matrix(rnorm(n*p), nrow=n, ncol=p)   #データの発生
b <- c(1.1, runif(pm, -1.0, 1.0), rep(0, length=(p-pm)))   #真の回帰係数

Xb <- b[1] + X %*% b[2:length(b)]   #線形結合
pr <- exp(Xb) / (1+exp(Xb))   #真の確率
y <- rbinom(n, 1, pr)   #観測データの生成


####L1正則化ロジスティック回帰モデルを推定####
##対数尤度関数の設定
#対数尤度の設定
fr <- function(b1, b2, lambda, x1, x2, y){
  #パラメータの設定
  alpha <- b1[1]
  beta <- b1[2]
  
  #尤度を定義して合計する
  Xb <- alpha + as.matrix(x1) %*% beta + as.matrix(x2) %*% b2 
  p <- exp(Xb) / (1 + exp(Xb))
  LLS <- y*log(p) + (1-y)*log(1-p) - lambda*abs(sum(b2)) - lambda*abs(beta)
  LL <- sum(LLS)
  return(LL)
}

##初期値の設定
b0 <- c(0.5, runif(p, -1, 1))
lambda <- 0.005 

##アルゴリズムの設定
max.iter <- 30
iter <- 1
tol <- 1
diff <- 100
L1 <- 0

##L1正則化ロジスティック回帰モデルを推定
while(diff >= tol & iter <= max.iter){
  for(i in 1:p){
    b1 <- b0[c(1, (i+1))]   #最適化するパラメータ
    b2 <- b0[c(-1,-(i+1))]   #固定するパラメータ
    x1 <- X[, i]   #最適化するパラメータに対応する説明変数
    x2 <- X[, -i]   #固定するパラメータに対応する説明変数
      
    #初期値をNelder-Mead法で決定する
    res1 <- optim(b1, fr, gr=NULL, b2, lambda, x1, x2, y,
                  method="Nelder-Mead", hessian=FALSE, control=list(fnscale=-1))
    
    #更新したパラメータの格納
    b0[1] <- res1$par[1]
    b0[(i+1)] <- res1$par[2]
  }
  #アルゴリズムのパラメータの更新
  iter <- iter+1
  LL <- res1$value
  diff <- abs(LL-L1)
  L1 <- LL
  
  #対数尤度の表示
  print(diff)
}

res1$value   #対数尤度
round(b0, 2)   #推定されたパラメータ
round(b, 2)   #真のパラメータ


####クロスバリテーションで最適なlambdaを決める####
##アルゴリズムとパラメータの設定
#パラメータの設定
b0 <- c(0.5, runif(p, -1, 1))   #初期値
lambdaE <- seq(0.001, 0.03, length=12)   #lambdaの候補

##アルゴリズムの設定
spl <- 5
len <- nrow(X)/spl   #サンプルを5分割
max.iter <- 30   #最大繰り返し数
iter <- 1   #繰り返し数の初期値
tol <- 1   #尤度の差のしきい値
diff <- 100   #尤度の差の初期値
L1 <- 0   #尤度の初期値
CR <- c()   #正答率の格納用

##5分割クロスバリデーションによる最適なlambdaの選択
for(lam in 1:length(lambdaE)){
  lambda <- lambdaE[lam]
  cr <- c()   
  
  for(k in 1:spl){
    l <- ((k-1)*len+1):(k*len)
    index <- subset(l, l <= nrow(X))   #デザイン行列の行数以上のインデックスは削除
    x.cv <- X[-index, ]
    y.cv <- y[-index]
    diff <- 100   #diffの初期化
    iter <- 30   #iterの初期化
    
    ##coordinate Descent法による推定
    while(diff >= tol & iter <= max.iter){
      for(i in 1:p){
        b1 <- b0[c(1, (i+1))]   #最適化するパラメータ
        b2 <- b0[c(-1,-(i+1))]   #固定するパラメータ
        x1 <- x.cv[, i]   #最適化するパラメータに対応する説明変数
        x2 <- x.cv[, -i]   #固定するパラメータに対応する説明変数
        
        #初期値をNelder-Mead法で決定する
        res1 <- optim(b1, fr, gr=NULL, b2, lambda, x1, x2, y.cv,
                      method="Nelder-Mead", hessian=FALSE, control=list(fnscale=-1))
        
        #更新したパラメータの格納
        b0[1] <- res1$par[1]
        b0[(i+1)] <- res1$par[2]
      }
      #アルゴリズムのパラメータの更新
      iter <- iter+1
      LL <- res1$value
      diff <- abs(LL-L1)
      L1 <- LL
    }
    
    ##推定されたパラメータの適合度の評価
    #パラメータの格納
    beta0 <- b0[1]   #切片の推定値
    beta1 <- b0[2:length(b0)]   #回帰パラメータの推定値
    
    #確率の計算
    logi <- beta0 + X[index, ] %*% beta1   #テストデータのロジット
    tf <- as.numeric(exp(logi) / (1 + exp(logi)) > 0.5)   #クラス判別
    
    #正答率の計算
    cr <- c(cr, sum(diag(table(y[index], tf))) / sum(table(y[index], tf)))
    print(cr)
  }
  CR <- c(CR, mean(cr))
  cat("ζ*'ヮ')ζうっうー 現在lambdaは \n", 
      round(lambdaE[lam], 5), "で結果は \n", 
      CR , "だよ \n")
}

####最適なlambdaでL1正則化ロジスティック回帰モデルを推定####
#最適な正則化パラメータを選択
plot(lambdaE, CR, type="l", lwd=2)
round(lambdaE, 3); round(CR, 3)   #lambdaと正答率の表示
(opt.lambda <- lambdaE[which.max(CR)])   #最適なlambdaを選択
b0 <- b0 + runif(length(b0), -0.2, 0.2)   #初期パラメータ


##アルゴリズムの設定
max.iter <- 30
iter <- 1
tol <- 1
diff <- 100
L1 <- 0

##L1正則化ロジスティック回帰モデルを推定
while(diff >= tol & iter <= max.iter){
  for(i in 1:p){
    b1 <- b0[c(1, (i+1))]   #最適化するパラメータ
    b2 <- b0[c(-1,-(i+1))]   #固定するパラメータ
    x1 <- X[, i]   #最適化するパラメータに対応する説明変数
    x2 <- X[, -i]   #固定するパラメータに対応する説明変数
    
    #初期値をNelder-Mead法で決定する
    res1 <- optim(b1, fr, gr=NULL, b2, opt.lambda, x1, x2, y,
                  method="Nelder-Mead", hessian=FALSE, control=list(fnscale=-1))
    
    #更新したパラメータの格納
    b0[1] <- res1$par[1]
    b0[(i+1)] <- res1$par[2]
  }
  #アルゴリズムのパラメータの更新
  iter <- iter+1
  LL <- res1$value
  diff <- abs(LL-L1)
  L1 <- LL
  
  #対数尤度の表示
  print(diff)
}

##関数を使った場合で推定
fit <- glmnet(x=X, y=y, lambda=opt.lambda, family="binomial", alpha=1)
b.func <- c(as.numeric(fit$a0), as.numeric(fit$beta))

##推定された結果と真のパラメータを比較
res1$value   #対数尤度
round(b0, 2)   #推定されたパラメータ
round(b.func, 2)   #関数で推定されたパラメータ
round(b, 2)   #真のパラメータ


