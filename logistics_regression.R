#####ロジスティック回帰モデル#####
####データの発生####
#set.seed(4543)
#説明変数とパラメータの設定
col <- 10   #パラメータ数
x <- rnorm(5000, 0, 3)
alpha <- 0.5
beta1 <- rnorm(10, 0, 0.7)
(beta <- c(alpha, beta1))
Xm <- matrix(x, nrow=length(x)/col, ncol=col, byrow=T)  
X <- as.data.frame(Xm)
X1 <- cbind(1, Xm)
n <- nrow(X)   #サンプル数

#選択確率の計算
p <- exp(X1 %*% as.vector(beta)) / (1 + exp(X1 %*% as.vector(beta)))
#p <- plogis(alpha + beta1 * x)   #これでもOK
mean(p)

#選択結果を発生させる
choice <- rbinom(n, 1, p)
mean(choice)

#データを結合
Xy <- cbind(choice, X) 
Xy <- as.data.frame(Xy)
round(head(Xy, 20), 4)

####ロジスティック回帰モデルを推定####
##対数尤度の設定
fr <- function(b, x, y){
  #パラメータの設定
  alpha <- b[1]
  beta <- b[2:11]
  
  #尤度を定義して合計する
  Xb <- alpha + as.matrix(X) %*% as.vector(beta) 
  p <- exp(Xb) / (1 + exp(Xb))
  LLS <- y*log(p) + (1-y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

##対数尤度を最大化する
b0 <- c(rep(0, 11))   #初期パラメータの設定
res <- optim(b0, fr, gr=NULL, x=X, y=choice, method="BFGS", hessian=TRUE, control=list(fnscale=-1))

#結果
(b <- res$par)   #推定されたパラメータ
beta   #真のパラメータ
(tval <- b/sqrt(-diag(solve(res$hessian))))   #t値
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(n)*length(b))   #BIC

#適合度
prob <- round(exp(as.matrix(X1) %*% as.vector(b)) / (1 + exp(as.matrix(X1) %*% as.vector(b))), 3)   #推定された確率
cbind(round(p, 3), prob , round(p-prob, 3))   #真の確率との誤差
rbind(beta, b, beta-b)   #パラメータの誤差

#関数を使うなら
res2 <- glm(choice ~ V1 + V2 + V3 + V4 +V5 + V6 + V7 + V8 + V9 + V10, data = Xy, family=binomial(link=logit))
summary(res2)
glmb <- coef(res2)
round(rbind(beta, b, glmb, beta-b, glmb-b), 3)   #パラメータの誤差
