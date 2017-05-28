#' # 自動レポート作成
# /*
library(knitr)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)
# */
#'
#' ## データの発生

# /* 
#set.seed(4543) 
# */

#'
#' ### 説明変数とパラメータの設定
#' #### データとパラメータの設定
#+ data
col <- 10   #パラメータ数
x <- rnorm(5000, 0, 3)
alpha <- 0.5
beta1 <- rnorm(10, 0, 0.7)
Xm <- matrix(x, nrow=length(x)/col, ncol=col, byrow=T)  
X <- as.data.frame(Xm)
X1 <- cbind(1, Xm)
n <- nrow(X)   #サンプル数

#' ### 選択確率の計算と応答変数の発生
#' #### 選択確率の計算
p <- exp(X1 %*% as.vector(beta)) / (1 + exp(X1 %*% as.vector(beta)))

# /* 
#p <- plogis(alpha + beta1 * x)   #これでもOK
# */

#' - 選択確率の平均
round(mean(p), 3)
#' - 選択確率の分位点
round(quantile(p), 3)

#' - 発生させた確率確率の分布
hist(p, breaks=25 , col="#0000ff40", border = "#0000ff", main="発生させた確率の分布", xlab="確率")


#' ### 選択結果を発生させる
choice <- rbinom(n, 1, p)
#' - 選択結果の単純集計
mean(choice)
table(choice)

# /*
#データを結合
Xy <- cbind(choice, X) 
Xy <- as.data.frame(Xy)
round(head(Xy, 10), 3)
# */

#' ## ロジスティック回帰モデルを推定
#' ### 対数尤度の設定
#+ analysis
fr <- function(b, x, y){
  #' - パラメータの設定
  alpha <- b[1]
  beta <- b[2:11]
  
  #' - 尤度を定義して合計する
  Xb <- alpha + as.matrix(X) %*% as.vector(beta) 
  p <- exp(Xb) / (1 + exp(Xb))
  LLS <- y*log(p) + (1-y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

#' ### 対数尤度を最大化する
#' #### 準ニュートン法で推定
b0 <- c(rep(0, ncol(X)+1))   #初期パラメータの設定
res <- optim(b0, fr, gr=NULL, x=X, y=choice, method="BFGS", hessian=TRUE, control=list(fnscale=-1))

#' #### 推定結果と適合度
#' - 推定されたパラメータ
round(b <- res$par, 3)   
#' - 真のパラメータ
round(beta, 3)   
#' - t値
round(tval <- b/sqrt(-diag(solve(res$hessian))), 3)   
#' - AIC
round(AIC <- -2*res$value + 2*length(res$par), 3)   
#' -BIC
round(BIC <- -2*res$value + log(n)*length(b), 3)   #BIC

#' #### 適合度
prob <- round(exp(as.matrix(X1) %*% as.vector(b)) / (1 + exp(as.matrix(X1) %*% as.vector(b))), 3)   #推定された確率
#' - 真の確率との誤差
cbind(round(p, 3), prob , round(p-prob, 3))[1:10, ]   

# /*
rbind(beta, b, beta-b)   #パラメータの誤差
# */

#' ## 関数を使うなら
#+ func
res2 <- glm(choice ~ V1 + V2 + V3 + V4 +V5 + V6 + V7 + V8 + V9 + V10, data = Xy, family=binomial(link=logit))
#' #### 関数glmの推定結果
summary(res2)

# /*
#ここから出力されません
glmb <- coef(res2)
error <- rbind(b, glmb, beta, b-beta, glmb-beta)   #パラメータの誤差
rownames(error) <- c("beta", "glmb", "btrue", "beta-btrue", "glmb-btrue")
#ここまで出力されません
# */

#' ## 推定結果の表
#+ results="asis", echo=FALSE
kable(round(error, 3))

# /*
spin("D:/Statistics/statistics_master/Simulation/auto_report_create.R")
# */


