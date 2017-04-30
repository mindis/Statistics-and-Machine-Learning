####繰り返しのないコックス回帰モデル(セミパラメトリック生存モデル)####
library(MASS)
library(survival)
library(reshape2)
library(plyr)

####データの発生####
#set.seed(9843)
t <- 100   #観測期間
n <- 1000   #サンプル数
col <- 10   #変数数

##ベースラインハザードの設定
tb <- -4.7   #初期値
trend <- numeric()
s <- seq(0.5, 0.9, length=t)
for(i in 1:t){
  r <- rnorm(5, tb, 0.07)
  sort <- sort(r)
  bi <- rbinom(1, 1, s[i])
  bb <- ifelse(bi == 1, sort[4], sort[2])
  tb <- bb
  trend <- c(trend, bb)
}
plot(1:length(trend), trend, type="l", xlab="time", lwd=2)   #トレンドのプロット
round(exp(trend)/(1+exp(trend)), 3) 

##説明変数の発生と回帰成分の設定
##静的連続変数と静的二値変数の発生
#静的変数の数
s.cont <- 5
s.bin <- 2

#静的連続変数の発生
X.scont <- matrix(rnorm(n*s.cont, 0, 1), nrow=n, ncol=s.cont)

#動的二値変数の発生
X.sbin <- matrix(0, n, s.bin)
pr_b <- runif(s.bin, 0.3, 0.7)
for(i in 1:s.bin){
  bin <- rbinom(n, 1, pr_b[i])
  X.sbin[, i] <- bin
}

##動的連続変数と動的二値変数の発生
#動的変数
d.cont <- 1
d.bin <- 2

##動的連続変数の発生
type1 <- 10   #連続変数(外生変数)のタイプ数
X.d <- matrix(0, t, 10)
for(i in 1:type1){
  X.d[, i] <- matrix(rnorm(t, 0, 1), t, 1)
}

#動的連続変数の個人ごとに割り当てる
r <- runif(type1, 0.2, 0.7)
X.dcont <- matrix(0, n, t)
X.dd <- t(X.d)
for(i in 1:n){
  m <- which.max(t(rmultinom(1, 1, r)))
  X.dcont[i, ] <- X.dd[m, ]
  print(i)
}

##動的二値変数の発生
type2 <- 10   #二値変数(二値変数)のタイプ数
X.db <- matrix(0, t, 10)

for(i in 1:type2){
  X.db[, i] <- matrix(rbinom(t, 1, runif(1, 0.3, 0.7)), t, 1)
}

#動的二値変数を個人ごとに割り当てる
r <- runif(type, 0.2, 0.7)
X.dbin1 <- matrix(0, n, t)
X.dd <- t(X.db)

for(i in 1:n){
  m <- which.max(t(rmultinom(1, 1, r)))
  X.dbin1[i, ] <- X.dd[m, ]
  print(i)
}

##内的な動的二値変数の発生
X.dbin2 <- matrix(0, n, t)

for(i in 1:n){
  r <- runif(1, -1.5, 1.5)
  for(j in 1:t){
    p <- exp(r)/(1+exp(r))
    X.dbin2[i, j] <- rbinom(1, 1, p)
    r <- r + rnorm(1, 0, 0.025)
  }
}

##データをパネルデータ形式に変更
#IDと時間および識別番号
ID <- rep(1:n, rep(t, n))
time <- rep(1:t, n)
No. <- 1:(n*t)

#ベースラインハザードをパネル形式に変換
TREND <- rep(trend, n)

#静的変数をパネル形式に変換
X.s <- matrix(0, nrow=n*t, ncol=s.cont+s.bin)

for(i in 1:n){
  x <- cbind(X.scont, X.sbin)[i, ]
  X.s[ID==i, ] <- matrix(x, nrow=t, ncol=s.cont+s.bin, byrow=T)
}

#動的変数をパネル形式に変換
xdc <- matrix(t(X.dcont), nrow=n*t, 1, byrow=T)
xdb1 <- matrix(t(X.dbin1), nrow=n*t, 1, byrow=T)
xdb2 <- matrix(t(X.dbin2), nrow=n*t, 1, byrow=T)

X.d <- cbind(xdc, xdb1, xdb2)   #データの結合

#すべてのデータを結合
round(X <- data.frame(No., ID, time, S=X.s, D=X.d), 2)

##回帰係数の設定
betasc <- runif(s.cont, -0.15, 0.15) #静的連続変数の回帰係数
betasb <- runif(s.bin, -0.2, 0.2)   #静的二値変数の回帰係数
betadc <- runif(d.cont, -0.15, 0.15)   #動的連続変数の回帰係数
betadb <- runif(d.bin, -0.2, 0.2)   #動的二値変数の回帰係数
beta <- c(betasc, betasb, betadc, betadb)   #回帰係数を結合

##リンク関数を計算
s1 <- which.max(colnames(X)=="S.1")
logit <- TREND + as.matrix(X[, s1:ncol(X)]) %*% beta
Pr <- exp(logit)/(1+exp(logit))

##イベント時間を発生させる(イベントが発生したら打ち切り)
Y <- c()
for(i in 1:n){
  yr <- c()
  for(j in 1:t){
    yr <- c(yr, rbinom(1, 1, Pr[X$ID==i & time==j, ]))
    if(max(yr[1:j])>0) break
  }
  yh <- if(max(yr)>0) which.max(yr)+round(runif(1, 0, 1), 1) else {0}
  Y <- c(Y, yh)
  print(i)
}
hist(Y, breaks=20, col="grey")   #結果の頻度を見る
table(Y)

##Cox回帰用のデータセットを作成
yt <- rep(Y, rep(t, n))
YXt <- data.frame(X, yt)

#Cox回帰のイベント時間ごとに入る集合ごとにグルーピングする
YXl <- list()
for(i in 1:n){
  if(Y[i]==0) next   #イベントが発生していないなら次のサンプルへ
  index_y <- trunc(Y[i])   #イベント時間を抽出 
  index_x <- subset(1:length(Y), Y[i] <= Y | Y == 0)   #イベントが発生していないサンプルを抽出
  YXz <- data.frame(YXt[YXt$ID %in% index_x & YXt$time %in% index_y, ], T=i)   #インデックスで絞った変数をデータを抽出
  YXl[[i]] <- YXz
  print(i)
}
YX <- do.call(rbind, YXl)
index_z <- unique(YX$T)

Z <- list()
for(i in 1:length(index_z)){
  h <- subset(1:length(YX$yt[YX$T==index_z[i]]), YX$yt[YX$T==index_z[i]]>0)
  m <- min(subset(YX$yt[YX$T==index_z[i]], YX$yt[YX$T==index_z[i]]>0))
  z <- as.numeric(YX$yt[YX$T==index_z[i]]==m)
  Z[[i]] <- z
  print(i)
}

##データセットをまとめる
Zn <- unlist(Z)   #リストをベクトル化
YXz <- data.frame(YX[2:length(YX)], Zn, Z=ifelse(YX$yt>0, 1, 0))   #データを結合
rownames(YXz) <- 1:nrow(YXz)

#説明変数だけを取り出す
val <- subset(1:ncol(YXz), colnames(YXz)==c("S.1", "D.3"))
XX <- YXz[, val[1]:val[2]]

####Cox回帰モデルを推定する####
##タイデータを特定
#イベント時間ごとの出現数を数える
(tac <- table(Y[Y!=0]))   

#イベント時間のユニークを特定し、昇順に並べ替え
yu <- unique(YXz$yt[YXz$yt!=0])   
sortlist <- order(yu)   
(y.uni <- yu[sortlist])

sum(as.numeric(rownames(tac))-y.uni)   #数値が一致しているかチェック

#イベント時間のユニークごとにタイデータを記録
index_d <- list()
index_m <- list()
for(i in 1:length(y.uni)){
  ZT <- YXz[YXz$Zn==1 & YXz$yt==y.uni[i], "T"]
  ind <- subset(1:nrow(YXz), YXz$T==ZT[1])
  mi <- subset(1:nrow(YXz), YXz$Zn==1 & YXz$yt==y.uni[i])
  index_d[[i]] <- ind
  index_m[[i]] <- mi[1:tac[i]]
  print(i)
}

##対数部分尤度を定義
fr <- function(b, D, M, X, uni){
  beta <- b[1:ncol(X)]
  link.f <- exp(as.matrix(X) %*% beta)
  LLs <- c()
  
  #タイデータを含む対数部分尤度(Efron法)
  for(i in 1:length(uni)){
    LLd <- log(prod(sum(link.f[D[[i]]]) - (1:length(M[[i]])-1)/length(M[[i]])*sum(link.f[M[[i]]])))
    LLm <- colSums(X[M[[i]], ]) %*% beta
    LLs <- c(LLs, LLm-LLd)
  }
  LL <- sum(LLs)
  return(LL)
}

##部分尤度を最大化して、Cox回帰を推定
D <- index_d   #生存時間がtのリスク集合イベント
M <- index_m   #イベント時間tに発生したサンプル
b0 <- runif(10, -0.5, 0.5)   #パラメータの初期値

#部分尤度を準ニュートン法で最大化
fit <- optim(b0, fr, gr=NULL, D=D, M=M, X=XX, uni=y.uni,
             method="BFGS", hessian=T, control=list(fnscale=-1))

fit$value   #最大化された対数尤度
round(betan <- fit$par, 2)    #推定されたパラメータ
round(exp(betan), 2)   #推定されたパラメータ(ハザード比基準)
round(exp(beta), 2)   #真のパラメータ(オッズ比基準)

betan/sqrt(-diag(solve(fit$hessian)))   #t値
(AIC <- -2*fit$value + 2*length(fit$par))   #AIC


####ベースラインハザード関数と生存関数の推定####


