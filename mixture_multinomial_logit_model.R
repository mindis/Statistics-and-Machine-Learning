#####混合多項ロジットモデル#####
library(MASS)
library(mlogit)
library(flexmix)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(8437)
##データの設定
sg <- 4
hh <- 2000   #サンプル数
pt <- rpois(hh, 5); pt <- ifelse(pt==0, 1, pt)   #購買機会(購買機会数が0なら1に置き換え)
hhpt <- sum(pt)
member <- 10   #選択可能メンバー数
st <- 10   #基準メンバー
k <- 5   #説明変数の数


##IDとセグメントの設定
id <- rep(1:hh, pt)
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}
seg <- rep(0, hhpt)
for(i in 1:sg){
  r <- ((i-1)*(hh/sg)+1):((i-1)*(hh/sg)+hh/sg)
  seg[id %in% r] <- i
}
ID <- data.frame(no=1:hhpt, id=id, t=t, seg=seg)   #データの結合

##説明変数の発生
#衣装の設定
c.num <- 8
CLOTH <- list()
for(i in 1:(member-1)){
  CLOTH[[i]] <- t(rmultinom(hhpt, 1, runif(c.num)))
}
CLOTH[[member]] <- matrix(0, nrow=hhpt, ncol=c.num)

#レベルの対数
lv.weib <- round(rweibull(hh*2, 1.8, 280), 0)
index.lv <- sample(subset(1:length(lv.weib), lv.weib > 80), hh)
lv <- log(lv.weib[index.lv])

#パネルに変更
LV <- c()
for(i in 1:hh){
  LV <- c(LV, rep(lv[i], pt[i]))
}

#スコアの対数
score.norm <- exp(rnorm(hh*2, 12.5, 0.5))
index.score <- sample(subset(1:length(score.norm), score.norm > 150000), hh)
score <- log(score.norm[index.score])

#パネルに変更
SCORE <- c()
for(i in 1:hh){
  SCORE <- c(SCORE, rep(score[i], pt[i]))
}

#どのメンバーの勧誘回だったか
prob <- 1/(choise-1)
scout <- t(rmultinom(hhpt, 2, rep(prob, member-1)))

#メンバーで勧誘が重複しなくなるまで乱数を発生させ続ける
for(i in 1:10000){
  if(max(scout)==1) break
  index.scout <- subset(1:nrow(scout), apply(scout, 1, max) > 1)
  scout[index.scout, ] <- t(rmultinom(length(index.scout), 2, rep(prob, member-1)))
  print(i)
}
SCOUT <- cbind(scout, 0)


##パラメータの設定
#切片の設定
beta0 <- matrix(0, nrow=sg, ncol=member-1)
for(i in 1:sg){
  beta0[i, ] <- runif(member-1, 0.8, 5.0)
}
beta0 <- cbind(beta0, 0)

#衣装の回帰係数の設定
beta1 <- matrix(0, nrow=sg, ncol=8)
for(i in 1:sg){
  beta1[i, ] <- runif(c.num, -2.0, 3.0)
}

beta2 <- runif(sg, 0.6, 4.0)   #勧誘の回帰係数
beta3 <- c(runif(member-1, -0.4, 0.4), 0)   #レベルの回帰係数
beta4 <- c(runif(member-1, -0.2, 0.2), 0)   #スコアの回帰係数

##応答変数の発生
#ロジットの計算
U <- list()
for(s in 1:sg){
  u <- c()
  index <- subset(1:nrow(ID), ID$seg==s)
  for(m in 1:member){
    betan <- c(beta1[s, ], beta2[s], beta3[m], beta4[m])
    u <- cbind(u, beta0[s, m] + cbind(CLOTH[[m]][index, ], SCOUT[index, m], LV[index], SCORE[index]) %*% betan)
  }
  U[[s]] <- u
}
U <- do.call(rbind, U)   #リストから行列に変換

#確率の計算
Pr <- exp(U)/rowSums(exp(U))
round(Pr, 3)

#応答変数の発生
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
round(cbind(Y, Pr), 3)
round(colMeans(Y), 3); colSums(Y)

####EMアルゴリズムで有限混合ロジットモデルを推定####
par.cnt <- (member-1)*sg + c.num*sg + sg + member-1 + member-1

length(b0)+length(beta1)+length(beta2)+length(beta3[-member])+length(beta4[-member])

zpt <- matrix(0, hhpt, sg)
for(i in 1:hhpt){
  s <- ID$seg[i]
  zpt[i, s] <- 1
}
x <- c(as.numeric(t(beta0[, -10])), as.numeric(t(beta1)), beta2, beta3[-member], beta4[-member])
ones <- rep(1, hhpt)   #切片

#パラメータベクトルの長さを設定
len <- cumsum(c(length(beta0[, -10]), length(beta1), length(beta2), length(beta3[-member]), length(beta4[-member])))
l <- as.numeric(rbind(c(1, (len[1:4]+1)), len))

##完全データのロジットモデルの尤度
cll <- function(x, Y, ones, CLOTH, SCOUT, LV, SCORE, zpt, hhpt, sg, member, c.num, l){
  b0 <- matrix(x[l[1]:l[2]], nrow=member-1, ncol=sg)
  b1 <- matrix(x[l[3]:l[4]], nrow=c.num, ncol=sg)
  b2 <- x[l[5]:l[6]]
  b3 <- x[l[7]:l[8]]
  b4 <- x[l[9]:l[10]]
  
  #完全データでのセグメント別の尤度を計算して和を取る
  U <- array(0, dim=c(hhpt, sg, member))
  for(i in 1:(member-1)){
    U[, , i] <- outer(ones, b0[i, ]) + CLOTH[[i]] %*% b1 + outer(SCOUT[, i], b2) + 
                                     matrix(LV*b3[i], hhpt, sg) + matrix(SCORE*b4[i], hhpt, sg)
  }
  U[, , member] <- CLOTH[[member]] %*% b1 + outer(SCOUT[, member], b2)
  
  #対数尤度を計算
  d <- apply(exp(U[, , ]), 2, rowSums)
  
  LLS <- matrix(0, nrow=hhpt, ncol=sg)
  for(s in 1:sg){
    LLS[, s] <- rowSums(Y * U[, s, ]) - log(d[, s])
  }
  LL <- sum(zpt * LLS)
  return(LL)
}


##観測データでの尤度と潜在変数zの計算
ollz <- function(x, Y, r, ones, CLOTH, SCOUT, LV, SCORE, hhpt, hh, sg, member, c.num, l){
  b0 <- matrix(x[l[1]:l[2]], nrow=member-1, ncol=sg)
  b1 <- matrix(x[l[3]:l[4]], nrow=c.num, ncol=sg)
  b2 <- x[l[5]:l[6]]
  b3 <- x[l[7]:l[8]]
  b4 <- x[l[9]:l[10]]
  
  #効用を計算
  U <- array(0, dim=c(hhpt, sg, member))
  for(i in 1:(member-1)){
    U[, , i] <- outer(ones, b0[i, ]) + CLOTH[[i]] %*% b1 + outer(SCOUT[, i], b2) + 
      matrix(LV*b3[i], hhpt, sg) + matrix(SCORE*b4[i], hhpt, sg)
  }
  U[, , member] <- CLOTH[[member]] %*% b1 + outer(SCOUT[, member], b2)
  
  #確率を計算
  LCo <- matrix(0, nrow=hhpt, ncol=sg)
  d <- apply(exp(U[, , ]), 2, rowSums)   #多項ロジットモデルの分母
  
  for(s in 1:sg){
    P <- exp(U[, s, ]) / matrix(d[, s], nrow=hhpt, ncol=member)
    LCo[, s] <- apply(P^Y, 1, prod)
  }
  
  #ID別に尤度の積を取る
  LLho <- matrix(0, nrow=hh, ncol=sg)
  for(i in 1:hh){
    if(length(ID$id[ID$id==i])==1){
      LLho[i, ] <- LCo[ID$id==i, ] 
    } else {
      LLho[i, ] <- apply(LCo[ID$id==i, ], 2, prod)
    }
  }
  
  #観測データでの対数尤度
  LLo <- sum(log(apply(matrix(r, nrow=hh, ncol=sg, byrow=T) * LLho, 1, sum)))
  
  #潜在変数zの計算
  z0 <- matrix(r, nrow=hh, ncol=sg, byrow=T) * LLho   #潜在変数zの分子
  z1 <- z0 / matrix(rowSums(z0), nrow=hh, ncol=sg)
  
  rval <- list(LLo=LLo, z1=z1)
  return(rval)
}
