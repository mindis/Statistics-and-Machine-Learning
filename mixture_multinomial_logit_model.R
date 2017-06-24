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
  CLOTH[[i]] <- CLOTH[[i]][, -c.num]
}
CLOTH[[member]] <- matrix(0, nrow=hhpt, ncol=c.num-1)


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
score.norm <- exp(rnorm(hhpt*2, 12.5, 0.5))
index.score <- sample(subset(1:length(score.norm), score.norm > 150000), hhpt)
score <- log(score.norm[index.score])
SCORE <- score

#どのメンバーの勧誘回だったか
prob <- 1/(member)
scout <- t(rmultinom(hhpt, 2, rep(prob, member)))

#メンバーで勧誘が重複しなくなるまで乱数を発生させ続ける
for(i in 1:10000){
  if(max(scout)==1) break
  index.scout <- subset(1:nrow(scout), apply(scout, 1, max) > 1)
  scout[index.scout, ] <- t(rmultinom(length(index.scout), 2, rep(prob, member)))
  print(i)
}
SCOUT <- scout

##パラメータの設定
#切片の設定
beta0 <- matrix(0, nrow=sg, ncol=member-1)
for(i in 1:sg){
  beta0[i, ] <- runif(member-1, 0.8, 5.0)
}
beta0 <- cbind(beta0, 0)

#衣装の回帰係数の設定
beta1 <- matrix(0, nrow=sg, ncol=c.num-1)
for(i in 1:sg){
  beta1[i, ] <- runif(c.num-1, -2.0, 3.0)
}

beta2 <- runif(sg, 0.6, 4.0)   #勧誘の回帰係数
beta3 <- c(runif(member-1, -0.3, 0.3), 0)   #レベルの回帰係数
beta4 <- c(runif(member-1, -0.15, 0.15), 0)   #スコアの回帰係数

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
##完全データのロジットモデルの尤度
cll <- function(x, Y, ones, CLOTH, SCOUT, LV, SCORE, zpt, hhpt, sg, member, c.num, l){
  b0 <- matrix(x[l[1]:l[2]], nrow=member-1, ncol=sg)
  b1 <- matrix(x[l[3]:l[4]], nrow=c.num-1, ncol=sg)
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
  b1 <- matrix(x[l[3]:l[4]], nrow=c.num-1, ncol=sg)
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

##EMアルゴリズムの設定
iter <- 0
par.cnt <- (member-1)*sg + (c.num-1)*sg + sg + member-1 + member-1   
cuml <- cumsum(c(length(beta0[, -10]), length(beta1), length(beta2), length(beta3[-member]), length(beta4[-member])))
p.len <- as.numeric(rbind(c(1, (cuml[1:4]+1)), cuml))   #パラメータベクトルの指示変数
ones <- rep(1, hhpt)   #切片の設定
dl <- 100   #EMステップでの対数尤度の差の初期値を設定
tol <- 1
maxit <- c(10, 20)   #準ニュートン法のステップ数

##EMアルゴリズムの初期値の設定
#ベストな初期パラメータを選択
r <- c(0.2, 0.2, 0.3, 0.3)
zpt <- matrix(0, nrow=hhpt, ncol=sg)

rp <- 200   #繰り返し数
Z <- list()
val <- c()
x <- matrix(0, nrow=rp, ncol=par.cnt)

for(i in 1:rp){
  #初期パラメータの設定
  x[i, ] <- c(runif((member-1)*sg, 0.5, 5), runif((c.num-1)*sg, -1, 2), runif(sg, 0.5, 2), 
              runif(2*(member-1), -0.2, 0.2))   
  
  #観測データの対数尤度の計算
  oll <- ollz(x=x[i, ], Y=Y, r=r, ones=ones, CLOTH=CLOTH, SCOUT=SCOUT, LV=LV, SCORE=SCORE, hhpt=hhpt, hh=hh,
              sg=sg, member=member, c.num=c.num, l=p.len)
  
  #パラメータの出力
  val <- c(val, oll$LLo)
  Z[[i]] <- oll$z1
  print(i)
}

#ベストな対数尤度でのパラメータ
opt <- which.max(val)
z <- Z[[opt]]
beta <- x[opt, ]
LL1 <- val[opt]

##EMアルゴリズムによる有限混合ロジットモデルの推定
for(j in 1:2){
  dl <- 10
 
  while(abs(dl) >= tol){   #dlがtol以上なら繰り返す
    for(i in 1:hh){
      zpt[ID$id==i, ] <- matrix(z[i, ], nrow=length(ID$id[ID$id==i]), ncol=sg, byrow=T)
    }
    #完全データでのロジットモデルの推定(Mステップ)
    res <- optim(beta, cll, Y=Y, ones=ones, CLOTH=CLOTH, SCOUT=SCOUT, LV=LV, SCORE=SCORE, zpt=zpt, hhpt=hhpt,
                 sg=sg, member=member, c.num=c.num, l=p.len, method="BFGS", hessian=FALSE, 
                 control=list(fnscale=-1, maxit=maxit[j]))
    beta <- res$par   #パラメータの更新
    r <- apply(z, 2, sum)/hh   #混合率の計算
    
    #Eステップでの対数尤度の期待値の計算
    obsllz <- ollz(x=beta, Y=Y, r=r, ones=ones, CLOTH=CLOTH, SCOUT=SCOUT, LV=LV, SCORE=SCORE, hhpt=hhpt, hh=hh,
                   sg=sg, member=member, c.num=c.num, l=p.len)
    LL <- obsllz$LLo
    z <- obsllz$z1
    
    #EMアルゴリズムのパラメータの更新
    iter <- iter+1
    dl <- LL-LL1
    LL1 <- LL
    print(LL)
  }
}

####推定結果と要約####
##推定されたパラメータと真のパラメータの比較
#推定されたパラメータ
#切片の回帰係数
round(t(matrix(beta[p.len[1]:p.len[2]], nrow=member-1, ncol=sg)), 3)  
round((beta0[, -member]), 3)

#衣装の回帰係数
round(t(matrix(beta[p.len[3]:p.len[4]], nrow=c.num-1, ncol=sg)), 3)  
round(beta1, 3)

#勧誘の回帰係数
round((beta[p.len[5]:p.len[6]]), 3)  
round(beta2, 3)

#LVの回帰係数
round(beta[p.len[7]:p.len[8]], 3)   
round(beta3, 3)

#スコアの回帰係数
round(beta[p.len[9]:p.len[10]], 3)   
round(beta4, 3)

##混合率とセグメントへの所属確率
round(r, 3)   #混合率
round(z, 3)   #潜在確率
apply(z, 1, which.max)   #セグメントへの所属
matplot(z[, ], ylab="セグメントへの所属確率", xlab="サンプルID", main="個人ごとのセグメント所属確率")

##AICとBICの計算
round(LL, 3)   #最大化された観測データの対数尤度
round(AIC <- -2*LL + 2*(length(res$par)+sg-1), 3)   #AIC
round(BIC <- -2*LL + log(hhpt)*length(res$par+sg-1), 3) #BIC
