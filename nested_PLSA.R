#####入れ子型トピックモデル#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
detach("package:bayesm", unload=TRUE)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####データの発生####
#set.seed(423943)
#データを仮設定
hh0 <- 500   #ユーザー数
item0 <- 50   #アイテム数

##IDとレビュー履歴を発生
#IDを仮設定
u.id0 <- rep(1:hh0, rep(item0, hh0))
i.id0 <- rep(1:item0, hh0)

#レビュー履歴を発生
hist <- rep(0, hh0*item0)
for(i in 1:item0){
  p <- runif(1, 0.25, 0.5)
  hist[i.id0==i] <- rbinom(hh0, 1, p)
}

#レビュー履歴からIDを再設定
index <- subset(1:length(hist), hist==1)
u.id <- u.id0[index]
i.id <- i.id0[index]
ID <- data.frame(no=1:length(u.id), id=u.id, item=i.id)

#データの再設定
k <- 8   #トピック数
hh <- length(unique(u.id))   #ユーザー数
item <- length(unique(i.id))   #アイテム数
d <- length(u.id)   #文書数
v <- 300   #語彙数

#1文書あたりの単語数
freq <- as.numeric(table(u.id))
w <- c()
for(i in 1:hh){
  par <- rgamma(1, 150, 1.10)
  w <- c(w, rpois(freq[i], par))   #購買数
}

####bag of word形式の文書行列を発生####
#パラメータの設定
alpha0 <- round(runif(k, 0.2, 0.3), 3)   #ユーザーのディレクリ事前分布のパラメータ
alpha1 <- round(runif(k, 0.15, 0.3), 3)   #アイテムのディレクリ事前分布のパラメータ
alpha2 <- rep(0.125, v)   #単語のディレクリ事前分布のパラメータ
alpha0
#ディレクリ乱数の発生
beta <- rdirichlet(hh, alpha0)   #ユーザーのトピック分布をディレクリ乱数から発生
gamma <- rdirichlet(item, alpha1)   #アイテムのトピック分布をディレクリ乱数から発生 
phi <- rdirichlet(k, alpha2)   #単語のトピック分布をディレクリ乱数から発生
pi <- rbeta(sum(w), 0.3, 0.35)

#スイッチングラベル選択のためのインデックスを作成
index_list <- list()
index_pi <- list()
for(i in 1:d){index_list[[i]] <- rep(i, w[i])}
index <- unlist(index_list)
for(i in 1:d){index_pi[[i]] <- subset(1:length(index), index==i)}

#多項分布の乱数からデータを発生
WX <- matrix(0, nrow=d, ncol=v)
Pi <- matrix(0, nrow=d, ncol=v)
theta1 <- matrix(0, nrow=d, ncol=k)
theta2 <- matrix(0, nrow=d, ncol=k)
Z <- list()
y <- list()

for(i in 1:hh){
  
  r1 <- i.id[u.id==i]
  u1 <- u.id[u.id==i]
  index1 <- subset(1:length(u.id), u.id==i)
  
  for(j in 1:sum(u.id==i)){
    
    #インデックスを設定
    r2 <- r1[j]   #アイテムインデックス
    u2 <- u1[j]   #ユーザーインデックス
    index2 <- index1[j]   #文書インデックス
    
    #多項分布より文書トピックを生成
    #文書のトピック分布を発生
    z1 <- t(rmultinom(w[index2], 1, beta[u2, ]))   #ユーザーのトピックを発生   
    z2 <- t(rmultinom(w[index2], 1, gamma[r2, ]))   #アイテムのトピックの発生
    
    #どちらのトピックから発生したかを決定
    y[[index2]] <- rbinom(w[index2], 1, pi[index_pi[[index2]]])
    Y <- matrix(y[[index2]], nrow=w[index2], ncol=k)
    z <- Y*z1 + (1-Y)*z2

    #トピックをベクトル形式に変換
    zn <- z %*% c(1:k)   #0,1を数値に置き換える
    zdn <- cbind(zn, z)   #apply関数で使えるように行列にしておく
    
    #文書トピックより単語を生成
    wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi[x[1], ])))   #文書のトピックから単語を生成
    wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
    WX[index2, ] <- wdn  
    Z[[index2]] <- zdn[, 1]
    
    #トピックの発生確率を格納しておく
    freq_w0 <- colSums(wn)
    freq_w <- ifelse(freq_w0==0, 1, freq_w0)
    Pi[index2, ] <- colSums(matrix(pi[index_pi[[index2]]], nrow=length(pi[index_pi[[index2]]]), ncol=v) * wn) / freq_w
  }
  print(i)
}

storage.mode(WX) <- "integer"
gc(); gc()

####EMアルゴリズムでトピックモデルを推定####
####トピックモデルのためのデータと関数の準備####
##それぞれの文書中の単語の出現をベクトルに並べる
##データ推定用IDを作成
ID_list <- list()
pi_list <- list()
user_list <- list()
item_list <- list()
wd_list <- list()

#文書ごとにユーザーIDおよび単語IDを作成
for(i in 1:nrow(WX)){
  print(i)
  
  #IDを格納
  ID_list[[i]] <- rep(i, w[i])
  user_list[[i]] <- rep(ID$id[i], w[i])
  item_list[[i]] <- rep(ID$item[i], w[i])
  
  #出現単語をベクトル化して格納
  num11 <- (WX[i, ] > 0) * c(1:v) 
  num12 <- subset(num11, num11 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number1 <- rep(num12, W1)
  wd_list[[i]] <- number1
  
  #トピック発生確率をベクトル化して格納
  num21 <- Pi[i, ] * rep(1, v) 
  num22 <- subset(num21, num21 > 0)
  P1 <- WX[i, (WX[i, ] > 0)]
  number2 <- rep(num22, P1)
  pi_list[[i]] <- number2
}

#リストをベクトルに変換
ID_d <- unlist(ID_list)
ID_u <- unlist(user_list)
ID_i <- unlist(item_list)
wd <- unlist(wd_list)
pd <- unlist(pi_list)

##インデックスを作成
doc_list <- list()
user_list <- list()
item_list <- list()
word_list <- list()

for(i in 1:length(unique(ID_d))) {doc_list[[i]] <- subset(1:length(ID_d), ID_d==i)}
for(i in 1:length(unique(ID_u))) {user_list[[i]] <- subset(1:length(ID_u), ID_u==i)}
for(i in 1:length(unique(ID_i))) {item_list[[i]] <- subset(1:length(ID_i), ID_i==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- subset(1:length(wd), wd==i)}
gc(); gc()

##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #負担係数の格納用
  for(kk in 1:k){
    #負担係数を計算
    Bi <- rep(theta[, kk], w) * phi[kk, c(wd)]   #尤度
    Bur[, kk] <- Bi   
  }
  Br <- Bur / rowSums(Bur)   #負担率の計算
  r <- colSums(Br) / sum(Br)   #混合率の計算
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}


Bur1 <- matrix(0, nrow=length(ID_i), ncol=k)
Bur2 <- matrix(0, nrow=length(ID_i), ncol=k)

for(kk in 1:k){
  #負担係数を計算
  Bi1 <- beta[ID_u, kk] * phi[kk, c(wd)]   #尤度
  Bi2 <- gamma[ID_i, kk] * phi[kk, c(wd)]   #尤度
  Bur1[, kk] <- Bi1   
  Bur2[, kk] <- Bi2
}

a1 <- rowSums(Bur1)
a2 <- rowSums(Bur2)

round(cbind(a1/(a1+a2), pd), 3)


####EMアルゴリズムの初期値を設定する####
##初期値をランダムに設定
#phiの初期値
freq_v <- matrix(colSums(WX), nrow=k, ncol=v, byrow=T)   #単語の出現率
rand_v <- matrix(trunc(rnorm(k*v, 0, (colSums(WX)/2))), nrow=k, ncol=v, byrow=T)   #ランダム化
phi_r <- abs(freq_v + rand_v) / rowSums(abs(freq_v + rand_v))   #トピックごとの出現率をランダムに初期化

#betaの初期値
beta_r <- rdirichlet(hh, runif(k, 0.2, 3))   #ディレクリ分布から初期値を設定

#gammaの初期値
gamma_r <- rdirichlet(item, runif(k, 0.2, 3))   #ディクレリ分布から初期値を設定

#ベイズの定理よりthetaを計算
par <- beta_r[ID$id, ] * gamma_r[ID$item, ]
theta_r <- par / matrix(rowSums(par), nrow=d, ncol=k)

###パラメータの更新
##負担率の計算
bfr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, w=w, k=k)
Br <- bfr$Br   #負担率
r <- bfr$r   #混合率

#betaの更新
usum <- (data.frame(id=ID_u, Br=Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
beta_r <- usum / matrix(rowSums(usum), nrow=hh, ncol=k)   #パラメータを計算

#gammaの更新
isum <- (data.frame(id=ID_i, Br=Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
gamma_r <- isum / matrix(rowSums(isum), nrow=item, ncol=k)   #パラメータを計算

#phiの更新
vf <- (data.frame(id=wd, Br=Br) %>%
         dplyr::group_by(id) %>%
         dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
phi_r <- t(vf) / matrix(colSums(vf), nrow=k, ncol=v)

#ベイズの定理よりthetaを計算
par <- beta_r[ID$id, ] * gamma_r[ID$item, ]
theta_r <- par / matrix(rowSums(par), nrow=d, ncol=k)

#対数尤度の計算
(LLS <- sum(log(rowSums(bfr$Bur))))

####EMアルゴリズムでパラメータを更新####
#更新ステータス
iter <- 1
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 0.1
LLo <- LLS   #対数尤度の初期値
LLw <- LLS


###パラメータの更新
##負担率の計算
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  ##負担率の計算
  bfr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, w=w, k=k)
  Br <- bfr$Br   #負担率
  r <- bfr$r   #混合率
  
  #betaの更新
  usum <- (data.frame(id=ID_u, Br=Br) %>%
             dplyr::group_by(id) %>%
             dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
  beta_r <- usum / matrix(rowSums(usum), nrow=hh, ncol=k)   #パラメータを計算
  
  #gammaの更新
  isum <- (data.frame(id=ID_i, Br=Br) %>%
             dplyr::group_by(id) %>%
             dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
  gamma_r <- isum / matrix(rowSums(isum), nrow=item, ncol=k)   #パラメータを計算
  
  #phiの更新
  vf <- (data.frame(id=wd, Br=Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
  phi_r <- t(vf) / matrix(colSums(vf), nrow=k, ncol=v)
  
  #ベイズの定理よりthetaを計算
  par <- beta_r[ID$id, ] * gamma_r[ID$item, ]
  theta_r <- par / matrix(rowSums(par), nrow=d, ncol=k)

  #対数尤度の計算
  (LLS <- sum(log(rowSums(bfr$Bur))))
  
  iter <- iter+1
  dl <- LLS-LLo
  LLo <- LLS
  LLw <- c(LLw, LLo)
  print(LLo)
}


