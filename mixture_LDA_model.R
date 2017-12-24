#####混合LDAモデル#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(2578)

####データの発生####
##データの設定
k1 <- 10
k2 <- 4
d <- 2000   #文書数
v1 <- 250   #内容に関係のある語彙数
v2 <- 100   #内容に関係のない語彙数
v <- v1 + v2   #語彙数 
s <- rpois(d, 14)   #文章数
s[s < 5] <- ceiling(runif(sum(s < 5), 5, 10))
a <- sum(s)   #総文章数
w <- rpois(a, 11.5)   #文章あたりの単語数
w[w < 5] <- ceiling(runif(sum(w < 5), 5, 10))
f <- sum(w)   #総単語数

#文書IDの設定
u_id <- rep(1:d, s)
t_id <- c()
for(i in 1:d){t_id <- c(t_id, 1:s[i])}
words <- as.numeric(tapply(w, u_id, sum))


##パラメータを設定
#ディレクリ分布のパラメータ
alpha0 <- rep(0.15, k1)
alpha1 <- c(rep(0.5, v1), rep(0.005, v2))
alpha2 <- rep(1, k2)
alpha3 <- c(rep(0.005, v1), rep(0.5, v2))

#ディレクリ分布よりパラメータを生成
thetat <- theta <- extraDistr::rdirichlet(d, alpha0)
phit <- phi <- extraDistr::rdirichlet(k1, alpha1)
omegat <- omega <- matrix(1/k2, nrow=d, ncol=k2)
gammat <- gamma <- extraDistr::rdirichlet(k2, alpha3)
delta <- deltat <- rbeta(d, 20, 15)


##文章ごとに単語を生成する
WX <- matrix(0, nrow=a, ncol=v)
y_list <- list()
Z1_list <- list()
Z2_list <- list()
index_v1 <- 1:v1
index_v2 <- (v1+1):v 

for(i in 1:a){
  
  ##文章ごとにトピックを生成
  id <- u_id[i]
  z1 <- rmnom(w[i], 1, theta[id, ])
  z1_vec <- as.numeric(z1 %*% 1:k1)
  
  ##一般語の判定と一般語トピックの生成
  #一般語かどうかの生成
  y0 <- rbinom(w[i], 1, delta[id])
  index <- which(y0==1)
  
  #文章のトピックを生成
  z2 <- rmnom(1, 1, omega[id, ])
  z2_vec <- as.numeric(z2 %*% 1:k2)
  
  #トピック分布に基づき単語を生成
  n1 <- length(index)
  n2 <- w[i]-length(index)
  w2 <- w1 <- matrix(0, nrow=1, ncol=v)
  
  if(n1 > 0){
    w1 <- rmnom(n1, 1, phi[z1_vec[index], ])   #トピック語を生成
  }
  if(n2 > 0){
    w2 <- rmnom(n2, 1, gamma[z2_vec, ])   #一般語を生成
  }
  
  #bag of words行列を作成
  wdn <- colSums(w1) + colSums(w2)
  WX[i, ] <- wdn
  
  #パラメータを格納
  Z1_list[[i]] <- z1
  Z2_list[[i]] <- z2
  y_list[[i]] <- y0 
}


#リスト形式を変換
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)
y_vec <- unlist(y_list)


####トピックモデル推定のためのデータと関数の準備####
##それぞれの文書中の単語の出現および補助情報の出現をベクトルに並べる
##データ推定用IDを作成
ID_list <- list()
td1_list <- list()
td2_list <- list()
wd_list <- list()

#文書idごとに文章idおよび単語idを作成
for(i in 1:a){
  
  #文書IDを記録
  ID_list[[i]] <- rep(u_id[i], w[i])
  td1_list[[i]] <- rep(i, w[i])
  td2_list[[i]] <- rep(t_id[i], w[i])
  
  #単語IDを記録
  num1 <- WX[i, ] * 1:v
  num2 <- which(num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
}

#リストをベクトルに変換
ID_d <- unlist(ID_list)
td1_d <- unlist(td1_list)
td2_d <- unlist(td2_list)
wd <- unlist(wd_list)

##インデックスを作成
doc_list <- list()
id_list <- list()
sent_list <- list()
word_list <- list()
for(i in 1:length(unique(ID_d))){doc_list[[i]] <- which(ID_d==i)}
for(i in 1:d){id_list[[i]] <- which(u_id==i)}
for(i in 1:length(unique(td1_d))){sent_list[[i]] <- which(td1_d==i)}
for(i in 1:length(unique(wd))){word_list[[i]] <- which(wd==i)}


####マルコフ連鎖モンテカルロ法で混合LDAモデルを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #負担係数の格納用
  for(j in 1:k){
    #負担係数を計算
    Bi <- rep(theta[, j], w) * phi[j, wd]   #尤度
    Bur[, j] <- Bi   
  }
  
  Br <- Bur / rowSums(Bur)   #負担率の計算
  r <- colSums(Br) / sum(Br)   #混合率の計算
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}


##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep

##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- 1  
beta01 <- 0.5
gamma01 <- 0.5 

##パラメータの初期値
#tfidfで初期値を設定
tf <- WX/rowSums(WX)
idf1 <- log(nrow(WX)/colSums(WX > 0))
idf2 <- log(nrow(WX)/colSums(WX==0))

#単語トピック単位のパラメータの初期値
theta <- extraDistr::rdirichlet(d, rep(1, k1))   #ユーザートピックの初期値
phi <- extraDistr::rdirichlet(k1, idf1*10)   #評価対象語の出現確率の初期値

#一般語トピック単位のパラメータの初期値
omega <- 1/k2
gamma <- extraDistr::rdirichlet(k2, idf2*100)   #一般語の出現確率の初期値
y <- rbinom(f, 1, 0.5)
r <- rep(0.5, f)


##パラメータの格納用配列
THETA <- array(0, dim=c(d, k1, R/keep))
PHI <- array(0, dim=c(k1, v, R/keep))
OMEGA <- array(0, dim=c(d, k2, R/keep))
GAMMA <- array(0, dim=c(k2, v, R/keep))
SEG1 <- matrix(0, nrow=f, ncol=k1)
SEG2 <- matrix(0, nrow=a, ncol=k2)
Y <- matrix(0, nrow=f, ncol=2)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"

##MCMC推定用配列
wsum0 <- matrix(0, nrow=d, ncol=k1)
vf0 <- matrix(0, nrow=k1, ncol=v)
dsum0 <- matrix(0, nrow=d, ncol=k2)
sf0 <- matrix(0, nrow=k2, ncol=v)


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##トピック関連語ごとにトピック分布のパラメータを推定
  #トピックの出現率と尤度を計算
  word_par1 <- burden_fr(theta, phi, wd, words, k1)
  word_rate1 <- word_par1$Br   #トピック割当確率
  LH1 <- rowSums(word_par1$Bur)   #尤度
  
  ##一般語かどうかをサンプリング
  #一般語のトピックの尤度と割当確率を計算
  word_par2 <- matrix(0, nrow=f, ncol=k2)
  for(j in 1:k2){
    word_par2[, j] <- omega * gamma[j, wd]
  }
  LH2 <- rowSums(word_par2)   #一般語の尤度
  
  
  #二項分布より一般語かどうかをサンプリング
  y_rate <- r*LH1 / (r*LH1 + (1-r)*LH2)   #潜在変数の割当確率
  y <- rbinom(f, 1, y_rate)   #二項分布より潜在変数をサンプリング
  index <- which(y==1)

  #混合率の更新
  r <- 0.5
  
  ##生成したスイッチング変数に基づきトピックを生成
  #単語単位のトピックをサンプリング
  Zi1 <- matrix(0, nrow=f, ncol=k1)
  zi1_vec <- rep(0, f)
  Zi1[index, ] <- rmnom(length(index), 1, word_rate1[index, ])
  zi1_vec[index] <- as.numeric(Zi1[index, ] %*% 1:k1)
  
  
  #文章単位のトピックをサンプリング
  LH2 <- word_par2 * 10^10   #尤度を桁落ちしないように定数をかける
  LH2[index, ] <- 1
  LL <- matrix(0, nrow=a, ncol=k2)
  for(i in 1:a){
    LL[i, ] <- colProds(LH2[sent_list[[i]], ])
  }
  index_ones <- which(LL[, 1]==1)
  
  
  #多項分布より文章単位のトピックをサンプリング
  sentence_rate <- LL / rowSums(LL)
  Zi2 <- rmnom(a, 1, sentence_rate)
  Zi2[index_ones, ] <- 0
  zi2_vec <- as.numeric(Zi2 %*% 1:k2)
  
  
  #トピックの混合率omegaをサンプリング
  #for(i in 1:d){
  #  dsum0[i, ] <- colSums(Zi2[sent_list[[i]], ])
  #}
  #dsum <- dsum0 + alpha01
  #omega <- extraDistr::rdirichlet(d, dsum)
  
  #単語分布をサンプリング
  zi2_word <- rep(zi2_vec, w)[-index]
  wd0 <- wd[-index]
  
  for(j in 1:k2){
    index_seg <- which(zi2_word==j)
    sf0[j, ] <- plyr::count(c(wd0[index_seg], 1:v))$freq
  }
  gamma <- extraDistr::rdirichlet(k2, sf0)
  
  
  ##単語単位のパラメータをサンプリング
  #トピック分布thetaをサンプリング
  for(i in 1:d){
    wsum0[i, ] <- colSums(Zi1[doc_list[[i]], ])
  }
  wsum <- wsum0 + alpha01
  theta <- extraDistr::rdirichlet(d, wsum)
  
  
  #単語分布phiをサンプリング
  for(j in 1:v){
    vf0[, j] <- colSums(Zi1[word_list[[j]], ] )
  }
  vf <- vf0 + beta01
  phi <- extraDistr::rdirichlet(k1, vf)
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    #OMEGA[, , mkeep] <- omega
    GAMMA[, , mkeep] <- gamma
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(mkeep >= burnin & rp%%keep==0){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
      Y <- Y + y
    }
    
    #サンプリング結果を確認
    print(rp)
    print(c(mean(y), mean(y_vec)))
    print(round(cbind(theta[1:7, ], thetat[1:7, ]), 3))
    print(round(cbind(phi[, 246:255], phit[, 246:255]), 3))
    print(round(cbind(gamma[, 246:255], gammat[, 246:255]), 3))
  }
}

