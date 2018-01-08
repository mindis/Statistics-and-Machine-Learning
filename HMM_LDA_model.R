#####HMM-LDAモデル#####
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

#set.seed(5723)

####データの発生####
##データの設定
k1 <- 8   #syntax数
k2 <- 10   #トピック数
d <- 2000   #文書数
v1 <- 300   #トピックに関係のある語彙数
v2 <- 300   #トピック以外のsyntaxに関係のある語彙数
v <- v1 + v2   #総語彙数
w <- rpois(d, rgamma(d, 100, 0.45))   #文書あたりの単語数
f <- sum(w)   #総単語数

##IDの設定
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d) {t_id <- c(t_id, 1:w[i])}

##パラメータを設定
#ディレクリ分布のパラメータ
alpha01 <- seq(3.0, 0.2, length=k1*5)[((1:(k1*5))%%5)==0]
alpha02 <- matrix(0.5, nrow=k1, ncol=k1)
alpha02[-k1, k1] <- 4.0   #k1の行および列はトピックsyntax 
alpha02[k1, k1] <- 1.25
alpha11 <- rep(0.3, k2)
alpha21 <- c(rep(0.15, v1), rep(0.0005, v2))

alloc <- as.numeric(rmnom(v2, 1, runif(k1-1, 5, 10)) %*% 1:(k1-1))
alpha22 <- cbind(matrix(0.001, nrow=k1-1, ncol=v1), matrix(0.005, nrow=k1-1, ncol=v2))
for(j in 1:(k1-1)){
  index <- which(alloc==j) + v1
  alpha22[j, index] <- 3.0              
}

#ディレクリ分布よりパラメータを生成
pi1 <- pit1 <- extraDistr::rdirichlet(1, alpha01)
pi2 <- pit2 <- extraDistr::rdirichlet(k1, alpha02)
theta <- thetat <- extraDistr::rdirichlet(d, alpha11)
phi <- phit <- extraDistr::rdirichlet(k2, alpha21)
psi <- psit <- extraDistr::rdirichlet(k1-1, alpha22)

##HMM-LDAモデルに基づき多項分布から単語を生成
wd_list <- list()
ID_list <- list()
Z1_list <- list()
Z2_list <- list()

for(i in 1:d){
  if(i%%100==0){
    print(i)
  }
  z1_vec <- rep(0, w[i])
  z2_vec <- rep(0, w[i])
  words <- rep(0, w[i])
  
  for(j in 1:w[i]){
    if(j==1){
      #文書の先頭単語のsyntaxを生成
      z1 <- rmnom(1, 1, pi1)
      z1_vec[j] <- as.numeric(z1 %*% 1:k1)
    } else {
      #先頭以降はマルコフ推移に基づきsyntaxを生成
      z1 <- rmnom(1, 1, pi2[z1_vec[j-1], ])
      z1_vec[j] <- as.numeric(z1 %*% 1:k1)
    }
  }
  
  #syntaxにもとづき単語を生成
  index_topic <- which(z1_vec==k1)
  words[-index_topic] 
  wn1 <- rmnom(w[i]-length(index_topic), 1, psi[z1_vec[-index_topic], ])
  words[-index_topic] <- as.numeric(wn1 %*% 1:v)
  
  #トピック分布を生成
  z2 <- rmnom(length(index_topic), 1, theta[i, ])
  z2_vec[index_topic] <- as.numeric(z2 %*% 1:k2)
  
  #トピック分布に基づき単語を生成
  wn2 <- rmnom(length(index_topic), 1, phi[z2_vec[index_topic], ])
  words[index_topic] <- as.numeric(wn2 %*% 1:v)
  
  #データを格納
  wd_list[[i]] <- words
  ID_list[[i]] <- rep(i, w[i]) 
  Z1_list[[i]] <- z1_vec
  Z2_list[[i]] <- z2_vec
}

#リスト形式をベクトル形式に変換
z1 <- unlist(Z1_list)
z2 <- unlist(Z2_list)
wd <- unlist(wd_list)
ID_d <- unlist(ID_list)
WX <- sparseMatrix(i=(1:f), j=wd, x=rep(1, f), dims=c(f, v))   #単語のスパース行列


##インデックスを作成
doc_list <- list()
word_list <- list()
for(i in 1:d){doc_list[[i]] <- which(ID_d==i)}
for(i in 1:v){word_list[[i]] <- which(wd==i)}


####マルコフ連鎖モンテカルロ法でHMM-LDAモデルを推定####
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


####HMM-LDAモデルのパラメータをサンプリング####
##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

#パラメータの真値
pi1 <- as.numeric(pit1)
pi2 <- pit2
theta <- thetat
phi <- phit
psi <- psit

#tfidfで初期値を設定
tf0 <- colMeans(WX)*10
idf1 <- log(nrow(WX)/colSums(WX > 0))
idf2 <- log(nrow(WX)/colSums(WX==0))

#単語トピック単位のパラメータの初期値
word_data <- matrix(0, nrow=d, ncol=v)
for(i in 1:d){
  word_data[i, ] <- colSums(WX[doc_list[[i]], ])
}
tf0 <- colSums(word_data)/sum(word_data)
idf0 <- log(d / colSums(word_data > 0))

theta <- extraDistr::rdirichlet(d, rep(0.2, k2))   #文書単位のトピックの初期値
phi <- extraDistr::rdirichlet(k2, c(rep(0.1, v1), rep(0.025, v2)))   #文書単位の出現確率の初期値
psi <- extraDistr::rdirichlet(k1-1, c(rep(0.025, v1), rep(0.3, v2)))   #文章セグメントの出現確率の初期値
pi1 <- rep(1/k1, k1)
pi2 <- extraDistr::rdirichlet(k1, c(rep(1, k1-1), 5))

##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- 0.01 
alpha02 <- 0.01
alpha11 <- 0.01
beta01 <- 1
beta02 <- 1


##パラメータの格納用配列
PI1 <- matrix(0, nrow=R/keep, ncol=k1)
PI2 <- array(0, dim=c(k1, k1, R/keep))
THETA <- array(0, dim=c(d, k2, R/keep))
PHI <- array(0, dim=c(k2, v, R/keep))
PSI <- array(0, dim=c(k1-1, v, R/keep))
SEG1 <- matrix(0, nrow=f, ncol=k1)
SEG2 <- matrix(0, nrow=f, ncol=k2)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"

##MCMC推定用配列
max_word <- max(t_id)
index_t11 <- which(t_id==1)
index_t21 <- list()
index_t22 <- list()
for(j in 2:max_word){
  index_t21[[j]] <- which(t_id==j)-1
  index_t22[[j]] <- which(t_id==j)
}

#対数尤度の基準値
LLst <- sum(WX %*% log(colSums(WX)/sum(WX)))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##単語ごとのsyntax尤度の推定とsyntaxの生成
  #syntaxごとの尤度を推定
  Li01 <- matrix(0, nrow=f, ncol=k1-1)
  for(j in 1:(k1-1)){
    Li01[, j] <- psi[j, wd]
  }
  #トピックモデルの尤度
  word_par <- burden_fr(theta, phi, wd, w, k2)   
  Li02 <- rowSums(word_par$Bur)   #トピックモデルの期待尤度
  
  #尤度の結合
  Li0 <- cbind(Li01, Li02)
  
  #マルコフ推移モデルの尤度と割当確率を逐次的に推定
  Zi1 <- matrix(0, nrow=f, ncol=k1)
  z1_rate <- matrix(0, nrow=f, ncol=k1)
  rf02 <- matrix(0, nrow=k1, ncol=k1)
  
  for(j in 1:max_word){
    if(j==1){
      #セグメントの割当確率
      Li <- matrix(pi1, nrow=length(index_t11), ncol=k1, byrow=T) * Li0[index_t11, ]   #重み付き尤度
      z1_rate[index_t11, ] <- Li / rowSums(Li)   #割当確率
      
      #多項分布よりセグメントを生成
      Zi1[index_t11, ] <- rmnom(length(index_t11), 1, z1_rate[index_t11, ])
      z1_vec[index_t11] <- as.numeric(Zi1[index_t11, ] %*% 1:k1)
      
      #混合率のパラメータを更新
      rf01 <- colSums(Zi1[index_t11, ])
      
    } else {
      
      #セグメントの割当確率
      index <- index_t22[[j]]
      Li <- pi2[z1_vec[index_t21[[j]]], , drop=FALSE] * Li0[index, , drop=FALSE]   #重み付き尤度
      z1_rate[index, ] <- Li / rowSums(Li)   #割当確率
      
      #多項分布よりセグメントを生成
      Zi1[index, ] <- rmnom(length(index), 1, z1_rate[index, ])
      z1_vec[index] <- as.numeric(Zi1[index, ] %*% 1:k1)
      
      #混合率のパラメータを更新
      rf02 <- rf02 + t(Zi1[index_t21[[j]], , drop=FALSE]) %*% Zi1[index, , drop=FALSE]   #マルコフ推移
    }
  }
  index_topic <- which(z1_vec==k1)
  n <- length(index_topic)
  
  #ディクレリ分布からHMMの混合率をサンプリング
  rf11 <- rf01 + alpha01
  rf12 <- rf02 + alpha01
  pi1 <- extraDistr::rdirichlet(1, rf11)
  pi2 <- extraDistr::rdirichlet(k1, rf12)
  
  #ディクレリ分布からsyntaxのパラメータをサンプリング
  Zi1_syntax <- Zi1[-index_topic, ]
  WX_syntax <- WX[-index_topic, ]
  df0 <- matrix(0, nrow=k1-1, ncol=v)
  for(j in 1:(k1-1)){
    df0[j, ] <- colSums(WX_syntax * Zi1_syntax[, j])
  }
  df <- df0 + alpha02 
  psi <- extraDistr::rdirichlet(k1-1, df)
  
  ##単語ごとにトピックをサンプリング
  Zi2 <- matrix(0, nrow=f, ncol=k2)
  Zi2[index_topic, ] <- rmnom(n, 1, word_par$Br[index_topic, ])   #多項分布よりトピックをサンプリング
  z2_vec <- as.numeric(Zi2 %*% 1:k2)
  
  
  ##トピックモデルのパラメータをサンプリング
  #トピック分布thetaをサンプリング
  wsum0 <- matrix(0, nrow=d, ncol=k2)
  for(i in 1:d){
    wsum0[i, ] <- colSums(Zi2[doc_list[[i]], ])
  }
  wsum <- wsum0 + beta01
  theta <- extraDistr::rdirichlet(d, wsum)
  
  #単語分布phiをサンプリング
  vf0 <- matrix(0, nrow=k2, ncol=v)
  for(j in 1:v){
    vf0[, j] <- colSums(Zi2[word_list[[j]], ])
  }
  vf <- vf0 + alpha11
  phi <- extraDistr::rdirichlet(k2, vf)
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    PI1[mkeep, ] <- pi1
    PI2[, , mkeep] <- pi2
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    PSI[, , mkeep] <- psi
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(mkeep >= burnin & rp%%keep==0){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
    }
    
    #サンプリング結果を確認
    if(rp%%disp==0){
      LL1 <- sum(log(rowSums(word_par$Bur[index_topic, ] * Zi2[index_topic, ])))
      LL2 <- sum(log(rowSums(Li01[-index_topic, ] * Zi1[-index_topic, -k1])))
      print(rp)
      print(c(LL1 + LL2, LLst))
      print(round(cbind(pi2, pit2), 3))
      print(round(cbind(theta[1:6, ], thetat[1:6, ]), 3))
      print(round(cbind(psi[, 296:305], psit[, 296:305]), 3))
    }
  }
}


####サンプリング結果の可視化と要約####
burnin <- 1000/keep   #バーンイン期間
RS <- R/keep

##サンプリング結果の可視化
#文書のトピック分布のサンプリング結果
matplot(t(THETA[1, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA[100, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA[1000, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA[2000, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")

#単語の出現確率のサンプリング結果
matplot(t(PHI[, 1, ]), type="l", ylab="パラメータ", main="トピック1の単語の出現率のサンプリング結果")
matplot(t(PHI[, 200, ]), type="l", ylab="パラメータ", main="トピック2の単語の出現率のサンプリング結果")
matplot(t(PHI[, 400, ]), type="l", ylab="パラメータ", main="トピック3の単語の出現率のサンプリング結果")
matplot(t(PHI[, 500, ]), type="l", ylab="パラメータ", main="トピック4の単語の出現率のサンプリング結果")
matplot(t(PSI[, 1, ]), type="l", ylab="パラメータ", main="トピック1の単語の出現率のサンプリング結果")
matplot(t(PSI[, 200, ]), type="l", ylab="パラメータ", main="トピック2の単語の出現率のサンプリング結果")
matplot(t(PSI[, 400, ]), type="l", ylab="パラメータ", main="トピック3の単語の出現率のサンプリング結果")
matplot(t(PSI[, 500, ]), type="l", ylab="パラメータ", main="トピック4の単語の出現率のサンプリング結果")


##サンプリング結果の要約推定量
#トピック分布の事後推定量
topic_mu <- apply(THETA[, , burnin:(R/keep)], c(1, 2), mean)   #トピック分布の事後平均
round(cbind(topic_mu, thetat), 3)
round(topic_sd <- apply(THETA[, , burnin:(R/keep)], c(1, 2), sd), 3)   #トピック分布の事後標準偏差

#単語出現確率の事後推定量
word_mu1 <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #単語の出現率の事後平均
word1 <- round(t(rbind(word_mu1, phit)), 3)
word_mu2 <- apply(OMEGA[, , burnin:(R/keep)], c(1, 2), mean)   #単語の出現率の事後平均
word2 <- round(t(rbind(word_mu2, omegat)), 3)
word <- round(t(rbind(word_mu1, word_mu2, phit, omegat)), 3)
colnames(word) <- 1:ncol(word)

word_mu3 <- apply(GAMMA[burnin:(R/keep), ], 2, mean)   #単語の出現率の事後平均
round(rbind(word_mu3, gamma=gammat), 3)


##トピックの事後分布の要約
round(seg1_mu <- SEG1 / rowSums(SEG1), 3)
round(seg2_mu <- SEG2 / rowSums(SEG2), 3)
