#####パチンコ配分モデル#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
detach("package:gtools", unload=TRUE)
library(bayesm)
library(ExtDist)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(76432)

####データの発生####
#set.seed(423943)
#文書データの設定
d <- 2500   #文書数
v <- 250   #語彙数
w <- rpois(d, rgamma(d, 160, 1.0))   #1文書あたりの単語数
hist(w, col="grey")

#トピックを設定
k1 <- 3   #上位トピック数
k2 <- 11   #下位トピック総数
k21 <- 5   #下位トピック数1
k22 <- 4   #下位トピック数2
k23 <- 4   #下位トピック数3
k2v <- c(k21, k22, k23)

#ネスト構造を設定
nest <- rbind(c(rep(1, k21), rep(0, k2-k21)),
              c(rep(0, 4), rep(1, k22), rep(0, 3)),
              c(rep(0, k2-k23), rep(1, k23)))


#IDの設定
word_id <- rep(1:d, w)

##パラメータの設定
alpha01 <- runif(k1, 0.4, 0.7)   #文書の上位トピックのディレクリ事前分布のパラメータ
alpha021 <- rep(0.25, k21)   #文書の下位トピック1のディレクリ事前分布のパラメータ
alpha022 <- rep(0.3, k22)   #文書の下位トピック2のディレクリ事前分布のパラメータ
alpha023 <- rep(0.3, k23)   #文書の下位トピック2のディレクリ事前分布のパラメータ
alpha1 <- rep(0.15, v)   #単語のディレクリ事前分布のパラメータ

#ディレクリ乱数の発生
theta01 <- extraDistr::rdirichlet(d, alpha01)
theta021 <- extraDistr::rdirichlet(d, alpha021)
theta022 <- extraDistr::rdirichlet(d, alpha022)
theta023 <- extraDistr::rdirichlet(d, alpha023)
theta02 <- list(theta021, theta022, theta023)
phi0 <- extraDistr::rdirichlet(k2, alpha1)   #単語のトピック分をディレクリ分布から発生

##多項分布からトピックおよび単語データを発生
WX <- matrix(0, nrow=d, ncol=v)
Z1 <- list()
Z2 <- list()

#文書ごとにトピックと単語を逐次生成
for(i in 1:d){
  print(i)
  
  #文書の上位トピック分布を発生
  z1 <- t(rmultinom(w[i], 1, theta01[i, ]))   #文書の上位トピック分布を発生
  z1_vec <- as.numeric(z1 %*% 1:k1)
  
  #文書の下位トピック分布を発生
  z2 <- matrix(0, nrow=length(z1_vec), ncol=k2)
  for(j in 1:length(z1_vec)){
    z2[j, nest[z1_vec[j], ]==1] <- as.numeric(t(rmultinom(1, 1, theta02[[z1_vec[j]]][i, ])))
  }
  z2
  #文書の下位トピック分布から単語を発生
  zn <- z2 %*% c(1:k2)   #0,1を数値に置き換える
  zdn <- cbind(zn, z2)   #apply関数で使えるように行列にしておく
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi0[x[1], ])))   #文書のトピックから単語を生成
  wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
  WX[i, ] <- wdn  
  
  #発生させたトピックを格納
  Z1[[i]] <- z1_vec
  Z2[[i]] <- zdn[, 1]
}

#発生させたデータの確認
colSums(WX); round(colMeans(WX), 3)
table(unlist(Z1)); table(unlist(Z2))
storage.mode(WX) <- "integer"   #データ行列を整数型行列に変更


####トピックモデルのためのデータと関数の準備####
##それぞれの文書中の単語の出現および補助情報の出現をベクトルに並べる
##データ推定用IDを作成
ID_list <- list()
wd_list <- list()

#求人ごとに求人IDおよび単語IDを作成
for(i in 1:nrow(WX)){
  print(i)
  
  #単語のIDベクトルを作成
  ID_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
  num2 <- subset(num1, num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
}

#リストをベクトルに変換
ID_d <- unlist(ID_list)
wd <- unlist(wd_list)

##インデックスを作成
doc_list <- list()
word_list <- list()
for(i in 1:length(unique(ID_d))) {doc_list[[i]] <- which(ID_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- which(wd==i)}
gc(); gc()


####ネストを識別するための初期値を設定する####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, k){
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

##EMアルゴリズムの初期値を設定する
##初期値をランダムに設定
#phiの初期値
freq_v <- matrix(colSums(WX), nrow=k2, ncol=v, byrow=T)   #単語の出現率
rand_v <- matrix(trunc(rnorm(k2*v, 0, (colSums(WX)/2))), nrow=k2, ncol=v, byrow=T)   #ランダム化
phi_r <- abs(freq_v + rand_v) / rowSums(abs(freq_v + rand_v))   #トピックごとの出現率をランダムに初期化

#thetaの初期値
theta_r <- rdirichlet(d, runif(k2, 0.2, 4))   #ディレクリ分布から初期値を設定


###パラメータの更新
##負担率の計算
bfr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, k=k2)
Br <- bfr$Br   #負担率
r <- bfr$r   #混合率

##thetaの更新
tsum <- (data.frame(id=ID_d, Br=Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_all(funs(sum)))[, 2:(k2+1)]
theta_r <- tsum / matrix(w, nrow=d, ncol=k2, byrow=T)   #パラメータを計算

##phiの更新
vf <- (data.frame(id=wd, Br=Br) %>%
         dplyr::group_by(id) %>%
         dplyr::summarize_all(funs(sum)))[, 2:(k2+1)]
phi_r <- t(vf) / matrix(colSums(vf), nrow=k2, ncol=v)

#対数尤度の計算
(LLS <- sum(log(rowSums(bfr$Bur))))


####EMアルゴリズムでパラメータを更新####
#更新ステータス
iter <- 1
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 0.5
LLo <- LLS   #対数尤度の初期値
LLw <- LLS

###パラメータの更新
##負担率の計算
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  bfr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, k=k2)
  Br <- bfr$Br   #負担率
  r <- bfr$r   #混合率
  
  ##thetaの更新
  tsum <- (data.frame(id=ID_d, Br=Br) %>%
             dplyr::group_by(id) %>%
             dplyr::summarize_all(funs(sum)))[, 2:(k2+1)]
  theta_r <- tsum / matrix(w, nrow=d, ncol=k2, byrow=T)   #パラメータを計算
  
  ##phiの更新
  vf <- (data.frame(id=wd, Br=Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_all(funs(sum)))[, 2:(k2+1)]
  phi_r <- t(vf) / matrix(colSums(vf), nrow=k2, ncol=v)
  
  #対数尤度の計算
  LLS <- sum(log(rowSums(bfr$Bur)))
  
  iter <- iter+1
  dl <- LLS-LLo
  LLo <- LLS
  LLw <- c(LLw, LLo)
  print(LLo)
}


####マルコフ連鎖モンテカルロ法でパチンコ配分モデルを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta1, theta2, phi, wd, w, k1, k2, nest){
  Bur <- matrix(0, nrow=length(wd), ncol=k2)   #負担係数の格納用配列
  for(j in 1:k1){
    for(k in 1:k2){
      if(nest[j, k]==0) next
      theta_nest <- theta2[[j]]
      Bur[, k] <- Bur[, k] + rep(theta1[, j], w) * rep(theta_nest[, sum(nest[j, 1:k]==1)], w) * phi[k, wd]
    }
  }
  Br <- Bur / rowSums(Bur)   #負担率の計算
  r <- colSums(Br) / sum(Br)   #混合率の計算
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}


##アルゴリズムの設定
R <- 10000   #サンプリング回数
keep <- 2   #2回に1回の割合でサンプリング結果を格納
iter <- 0

##ネストの設定
sort_vec <- rep(0, k2)
for(i in 1:k2){
  er_sq <- rep(0, k2)
  for(j in 1:k2){
    er_sq[j] <- sum((phi_r[j, ] - phi0[i, ])^2)
  }
  er_sq
  sort_vec[i] <- which.min(er_sq)
}


##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- rep(1.0, k1)
alpha02 <- rep(1.0, k2)
alpha01m <- matrix(alpha01, nrow=d, ncol=k1, byrow=T)
alpha02m <- matrix(alpha02, nrow=d, ncol=k2, byrow=T)
beta0 <- rep(0.5, v)
beta0m <- matrix(beta0, nrow=v, ncol=k2)


##パラメータの初期値
#上位トピックおよび下位トピックの初期パラメータを設定
Zi0 <- Br[, sort_vec]
theta_ho <- matrix(0, nrow=length(wd), ncol=k1)
wsum1 <- matrix(0, nrow=d, ncol=k1)
theta2 <- list()
theta1 <- matrix(1/k1, nrow=d, ncol=k1)   #文書の上位トピックのパラメータの初期値

for(j in 1:k1){
  theta2[[j]] <- (theta_r[, sort_vec])[, nest[j, ]==1]
  theta_ho[, j] <- rowSums((matrix(theta1[, j], nrow=d, ncol=sum(nest[j, ])) * theta2[[j]])[ID_d, ] * Zi0[, nest[j, ]==1])
}

theta_rate <- theta_ho / rowSums(theta_ho)   #所属確率

#多項分布から上位トピックをサンプリング
Zi1 <- rmnom(nrow(theta_rate), 1, theta_rate)

#上位トピック分布のパラメータを更新
for(i in 1:d){wsum1[i, ] <- colSums(Zi1[doc_list[[i]], ])}
theta1 <- extraDistr::rdirichlet(d, wsum1 + alpha01m)

#単語トピックのパラメータの初期値
phi <- phi_r[sort_vec, ]   


##パラメータの格納用配列
THETA1 <- array(0, dim=c(d, k1, R/keep))
THETA21 <- array(0, dim=c(d, k21, R/keep))
THETA22 <- array(0, dim=c(d, k22, R/keep))
THETA23 <- array(0, dim=c(d, k23, R/keep))
PHI <- array(0, dim=c(k2, v, R/keep))
Z0_SEG <- matrix(0, nrow=length(wd), ncol=k2)
Z1_SEG <- matrix(0, nrow=length(wd), ncol=k1)
storage.mode(Z0_SEG) <- "integer"
storage.mode(Z1_SEG) <- "integer"
gc(); gc()


##MCMC推定用配列
wsum0 <- matrix(0, nrow=d, ncol=k2)
wsum1 <- matrix(0, nrow=d, ncol=k1)
wsum2 <- list()
for(j in 1:k1){
  wsum2[[j]] <- matrix(0, nrow=d, ncol=sum(nest[j, ]))
}
theta_ho <- matrix(0, nrow=length(ID_d), ncol=k1)
vf0 <- matrix(0, nrow=v, ncol=k2)


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##単語ごとにトピックをサンプリング
  #単語ごとにトピックの出現率を計算
  fr <- burden_fr(theta1, theta2, phi, wd, w, k1, k2, nest)
  word_rate <- fr$Br
  
  #多項分布から単語トピックをサンプリング
  Zi0 <- rmnom(nrow(word_rate), 1, word_rate)
  
  #文書ごとにトピック出現確率を計算
  for(i in 1:d){wsum0[i, ] <- colSums(Zi0[doc_list[[i]], ])}
  
  
  ##上位トピックパラメータをサンプリング
  #下位トピックの出現率を周辺化した負担率
  for(j in 1:k1){
    theta_ho[, j] <- rowSums((matrix(theta1[, j], nrow=d, ncol=sum(nest[j, ])) * theta2[[j]])[ID_d, ] * Zi0[, nest[j, ]==1])
  }
  theta_rate <- theta_ho / rowSums(theta_ho)
  
  #多項分布から上位トピックをサンプリング
  Zi1 <- rmnom(nrow(theta_rate), 1, theta_rate)
  
  #上位トピック分布のパラメータを更新
  for(i in 1:d){wsum1[i, ] <- colSums(Zi1[doc_list[[i]], ])}
  theta1 <- extraDistr::rdirichlet(d, wsum1 + alpha01m)
  
  ##下位トピックパラメータをサンプリング
  for(j in 1:k1){
    for(i in 1:d){
      #ディクレリ分布のパラメータ
      wsum2[[j]][i, ] <- colSums(Zi0[doc_list[[i]], nest[j, ]==1] * 
                                   matrix(Zi1[doc_list[[i]], j], nrow=length(doc_list[[i]]), ncol=sum(nest[j, ]))) + 1
    }
    
    #下位トピック分布のパラメータを更新
    theta2[[j]] <- extraDistr::rdirichlet(d, wsum2[[j]])
  }
  
  ##トピックごとに単語分布のパラメータをサンプリング
  for(i in 1:v){vf0[i, ] <- colSums(Zi0[word_list[[i]], ])}
  vf <- vf0 + beta0m
  phi <- extraDistr::rdirichlet(k2, t(vf))
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA1[, , mkeep] <- theta1
    THETA21[, , mkeep] <- theta2[[1]]
    THETA22[, , mkeep] <- theta2[[2]]
    THETA23[, , mkeep] <- theta2[[3]]
    PHI[, , mkeep] <- phi
    
    #トピック割当はサンプリング期間の半分を超えたら格納する
    if(rp >= R/2){
      Z0_SEG <- Z0_SEG + Zi0
      Z1_SEG <- Z1_SEG + Zi1
    }
    
    #サンプリング結果を確認
    print(rp)
    print(round(cbind(theta1[1:12, ], theta01[1:12, ], wsum2[[1]][1:12, ], theta2[[1]][1:12, ], theta02[[1]][1:12, ]), 3))
    print(round(cbind(phi[, 1:10], phi0[, 1:10]), 3))
  }
}


####サンプリング結果の可視化と要約####
burnin <- 2000   #バーンイン期間

##サンプリング結果の可視化
#文書の上位トピック分布のサンプリング結果
matplot(t(THETA1[1, , ]), type="l", ylab="パラメータ", main="文書1の上位トピック分布のサンプリング結果")
matplot(t(THETA1[2, , ]), type="l", ylab="パラメータ", main="文書2の上位トピック分布のサンプリング結果")
matplot(t(THETA1[3, , ]), type="l", ylab="パラメータ", main="文書3の上位トピック分布のサンプリング結果")
matplot(t(THETA1[4, , ]), type="l", ylab="パラメータ", main="文書4の上位トピック分布のサンプリング結果")

#文書の下位トピック分布のサンプリング結果
matplot(t(THETA21[1, , ]), type="l", ylab="パラメータ", main="文書1の下位トピック分布のサンプリング結果")
matplot(t(THETA21[2, , ]), type="l", ylab="パラメータ", main="文書2の下位トピック分布のサンプリング結果")
matplot(t(THETA21[3, , ]), type="l", ylab="パラメータ", main="文書3の下位トピック分布のサンプリング結果")
matplot(t(THETA21[4, , ]), type="l", ylab="パラメータ", main="文書4の下位トピック分布のサンプリング結果")

#単語の出現確率のサンプリング結果
matplot(t(PHI[1, 1:10, ]), type="l", ylab="パラメータ", main="トピック1の単語の出現率のサンプリング結果")
matplot(t(PHI[2, 11:20, ]), type="l", ylab="パラメータ", main="トピック2の単語の出現率のサンプリング結果")
matplot(t(PHI[3, 21:30, ]), type="l", ylab="パラメータ", main="トピック3の単語の出現率のサンプリング結果")
matplot(t(PHI[4, 31:40, ]), type="l", ylab="パラメータ", main="トピック4の単語の出現率のサンプリング結果")
matplot(t(PHI[5, 41:50, ]), type="l", ylab="パラメータ", main="トピック5の単語の出現率のサンプリング結果")
matplot(t(PHI[6, 51:60, ]), type="l", ylab="パラメータ", main="トピック6の単語の出現率のサンプリング結果")
matplot(t(PHI[7, 61:70, ]), type="l", ylab="パラメータ", main="トピック7の単語の出現率のサンプリング結果")
matplot(t(PHI[8, 71:80, ]), type="l", ylab="パラメータ", main="トピック8の単語の出現率のサンプリング結果")


##サンプリング結果の要約推定量
#上位トピック分布の事後推定量
topic_mu1 <- apply(THETA1[, , burnin:(R/keep)], c(1, 2), mean)   #トピック分布の事後平均
round(cbind(topic_mu1, theta01), 3)
round(topic_sd1 <- apply(THETA1[, , burnin:(R/keep)], c(1, 2), sd), 3)   #トピック分布の事後標準偏差

#下位トピック分布の事後推定量
topic_mu21 <- apply(THETA21[, , burnin:(R/keep)], c(1, 2), mean)   #トピック分布の事後平均
topic_mu22 <- apply(THETA22[, , burnin:(R/keep)], c(1, 2), mean)   #トピック分布の事後平均
topic_mu23 <- apply(THETA23[, , burnin:(R/keep)], c(1, 2), mean)   #トピック分布の事後平均
round(cbind(topic_mu21, theta021), 3)
round(cbind(topic_mu22, theta022), 3)
round(cbind(topic_mu23, theta023), 3)
round(topic_sd21 <- apply(THETA21[, , burnin:(R/keep)], c(1, 2), sd), 3)   #トピック分布の事後標準偏差
round(topic_sd22 <- apply(THETA22[, , burnin:(R/keep)], c(1, 2), sd), 3)   #トピック分布の事後標準偏差
round(topic_sd23 <- apply(THETA23[, , burnin:(R/keep)], c(1, 2), sd), 3)   #トピック分布の事後標準偏差

#単語出現確率の事後推定量
word_mu <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #単語の出現率の事後平均
round(rbind(word_mu, phi0)[, 1:50], 3)



  