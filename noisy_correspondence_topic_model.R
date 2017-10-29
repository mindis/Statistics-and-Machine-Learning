#####ノイズあり対応トピックモデル#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
detach("package:gtools", unload=TRUE)
detach("package:bayesm", unload=TRUE)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(8079)

####データの発生####
#set.seed(423943)
#文書データの設定
k <- 8   #トピック数
d <- 2000   #文書数
v <- 300   #語彙数
w <- rpois(d, 200)   #1文書あたりの単語数
a1 <- 40   #トピックに関係のあるタグ数
a2 <- 20   #トピックに関係のないタグ数
a <- a1 + a2   #補助変数数
x0 <- rpois(d, 25)
x <- ifelse(x0 < 1, 1, x0)

#IDの設定
word_id <- rep(1:d, w)
aux_id <- rep(1:d, x)

#パラメータの設定
alpha0 <- round(runif(k, 0.1, 1.25), 3)   #文書のディレクリ事前分布のパラメータ
alpha1 <- rep(0.25, v)   #単語のディレクリ事前分布のパラメータ
alpha2 <- c(rep(0.3, a1), rep(0.025, a2))   #トピックに関係のあるタグのディクレリ事前分布のパラメータ
alpha3 <- c(rep(0.05, a1), rep(1, a2))   #トピックに関係のないタグのディクレリ事前分布のパラメータ
beta0 <- rbeta(sum(x), 0.55, 0.15)

#ディレクリ乱数の発生
thetat <- theta <- rdirichlet(d, alpha0)   #文書のトピック分布をディレクリ乱数から発生
phit <- phi <- rdirichlet(k, alpha1)   #単語のトピック分布をディレクリ乱数から発生
lambda <- matrix(0, nrow=d, ncol=k)   #文書に含むトピックだけを補助情報のトピックにするための確率を格納する行列
omegat <- omega <- rdirichlet(k, alpha2)   #補助情報のトピック分布をディクレリ乱数から発生
gammat <- gamma <- rdirichlet(1, alpha3)   #トピックに関係のないタグ
omega0 <- rbind(omega, gamma)   #単語分布の結合


##多項分布からトピックおよび単語データを発生
WX <- matrix(0, nrow=d, ncol=v)
AX <- matrix(0, nrow=d, ncol=a)
z0 <- rep(0, sum(x)) 
Z1 <- list()
Z2 <- list()

#文書ごとにトピックと単語を逐次生成
for(i in 1:d){
  print(i)
  
  #文書のトピック分布を発生
  z1 <- t(rmultinom(w[i], 1, theta[i, ]))   #文書のトピック分布を発生
  
  #文書のトピック分布から単語を発生
  zn <- z1 %*% c(1:k)   #0,1を数値に置き換える
  zdn <- cbind(zn, z1)   #apply関数で使えるように行列にしておく
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi[x[1], ])))   #文書のトピックから単語を生成
  wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
  WX[i, ] <- wdn  
  
  #文書のトピック分布から補助変数を発生
  #文書で発生させたトピックのみを補助情報のトピック分布とする
  rate <- rep(0, k)
  z_table <- table(zn)
  lambda[i, as.numeric(names(z_table))] <- z_table/sum(z_table)   #補助情報のトピック分布
  
  #ベルヌーイ分布からトピックに関係があるかどうかを発生
  index <- which(aux_id==i) 
  z0[index] <- rbinom(length(index), 1, beta0[index])
  
  #トピックを発生
  z02 <- t(rmultinom(x[i], 1, lambda[i, ]))  
  z2 <- cbind(z02 * matrix(z0[index], nrow=length(index), ncol=k), 1-z0[index])
  
  #発生させたトピックの単語分布に従い単語を発生
  zx <- z2 %*% 1:(k+1)
  zax <- cbind(zx, z2)
  an <- t(apply(zax, 1, function(x) rmultinom(1, 1, omega0[x[1], ])))
  adn <- colSums(an)
  AX[i, ] <- adn
  
  #文書トピックおよび補助情報トピックを格納
  Z1[[i]] <- zdn[, 1]
  Z2[[i]] <- zax[, 1]
}

#データ行列を整数型行列に変更
storage.mode(WX) <- "integer"
storage.mode(AX) <- "integer"
r0 <- c(mean(z0), 1-mean(z0))


####EMアルゴリズムでトピックモデルを推定####
####トピックモデルのためのデータと関数の準備####
##それぞれの文書中の単語の出現および補助情報の出現をベクトルに並べる
##データ推定用IDを作成
ID1_list <- list()
wd_list <- list()
ID2_list <- list()
ad_list <- list()

#求人ごとに求人IDおよび単語IDを作成
for(i in 1:nrow(WX)){
  print(i)
  
  #単語のIDベクトルを作成
  ID1_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
  num2 <- subset(num1, num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
  
  #補助情報のIDベクトルを作成
  ID2_list[[i]] <- rep(i, x[i])
  num1 <- (AX[i, ] > 0) * (1:a)
  num2 <- subset(num1, num1 > 0)
  A1 <- AX[i, (AX[i, ] > 0)]
  number <- rep(num2, A1)
  ad_list[[i]] <- number
}

#リストをベクトルに変換
ID1_d <- unlist(ID1_list)
ID2_d <- unlist(ID2_list)
wd <- unlist(wd_list)
ad <- unlist(ad_list)

##インデックスを作成
doc1_list <- list()
word_list <- list()
doc2_list <- list()
aux_list <- list()
for(i in 1:length(unique(ID1_d))) {doc1_list[[i]] <- subset(1:length(ID1_d), ID1_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- subset(1:length(wd), wd==i)}
for(i in 1:length(unique(ID2_d))) {doc2_list[[i]] <- subset(1:length(ID2_d), ID2_d==i)}
for(i in 1:length(unique(ad))) {aux_list[[i]] <- subset(1:length(ad), ad==i)}
gc(); gc()


####マルコフ連鎖モンテカルロ法で対応トピックモデルを推定####
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

##アルゴリズムの設定
R <- 10000   #サンプリング回数
keep <- 2   #2回に1回の割合でサンプリング結果を格納
iter <- 0

##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- rep(1.0, k)
beta0 <- rep(0.5, v)
gamma0 <- rep(0.5, a)
alpha01m <- matrix(alpha01, nrow=d, ncol=k, byrow=T)
beta0m <- matrix(beta0, nrow=v, ncol=k)
gamma0m <- matrix(gamma0, nrow=a, ncol=k)
delta0m <- gamma0

##パラメータの初期値
theta.ini <- runif(k, 0.5, 2)
phi.ini <- runif(v, 0.5, 1)
omega.ini <- runif(a, 0.5, 1)
theta <- rdirichlet(d, theta.ini)   #文書トピックのパラメータの初期値
phi <- rdirichlet(k, phi.ini)   #単語トピックのパラメータの初期値
omega <- rdirichlet(k, omega.ini)   #タグのトピックのパラメータの初期値
gamma <- rdirichlet(1, rep(1, a))   #内容と関係のタグのパラメータの初期値
r <- c(0.5, 0.5)   #内容に関係があるかどうかの混合率

##パラメータの格納用配列
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
OMEGA <- array(0, dim=c(k, a, R/keep))
GAMMA <- matrix(0, nrow=R/keep, ncol=a)
TAU_Z <- rep(0, length(ad))
W_SEG <- matrix(0, nrow=R/keep, ncol=k)
A_SEG <- matrix(0, nrow=R/keep, ncol=k)
storage.mode(W_SEG) <- "integer"
storage.mode(A_SEG) <- "integer"
storage.mode(TAU_Z) <- "integer"
gc(); gc()

##MCMC推定用配列
wsum0 <- matrix(0, nrow=d, ncol=k)
vf0 <- matrix(0, nrow=v, ncol=k)
af0 <- matrix(0, nrow=a, ncol=k)
aux_z <- rep(0, length(ad))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##単語トピックをサンプリング
  #単語ごとにトピックの出現率を計算
  word_rate <- burden_fr(theta, phi, wd, w, k)$Br
  
  #多項分布から単語トピックをサンプリング
  word_cumsums <- rowCumsums(word_rate)
  rand <- matrix(runif(nrow(word_rate)), nrow=nrow(word_rate), ncol=k)   #一様乱数
  word_z <- (k+1) - (word_cumsums > rand) %*% rep(1, k)   #トピックをサンプリング
  
  #発生させたトピックをダミー変数化
  Zi1 <- matrix(0, nrow(word_z), k)
  for(i in 1:k){
    index <- which(word_z==i)
    Zi1[index, i] <- 1
  }

  ##単語トピックのパラメータを更新
  #ディクレリ分布からthetaをサンプリング
  for(i in 1:d){wsum0[i, ] <- colSums(Zi1[doc1_list[[i]], ])}
  wsum <- wsum0 + alpha01m 
  theta <- t(apply(wsum, 1, function(x) rdirichlet(1, x)))
  
  #ディクレリ分布からphiをサンプリング
  for(i in 1:v){vf0[i, ] <- colSums(Zi1[word_list[[i]], ])}
  vf <- vf0 +  + beta0m
  phi <- t(apply(t(vf), 1, function(x) rdirichlet(1, x)))
  
  
  ##タグが文書のトピックと関連があるかどうかを決定
  #発生させた単語トピックからトピック抽出確率を計算
  aux_sums <- wsum - alpha01m
  theta_aux <- aux_sums/rowSums(aux_sums)
  
  #ベルヌーイ分布の確率を計算
  tau01 <- r[1]*rowSums(theta_aux[ID2_d, ]*t(omega[, ad]))
  tau02 <- r[2]*gamma[, ad]
  tau <- tau01 / (tau01 + tau02)
  
  #ベルヌーイ分布より潜在変数を発生
  z <- rbinom(length(tau), 1, tau)
  r <- c(mean(z), 1-mean(z))   #混合率
  
  
  ##補助情報トピックをサンプリング
  #z=1の場合、タグごとにトピックの出現率を計算
  index_z <- which(z==1)   #z=1のみ抽出
  aux_rate0 <- burden_fr(theta_aux, omega, ad, x, k)$Br
  aux_rate <- aux_rate0[index_z, ]  
  
  #多項分布から補助情報トピックをサンプリング
  aux_cumsums <- rowCumsums(aux_rate)
  rand <- matrix(runif(nrow(aux_rate)), nrow=nrow(aux_rate), ncol=k)   #一様乱数
  aux_z[index_z] <- (k+1) - (aux_cumsums > rand) %*% rep(1, k)   #トピックをサンプリング
  aux_z[-index_z] <- k+1
  
  #発生させたトピックをダミー変数化
  Zi2 <- matrix(0, length(aux_z), ncol=k+1)
  for(i in 1:(k+1)){
    index <- which(aux_z==i)
    Zi2[index, i] <- 1
  }
  
  ##タグトピックのパラメータを更新
  af <- as.matrix(data.frame(id=ad[index_z], Br=Zi2[index_z, -(k+1)]) %>%
                    dplyr::group_by(id) %>%
                    dplyr::summarize_all(funs(sum)))[, 2:(k+1)] + gamma0m
  omega <- t(apply(t(af), 1, function(x) rdirichlet(1, x)))
  
  ##内容に関係のないタグのパラメータの更新
  gamma0 <- delta0m
  gamma_sum <- tapply(Zi2[-index_z, k+1], ad[-index_z], sum)
  index_gamma <- as.numeric(names(gamma_sum))
  gamma0[index_gamma] <- gamma_sum   #ディクレリ分布のパラメータ
  gamma <- rdirichlet(1, gamma0)   #ディクレリ分布よりパラメータをサンプリング
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep1 <- rp/keep
    THETA[, , mkeep1] <- theta
    PHI[, , mkeep1] <- phi
    OMEGA[, , mkeep1] <- omega
    GAMMA[mkeep1, ] <- gamma
    
    #トピック割当はサンプリング期間の半分を超えたら格納する
    if(rp >= R/2){
      TAU_Z <- TAU_Z + z
      W_SEG <- W_SEG + Zi1
      A_SEG <- aux_z + Zi2
    }

    #サンプリング結果を確認
    print(rp)
    print(round(c(r, r0), 3))
    print(round(cbind(theta[1:5, ], thetat[1:5, ]), 3))
    print(round(rbind(gamma, gammat), 3))
    #print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
    #print(round(cbind(omega[, 1:10], omegat[, 1:10]), 3))
  }
}

####サンプリング結果の可視化と要約####
burnin <- 2000   #バーンイン期間

##サンプリング結果の可視化
#文書のトピック分布のサンプリング結果
matplot(t(THETA[1, , ]), type="l", ylab="パラメータ", main="文書1のトピック分布のサンプリング結果")
matplot(t(THETA[2, , ]), type="l", ylab="パラメータ", main="文書2のトピック分布のサンプリング結果")
matplot(t(THETA[3, , ]), type="l", ylab="パラメータ", main="文書3のトピック分布のサンプリング結果")
matplot(t(THETA[4, , ]), type="l", ylab="パラメータ", main="文書4のトピック分布のサンプリング結果")

#単語の出現確率のサンプリング結果
matplot(t(PHI[1, 1:10, ]), type="l", ylab="パラメータ", main="トピック1の単語の出現率のサンプリング結果")
matplot(t(PHI[2, 11:20, ]), type="l", ylab="パラメータ", main="トピック2の単語の出現率のサンプリング結果")
matplot(t(PHI[3, 21:30, ]), type="l", ylab="パラメータ", main="トピック3の単語の出現率のサンプリング結果")
matplot(t(PHI[4, 31:40, ]), type="l", ylab="パラメータ", main="トピック4の単語の出現率のサンプリング結果")
matplot(t(PHI[5, 41:50, ]), type="l", ylab="パラメータ", main="トピック5の単語の出現率のサンプリング結果")
matplot(t(PHI[6, 51:60, ]), type="l", ylab="パラメータ", main="トピック6の単語の出現率のサンプリング結果")
matplot(t(PHI[7, 61:70, ]), type="l", ylab="パラメータ", main="トピック7の単語の出現率のサンプリング結果")
matplot(t(PHI[8, 71:80, ]), type="l", ylab="パラメータ", main="トピック8の単語の出現率のサンプリング結果")

#タグの出現確率のサンプリング結果
matplot(t(OMEGA[1, 1:10, ]), type="l", ylab="パラメータ", main="トピック1のタグの出現率のパラメータのサンプリング結果")
matplot(t(OMEGA[2, 6:15, ]), type="l", ylab="パラメータ", main="トピック2のタグの出現率のパラメータのサンプリング結果")
matplot(t(OMEGA[3, 16:25, ]), type="l", ylab="パラメータ", main="トピック3のタグの出現率のパラメータのサンプリング結果")
matplot(t(OMEGA[4, 21:30, ]), type="l", ylab="パラメータ", main="トピック4のタグの出現率のパラメータのサンプリング結果")
matplot(t(OMEGA[5, 26:35, ]), type="l", ylab="パラメータ", main="トピック5のタグの出現率のパラメータのサンプリング結果")
matplot(t(OMEGA[6, 31:40, ]), type="l", ylab="パラメータ", main="トピック6のタグの出現率のパラメータのサンプリング結果")
matplot(t(OMEGA[7, 36:45, ]), type="l", ylab="パラメータ", main="トピック7のタグの出現率のパラメータのサンプリング結果")
matplot(t(OMEGA[8, 41:50, ]), type="l", ylab="パラメータ", main="トピック8のタグの出現率のパラメータのサンプリング結果")
matplot(GAMMA[, 41:50], type="l", ylab="パラメータ", main="トピックと無関係のタグの出現率のパラメータのサンプリング結果1")
matplot(GAMMA[, 51:60], type="l", ylab="パラメータ", main="トピックと無関係のタグの出現率のパラメータのサンプリング結果2")

##サンプリング結果の要約推定量
#トピック分布の事後推定量
topic_mu <- apply(THETA[, , burnin:(R/keep)], c(1, 2), mean)   #トピック分布の事後平均
round(cbind(topic_mu, thetat), 3)
round(topic_sd <- apply(THETA[, , burnin:(R/keep)], c(1, 2), sd), 3)   #トピック分布の事後標準偏差

#単語出現確率の事後推定量
word_mu <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #単語の出現率の事後平均
round(rbind(word_mu, phit)[, 1:50], 3)

#タグ出現率の事後推定量
tag_mu1 <- apply(OMEGA[, , burnin:(R/keep)], c(1, 2), mean)   #タグの出現率の事後平均
round(rbind(tag_mu1, omegat), 3)

#トピックと無関係のタグの事後推定量
round(rbind(colMeans(GAMMA[burnin:(R/keep), ]), gammat), 3)   #無関係タグの事後平均

