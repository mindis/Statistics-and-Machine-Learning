#####ノイズあり対応トピックモデル#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(8079)

####データの生成####
#set.seed(423943)
#文書データの設定
k <- 15   #トピック数
d <- 2500   #文書数
v <- 1000   #語彙数
w <- rpois(d, rgamma(d, 75, 0.5))   #1文書あたりの単語数
f <- sum(w)   #総単語数
a1 <- 100   #トピックに関係のあるタグ数
a2 <- 50   #トピックに関係のないタグ数
a <- a1 + a2   #補助変数数
x <- rtpois(d, 30, 2, Inf)
f1 <- sum(w)
f2 <- sum(x)

#IDの設定
w_id <- rep(1:d, w)
a_id <- rep(1:d, x)

#パラメータの設定
alpha0 <- rep(0.15, k)   #文書のディレクリ事前分布のパラメータ
alpha1 <- rep(0.15, v)   #単語のディレクリ事前分布のパラメータ
alpha2 <- c(rep(0.2, a1), rep(0.001, a2))   #トピックに関係のあるタグのディクレリ事前分布のパラメータ
alpha3 <- c(rep(0.001, a1), rep(1.0, a2))   #トピックに関係のないタグのディクレリ事前分布のパラメータ
beta0 <- rbeta(sum(x), 0.55, 0.175)

##モデルに基づき単語を生成
for(rp in 1:1000){
  print(rp)
  
  #ディレクリ分布からパラメータを生成
  thetat <- theta <- extraDistr::rdirichlet(d, alpha0)   #文書のトピック分布をディレクリ乱数から生成
  phit <- phi <- extraDistr::rdirichlet(k, alpha1)   #単語のトピック分布をディレクリ乱数から生成
  lambda <- matrix(0, nrow=d, ncol=k)   #文書に含むトピックだけを補助情報のトピックにするための確率を格納する行列
  omegat <- omega <- extraDistr::rdirichlet(k, alpha2)   #補助情報のトピック分布をディクレリ乱数から生成
  gammat <- gamma <- extraDistr::rdirichlet(1, alpha3)   #トピックに関係のないタグ
  omega0 <- rbind(omega, gamma)   #単語分布の結合
  
  ##多項分布からトピックおよび単語データを生成
  WX <- matrix(0, nrow=d, ncol=v)
  AX <- matrix(0, nrow=d, ncol=a)
  word_list <- list()
  aux_list <- list()
  Z0 <- rep(0, sum(x)) 
  Z1_list <- list()
  Z2_list <- list()
  
  #文書ごとにトピックと単語を逐次生成
  for(i in 1:d){
    
    #文書のトピック分布を生成
    z1 <- rmnom(w[i], 1, theta[i, ])   #文書のトピック分布を生成
    z1_vec <- as.numeric(z1 %*% 1:k)
    
    #文書のトピック分布から単語を生成
    word <- rmnom(w[i], 1, phi[z1_vec, ])   #文書のトピックから単語を生成
    word_vec <- colSums(word)   #単語ごとに合計して1行にまとめる
    WX[i, ] <- word_vec
    
    #文書のトピック分布から補助変数を生成
    #文書で生成させたトピックのみを補助情報のトピック分布とする
    lambda[i, ] <- colSums(z1) / w[i]   #補助情報のトピック分布
    
    #ベルヌーイ分布からトピックに関係があるかどうかを生成
    index <- which(a_id==i)
    Z0[index] <- rbinom(length(index), 1, beta0[index])
  
    #補助情報のトピックを生成
    z2_aux <- rmnom(x[i], 1, lambda[i, ])
    z2 <- cbind(z2_aux * Z0[index], 1-Z0[index])
    z2_vec <- as.numeric(z2 %*% 1:(k+1))
    
    #生成させたトピックの単語分布に従い単語を生成
    aux <- rmnom(x[i], 1, omega0[z2_vec, ])
    aux_vec <- colSums(aux)
    AX[i, ] <- aux_vec
    
    #文書トピックおよび補助情報トピックを格納
    Z1_list[[i]] <- z1
    Z2_list[[i]] <- z2
    word_list[[i]] <- as.numeric(word %*% 1:v)
    aux_list[[i]] <- as.numeric(aux %*% 1:a)
  }
  if(min(colSums(AX)) > 0 & min(colSums(WX)) > 0){
    break
  }
}

#データ行列を整数型行列に変更
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)
wd <- unlist(word_list)
ad <- unlist(aux_list)
storage.mode(WX) <- "integer"
storage.mode(AX) <- "integer"
r0 <- c(mean(Z0), 1-mean(Z0))


####トピックモデル推定のためのデータと関数の準備####
##インデックスを作成
doc1_list <- list()
doc1_vec <- list()
word_list <- list()
word_vec <- list()
doc2_list <- list()
doc2_vec <- list()
aux_list <- list()
aux_vec <- list()

for(i in 1:d){
  doc1_list[[i]] <- which(w_id==i)
  doc1_vec[[i]] <- rep(1, length(doc1_list[[i]]))
  doc2_list[[i]] <- which(a_id==i)
  doc2_vec[[i]] <- rep(1, length(doc2_list[[i]]))
}
for(j in 1:v){
  word_list[[j]] <- which(wd==j)
  word_vec[[j]] <- rep(1, length(word_list[[j]]))
}
for(j in 1:a){
  aux_list[[j]] <- which(ad==j)
  aux_vec[[j]] <- rep(1, length(aux_list[[j]]))
}
gc(); gc()

####マルコフ連鎖モンテカルロ法で対応トピックモデルを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k){
  #負担係数を計算
  Bur <- theta[w, ] * t(phi)[wd, ]   #尤度
  Br <- Bur / rowSums(Bur)   #負担率
  r <- colSums(Br) / sum(Br)   #混合率
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}


##アルゴリズムの設定
R <- 5000   #サンプリング回数
keep <- 2   #2回に1回の割合でサンプリング結果を格納
disp <- 10
iter <- 0
burnin <- 1000/keep

##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- 1
alpha02 <- 1
beta01 <- 0.1
beta02 <- 0.1
gamma01 <- 1
gamma02 <- 1

##パラメータの初期値
#tfidfで初期値を設定
tf <- AX / rowSums(AX)
idf1 <- log(nrow(AX) / colSums(AX > 0))
idf2 <- log(nrow(AX) / colSums(AX==0))

theta <- extraDistr::rdirichlet(d, rep(1, k))   #文書トピックのパラメータの初期値
phi <- extraDistr::rdirichlet(k, colSums(WX)/sum(WX)*100)   #単語トピックのパラメータの初期値
omega <- extraDistr::rdirichlet(k, idf1)   #タグのトピックのパラメータの初期値
gamma <- extraDistr::rdirichlet(1, idf2)   #内容と関係のタグのパラメータの初期値
r <- c(0.5, 0.5)   #内容に関係があるかどうかの混合率

##パラメータの格納用配列
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
OMEGA <- array(0, dim=c(k, a, R/keep))
GAMMA <- matrix(0, nrow=R/keep, ncol=a)
Z_SEG <- rep(0, f2)
W_SEG <- matrix(0, nrow=f1, ncol=k)
A_SEG <- matrix(0, nrow=f2, ncol=k+1)
storage.mode(W_SEG) <- "integer"
storage.mode(A_SEG) <- "integer"
storage.mode(Z_SEG) <- "integer"
gc(); gc()


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##単語トピックをサンプリング
  #単語ごとにトピックの出現率を計算
  word_rate <- burden_fr(theta, phi, wd, w, k)$Br
  
  #多項分布から単語トピックをサンプリング
  Zi1 <- rmnom(f, 1, word_rate)   
  z_vec1 <- Zi1 %*% 1:k
  
  ##単語トピックのパラメータを更新
  #ディクレリ分布からthetaをサンプリング
  for(i in 1:d){
    wsum0[i, ] <- colSums(Zi1[doc1_list[[i]], ])
  }
  wsum <- wsum0 + alpha01m 
  theta <- extraDistr::rdirichlet(d, wsum)
  
  #ディクレリ分布からphiをサンプリング
  for(j in 1:v){
    vf0[, j] <- colSums(Zi1[word_list[[j]], ])
  }
  vf <- vf0 + beta0m
  phi <- extraDistr::rdirichlet(k, vf)
  
  
  ##タグが文書のトピックと関連があるかどうかを決定
  #生成させた単語トピックからトピック抽出確率を計算
  aux_sums <- wsum - alpha01m
  theta_aux <- aux_sums/rowSums(aux_sums)
  
  #ベルヌーイ分布の確率を計算
  tau01 <- r[1] * rowSums(theta_aux[ID2_d, ]*t(omega[, ad]))
  tau02 <- r[2] * gamma[, ad]
  tau <- tau01 / (tau01 + tau02)
  
  #ベルヌーイ分布より潜在変数を生成
  z <- rbinom(e, 1, tau)
  r <- c(mean(z), 1-mean(z))   #混合率
  
  ##補助情報トピックをサンプリング
  #z=1の場合、タグごとにトピックの出現率を計算
  index_z <- which(z==1)   #z=1のみ抽出
  aux_rate <- burden_fr(theta_aux, omega, ad, x, k)$Br
  
  #多項分布から補助情報トピックをサンプリング
  Zi2 <- cbind(rmnom(e, 1, aux_rate), 0)   
  Zi2[-index_z, k+1] <- 1
  z_vec1 <- Zi2 %*% 1:(k+1)
  
  ##タグトピックのパラメータを更新
  af <- as.matrix(data.frame(id=ad[index_z], Br=Zi2[index_z, -(k+1)]) %>%
                    dplyr::group_by(id) %>%
                    dplyr::summarize_all(funs(sum)))[, 2:(k+1)] + gamma0m
  omega <- extraDistr::rdirichlet(k, t(af))
  
  ##内容に関係のないタグのパラメータの更新
  gamma0 <- delta0m
  par <- tapply(c(Zi2[-index_z, k+1], rep(1, a)), c(ad[-index_z], 1:a), sum) - 1 + gamma0
  gamma <- rdirichlet(1, par)   #ディクレリ分布よりパラメータをサンプリング
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep1 <- rp/keep
    THETA[, , mkeep1] <- theta
    PHI[, , mkeep1] <- phi
    OMEGA[, , mkeep1] <- omega
    GAMMA[mkeep1, ] <- gamma
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(rp >= burnin){
      TAU_Z <- TAU_Z + z
      W_SEG <- W_SEG + Zi1
      A_SEG <- A_SEG + Zi2
    }
    
    #サンプリング結果を確認
    print(rp)
    print(round(c(r, r0), 3))
    print(round(cbind(theta[1:10, ], thetat[1:10, ]), 3))
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
round(rbind(colMeans(GAMMA[burnin:(R/keep), ]), gammat), 3) #無関係タグの事後平均