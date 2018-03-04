#####Collaborative Topic Model#####
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

#set.seed(93441)

####データの発生####
##データの設定
k <- 10
hh <- 1500   #レビュー人数
item <- 500   #アイテム数
s <- rpois(hh, rgamma(hh, 7.5, 0.6))   #1人あたりのレビュー数
s[s==0] <- ceiling(runif(length(s[s==0]), 1, 5))
d <- sum(s)   #総文書数
v <- 1000   #語彙数
w <- rpois(d, rgamma(d, 40, 0.8))   #文書あたりの単語数
f <- sum(w)   #総単語数

##アイテム購買を設定
u_list <- list()
par <- as.numeric(extraDistr::rdirichlet(1, rep(3.0, item)))   #アイテム購買確率
for(rp in 1:1000){
  for(i in 1:hh){
    for(j in 1:1000){
      pi <- rmnom(s[i], 1, par)
      if(max(colSums(pi))==1){
        break
      }
    }
    u_list[[i]] <- pi
  }
  if(min(colSums(do.call(rbind, u_list)))==0){
    break
  }
}
U <- do.call(rbind, u_list)
u_vec <- as.numeric(U %*% 1:item)
colSums(U)

##IDを設定
u_id <- rep(1:hh, s)
w_id <- rep(u_vec, w)


##パラメータの設定
#トピックモデルのパラメータ
alpha01 <- rep(0.2, k)
alpha02 <- rep(0.15, v)
theta1 <- thetat1 <- extraDistr::rdirichlet(item, alpha01)
phi <- phit <- extraDistr::rdirichlet(k, alpha02)

#潜在因子モデルのパラメータ
tau1 <- 0.5
tau2 <- 0.025
sigma <- 0.2
theta2 <- thetat2 <- mvrnorm(hh, rep(0, k), diag(tau1, k))
psi <- psit <- t(theta1 + mvrnorm(item, rep(0, k), diag(tau2, k)))


##Collavorative topic modelのデータを生成
Z_list <- list()
WX <- matrix(0, nrow=d, ncol=v)
word_list <- list()
score <- rep(0, d)

for(i in 1:d){
  #アイテムトピックと単語を生成
  z <- rmnom(w[i], 1, theta1[u_vec[i], ])   #トピックを生成
  z_vec <- as.numeric(z %*% 1:k)
  words <- rmnom(w[i], 1, phi[z_vec, ])   #単語を生成
  words_vec <- colSums(words)
  
  #スコアを生成
  r <- as.numeric(theta2[u_id[i], ] %*% psi[, u_vec[i]]) + rnorm(1, 0, sigma)
  
  #パラメータを格納
  Z_list[[i]] <- z
  WX[i, ] <- words_vec
  word_list[[i]] <- as.numeric(words %*% 1:v)
  score[i] <- r
}

#リストを変換
Z <- do.call(rbind, Z_list)
words_vec <- unlist(word_list)
storage.mode(WX) <- "integer"


####EMアルゴリズムでCollaborative Topic Modelを推定####
##単語ごとに尤度と負担率を計算する関数
burden_fr <- function(theta, phi, wd, w, k){
  #負担係数を計算
  Bur <- theta1[w_id, ] * t(phi)[words_vec, ]   #尤度
  Br <- Bur / rowSums(Bur)   #負担率
  r <- colSums(Br) / sum(Br)   #混合率
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##インデックスの作成
#トピックモデルのインデックスを作成
doc_list <- list()
word_list <- list()
doc_ones <- list()
word_ones <- list()
doc_n <- rep(0, item)
word_n <- rep(0, v)

for(i in 1:item){
  doc_list[[i]] <- which(w_id==i)
  doc_n[i] <- length(doc_list[[i]])
  doc_ones[[i]] <- rep(1, doc_n[i])
}
for(i in 1:v){
  word_list[[i]] <- which(words_vec==i)
  word_n[i] <- length(word_list[[i]])
  word_ones[[i]] <- rep(1, word_n[i])
}

#潜在因子モデルのインデックスを作成
user_list <- list()
user_n <- rep(0, hh)
item_list <- list()
item_n <- rep(0, item)
for(i in 1:hh){
  user_list[[i]] <- which(u_id==i)
  user_n[i] <- length(user_list[[i]])
}
for(i in 1:item){
  item_list[[i]] <- which(u_vec==i)
  item_n[i] <- length(item_list[[i]])
}


#初期値の設定
theta1 <- extraDistr::rdirichlet(item, rep(0.5, k))
phi <- extraDistr::rdirichlet(k, rep(1.0, v))
theta2 <- extraDistr::rdirichlet(hh, rep(0.5, k))
sigma <- 0.2   #誤差標準偏差の初期値
tau1 <- 0.025   #アイテム因子のハイパーパラメータ
tau2 <- 0.25   #ユーザー因子のハイパーパラメータ


#EMアルゴリズムの更新ステータス
iter <- 1
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 0.1
LLo <- -1000000000   #対数尤度の初期値
LLw <- c()

##EMアルゴリズムでパラメータを更新
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  
  ##トピックモデルで単語のトピックを推定
  ##尤度と潜在変数zを推定
  par <- burden_fr(theta1, phi, words_vec, w_id, k)
  z <- par$Br   #潜在変数zを出力
  r <- par$r   #混合率rを出力
  
  ##トピックモデルのパラメータを更新
  #トピック分布theta1を更新
  wsum <- matrix(0, nrow=item, ncol=k)
  for(i in 1:item){
    wsum[i, ] <- t(z[doc_list[[i]], ]) %*% doc_ones[[i]]
  }
  theta1 <- wsum / doc_n
  theta1[is.nan(theta1)] <- rep(1/k, k)
  
  #単語分布phiを更新
  vsum <- matrix(0, nrow=k, ncol=v)
  for(j in 1:v){
    vsum[, j] <- t(z[word_list[[j]], ]) %*% word_ones[[j]]
  }
  phi <- vsum / matrix(rowSums(vsum), nrow=k, ncol=v)
  
  
  ##潜在因子モデルでスコア要因を推定
  ##アイテム因子psiを推定
  #データの設定
  x0 <- thetat2[u_id, ]
  lambda <- diag(tau1, k)
  
  #アイテムごとに因子行列を最適化
  for(j in 1:item){
    
    #アイテムがない場合は次のアイテムへ
    if(item_n[j]==0){
      psi[, j] <- theta1[j, ]
      next
    }
    #アイテムごとにデータを設定
    index <- item_list[[j]]
    x <- x0[index, ]
    y <- score[index]
    omega <- diag(sigma^2, item_n[j])
    
    #psiを正則化最小二乗法で最適化
    if(item_n[j] > 1){
      psi[, j] <- as.numeric(solve(t(x) %*% omega %*% x + lambda) %*% (t(x) %*% omega %*% y + lambda %*% theta1[j, ]))
    } else {
      psi[, j] <- as.numeric(solve(x %*% omega %*% x + lambda) %*% (x %*% omega %*% y + lambda %*% theta1[j, ]))
    }
  }
  psi[is.nan(psi)] <- 0
  
  ##ユーザー因子theta2を推定
  #データの設定
  x0 <- t(psi)[u_vec, ]
  lambda <- diag(tau2, k)
  
  #ユーザーごとに因子行列を最適化
  for(i in 1:hh){
    #ユーザーごとにデータを設定
    index <- user_list[[i]]
    x <- x0[index, ]
    y <- score[index]
    omega <- diag(sigma, user_n[i])

    #theta2を正則化最小二乗法で最適化
    if(user_n[i] > 1){
      theta2[i, ] <- as.numeric(solve(t(x) %*% omega %*% x + lambda) %*% t(x) %*% omega %*% y)
    } else {
      theta2[i, ] <- as.numeric(solve(x %*% omega %*% x + lambda) %*% x %*% omega %*% y)
    }
  }
  
  ##誤差sigmaを更新
  mu <- rowSums(theta2[u_id, ] * t(psi)[u_vec, ])
  er <- score - mu
  sigma <- sum(er^2) / (d-1)

  
  ##EMアルゴリズムの更新
  #対数尤度の計算
  LLS1 <- sum(log(rowSums(par$Bur)))   #トピックモデルの対数尤度
  LLS2 <- -sum(dnorm(score, mu, sqrt(sigma), log=TRUE))   #潜在因子モデルの対数尤度
  LLS <- LLS1 + LLS2
  
  #収束判定
  iter <- iter+1
  dl <- LLS-LLo
  LLo <- LLS
  LLw <- c(LLw, LLo)
  print(LLo)
}

#####推定結果の要約#####
##パラメータ推定値
round(cbind(par$Br, Z %*% 1:k), 3)   #生成トピック
round(cbind(theta1, thetat1), 2)   #トピック分布
round(cbind(t(phi), t(phit)), 3)   #単語分布
round(data.frame(n=user_n, x=theta2, y=thetat2), 2)   #ユーザー因子
round(data.frame(n=item_n, x=t(psi), y=t(psit)), 2)   #アイテム因子

##スコアの予測値
pred_mu0 <- theta2 %*% psi
pred_mu1 <- thetat2 %*% psit
round(cbind(pred_mu0[, 1:10], pred_mu1[, 1:10]), 2)



