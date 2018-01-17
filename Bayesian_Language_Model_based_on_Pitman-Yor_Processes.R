#####階層pitman-yor過程言語モデル#####
options(warn=2)
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

#set.seed(5723)
#データの設定
d <- 1000   #文書数
v <- 700   #語彙数
w <- rpois(d, rgamma(d, 80, 0.40))   #1文書あたりの単語数
f <- sum(w)

#IDの設定
u_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}
doc_list <- list()
for(i in 1:d){doc_list[[i]] <- which(u_id==i)}


##bi-gramモデルのパラメータを設定
alpha0 <- sort(rbeta(v, 0.65, 7.0), decreasing=TRUE)
theta01 <- thetat01 <- extraDistr::rdirichlet(1, alpha0)
theta02 <- thetat02 <- extraDistr::rdirichlet(v, alpha0) 
summary(as.numeric(theta02))

##単語を生成
wd0 <- rep(0, f)
WX0 <- matrix(0, nrow=f, ncol=v)
for(i in 1:d){
  if(i%%100==0){
    print(i)
  }
  index <- doc_list[[i]]
  freq <- length(index)
  x <- rep(0, freq)
  
  for(j in 1:freq){
    if(j==1){
      z <- rmnom(1, 1, theta01)
      x[j] <- as.numeric(z %*% 1:v)
      WX0[index[j], ] <- z
      
    } else {
      z <- rmnom(1, 1, theta02[x[j-1], ])
      x[j] <- as.numeric(z %*% 1:v)
      WX0[index[j], ] <- z
    }
  }
  wd0[index] <- x
}
word_freq <- as.numeric((table(c(wd0, 1:v))-1))   #単語の出現数
word_probs <- word_freq / f   #単語の出現確率
hist(word_probs, breaks=25, col="grey", xlab="出現確率", main="単語の出現確率の分布")

#頻度がゼロの単語は除く
index_zeros <- which(colSums(WX0)==0)
WX <- WX0[, -index_zeros]
v <- ncol(WX)
wd <- as.numeric(WX %*% 1:v)
theta1 <- thetat1 <- theta01[, -index_zeros] / sum(theta01[, -index_zeros])
theta2 <- thetat2 <- theta02[-index_zeros, -index_zeros] / rowSums(theta02[-index_zeros, -index_zeros])

##Bi-gramのデータ設計
Data <- matrix(0, nrow=v+1, ncol=v) 
index_ngram <- matrix(0, nrow=f, ncol=2)
for(i in 1:d){
  index <- doc_list[[i]]
  freq <- length(index)
  for(j in 1:freq){
    if(j==1){
      Data[v+1, wd[index[j]]] <- Data[v+1, wd[index[j]]] + 1
      index_ngram[index[j], ] <- c(v+1, wd[index[j]])
    } else {
      Data[wd[index[j-1]], wd[index[j]]] <- Data[wd[index[j-1]], wd[index[j]]] + 1
      index_ngram[index[j], ] <- c(wd[index[j-1]], wd[index[j]])
    }
  }
}

Data[index_ngram[1, 1], index_ngram[1, 2]]


