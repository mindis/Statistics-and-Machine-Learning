#####split hazard model#####
library(MASS)
library(survival)
library(reshape2)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

####データの発生####
##データの設定
zk <- 3   #潜在変数数
n <- 2000   #サンプル数
tmax <- 10   #最大観測数
nmax <- n*tmax   #最大サンプル数

##潜在変数の発生
#デモグラフィック変数の発生
#連続変数の発生
cont_z <- matrix(rnorm(n*2), nrow=n, ncol=2)

#二値変数の発生
bin_z <- matrix(0, nrow=n, ncol=3)
pbin <- c(0.4, 0.3, 0.6)
for(i in 1:length(pbin)){
  bin_z[, i] <- rbinom(n, 1, pbin[i])
}

#多値変数の発生
pmulti <- c(0.2, 0.1, 0.3, 0.3, 0.1)
multi_z <- t(rmultinom(n, 1, pmulti))
multi_z <- multi_z[, -5]

#パラメータの設定
bz1 <- matrix(runif((zk-1)*ncol(cont_z), -1.0, 1.0), nrow=2, ncol=ncol(cont_z))
bz2 <- matrix(runif((zk-1)*ncol(bin_z), -1.2, 1.2), nrow=2, ncol=ncol(bin_z))
bz3 <- matrix(runif((zk-1)*ncol(multi_z), -1.5, 1.5), nrow=2, ncol=ncol(multi_z))
bz0 <- c(-0.4, 0.7)
bz_t <- cbind(bz1, bz2, bz3)

#効用関数の定義
U1 <- bz0[1] + cbind(cont_z, bin_z, multi_z) %*% bz_t[1, ] 
U2 <- bz0[1] + cbind(cont_z, bin_z, multi_z) %*% bz_t[2, ]

#確率の計算と潜在変数zの発生
d <- cbind(0, U1, U2)
Pr_z <- exp(d) / rowSums(exp(d))
Z <- t(apply(Pr_z, 1, function(x) rmultinom(1, 1, x)))

#潜在変数の要約
colSums(Z)
round(colMeans(Z), 3)

#潜在変数の割当
z.index <- rep(Z %*% 1:zk, rep(tmax, n))
id <- rep(1:n, rep(tmax, n))
time <- rep(1:tmax, n)

##ハザードモデルのパラメータの設定


##ハザードモデルの説明変数の発生
page_cnt <- 10

##ページ閲覧回数と閲覧履歴の発生
#ページ閲覧回数の発生
lam_lower <- 5
lam_upper <- 9
p_cnt.zero <- rpois(nmax, runif(nmax, lam_lower, lam_upper))
p_cnt <- ifelse(p_cnt.zero==0, 1, p_cnt.zero)
hist(p_cnt, breaks=15, col="grey", xlab="ページ閲覧数", main="ページ閲覧数の分布")

#ページ閲覧履歴の発生
p_rate <- runif(page_cnt)
p_hist <- matrix(0, nmax, page_cnt)

for(i in 1:nmax){
  p_hist[i, ] <- t(rmultinom(1, p_cnt[i], p_rate))
}
p_hist
p_hist.r <- p_hist / rowSums(p_hist)

#離脱時のページの発生
p_last <- t(rmultinom(N, 1, p_rate))

#前回のページ閲覧でもっとも見ていたページ
index.m <- subset(1:nrow(p_hist), apply(p_hist, 1, max) > 1)

pm <- t(apply((p_hist-2), 1, function(x) x-abs(max(x))))
p_most1 <- ifelse(pm==0, 1, 0)

#1番目のアクセスはすべて0になる
p_most2 <- rbind(0, p_most1)
p_most2[time==1, ] <- rep(0, page_cnt)
p_most <- p_most2[1:nrow(p_hist), ]


##前回からのアクセス経過時間(単位=日)の発生
shape <- 1.65
rate <- 0.6
t <- round(rgamma(nmax, shape, rate), 0)
index.t <- subset(1:length(t), t == 0)
t[index.t] <- round(runif(length(index.t), 1, 5), 0)


##冗長な変数を削除してデータを結合
index.h <- which.min(colSums(p_hist))
ph_hist <- p_hist[, -index.h]
ph_hist.r <- p_hist.r[, -index.h]
ph_last <- p_last[, -index.h]

X <- data.frame(page=ph_hist.r, last_p=ph_last, cnt_p=p_cnt, t=t)   #データの結合
X
round(X, 3)
