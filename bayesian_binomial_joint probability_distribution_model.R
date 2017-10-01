####二項分布モデルのパラメータの同時推定モデル####

####データの発生####
n <- 800
p <- 0.6
x <- rbinom(1, n, p)

####マルコフ連鎖モンテカルロ法の設定####
##アルゴリズムの設定
R <- 10000
keep <- 2

##事前分布の設定
a1 <- 1
b1 <- 1
num <- (x+100):1000
p1 <- rep(1/length((x+100):1000), length((x+100):1000))

##初期値の設定
N <- x * 2
Pr <- x / N

##サンプリング結果の保存用配列
Sample <- rep(0, R/keep)
Prob <- rep(0, R/keep)

####MCMCでパラメータをサンプリング#####
for(rp in 1:R){

  ##サンプル数Nをサンプリング
  old_N <- N
  new_N <- sum((t(rmultinom(1, 1, p1)) * num))
  
  lognew <- dbinom(x, new_N, Pr, log=T)
  logold <- dbinom(x, old_N, Pr, log=T)
  
  #MHサンプリング
  alpha <- min(1, exp(lognew - logold))
  if(alpha == "NAN") alpha <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しい固定効果betaを採択
  if(u < alpha){
    N <- new_N
    
    #そうでないなら固定効果betaを更新しない
  } else {
    N <- old_N
  }
  
  ##ベータ分布から回収率をサンプリング
  Pr <- rbeta(1, a1+x, b1+N-x)

  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep   
    Sample[mkeep] <- N
    Prob[mkeep] <- Pr
    print(rp)
  }
}


####サンプリング結果の確認と可視化
burnin <- 1000

#MCMCサンプリング結果の可視化
matplot(1:(R/keep), Sample, type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(1:(R/keep), Prob, type="l", xlab="サンプリング回数", ylab="パラメータ")
plot(Sample[burnin:(R/keep)], Prob[burnin:(R/keep)], xlab="サンプル数", ylab="回収率")

#パラメータの事後平均
mean(Prob[burnin:(R/keep)])
mean(Sample[burnin:(R/keep)])


