library(MASS)

##相関の変化による選択確率の変動をシミュレーションで確認
n <- 500000
rank <- c(0.5, 0.5, 0.5)
cc <- seq(-1, 1, length=101)
Rank.ones <- matrix(0, nrow=n, ncol=length(cc))
Rank.agg <- matrix(0, nrow=length(cc), ncol=length(rank))

for(i in 1:length(cc)){
  print(i)
  CORM <- matrix(c(1, cc[i], 0, cc[i], 1, 0, 0, 0, 1), nrow=3, ncol=3)
  U <- mvrnorm(n, rank, CORM)
  colnames(U) <- c("にこ", "真姫", "ことり")
  
  first <- apply(U, 1, which.max)
  Rank.ones[, i]　<- ifelse(first==1, "にこ", ifelse(first==2, "真姫", "ことり"))
  Rank.agg[i, ] <- table(Rank.ones[, i])
}

#変数に名前をつける
colnames(Rank.agg) <- c("ことり", "にこ", "真姫")
rownames(Rank.agg) <- round(cc, 2)

#可視化と結果の確認
matplot(cc, Rank.agg/n, type="l", xlab="相関係数", ylab="選択確率", lwd=2, main="メンバー間の相関による選択確率の変化")
round(Rank.agg/n, 3)


