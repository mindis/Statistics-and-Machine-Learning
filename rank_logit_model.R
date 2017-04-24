#####ランクロジットモデル#####
library(MASS)
library(mlogit)
library(nnet)
library(VGAM)
library(reshape2)
library(plyr)

####データの発生####
#set.seed(43890)
##パラメータの設定
betacs <- runif(9, 0, 2.0)   #衣装の回帰係数
betasv <- 4.2   #share of voiceの回帰係数
betace <- 0.25   #メンバーごとのセンター設定回数
betatm <- runif(9, -0.25, 0.25)   #登録からの経過時間の回帰係数
betapt <- runif(9, -0.25, 0.25)   #月あたりの平均プレイ時間の回帰係数
betaho <- 1.7   #穂乃果ちゃんの回帰係数
betako <- 2.5   #ことりちゃんの回帰係数
betaum <- 1.6   #海未ちゃんの回帰係数
betama <- 2.3   #真姫ちゃんの回帰係数
betaha <- 1.5   #かよちんの回帰係数
betari <- 1.2   #凛ちゃんの回帰係数
betano <- 1.1   #のんたんの回帰係数
betaer <- 1.9   #えりちかの回帰係数
betani <- 2.0   #にこにーの回帰係数
betamus <- c(betaho, betako, betaum, betama, betaha, betari, betano, betaer, betani)   #メンバーの回帰係数ベクトル

##データの設定
hh <- 3000   #プレイヤー数
pt <- 3   #3位まで選択
hhpt <- hh*pt   #データの行数
cl <- 9   #衣装数
member <- 10   #選択可能数

ID <- rep(1:hh, rep(pt, hh))   #プレイヤーID
RANK <- rep(1:pt, hh)   #ランク
U <- matrix(0, hhpt, 10)   #効用関数を格納
CLOTH.ho <- matrix(0, hhpt, 9)   #穂乃果ちゃんの衣装を格納
CLOTH.ko <- matrix(0, hhpt, 9)   #ことりちゃんの衣装を格納
CLOTH.um <- matrix(0, hhpt, 9)   #海未ちゃんの衣装を格納
CLOTH.ma <- matrix(0, hhpt, 9)   #真姫ちゃんの衣装を格納
CLOTH.ha <- matrix(0, hhpt, 9)   #かよちんの衣装を格納
CLOTH.ri <- matrix(0, hhpt, 9)   #凛ちゃんの衣装を格納
CLOTH.no <- matrix(0, hhpt, 9)   #のんたんの衣装を格納
CLOTH.er <- matrix(0, hhpt, 9)   #えりちかの衣装を格納
CLOTH.ni <- matrix(0, hhpt, 9)   #にこにーの衣装を格納
SoV <- matrix(0, hhpt, 10)   #メンバーごとのshare of voiceを格納
CENTER <- matrix(0, hhpt, 10)   #メンバーごとのセンター設定回数を格納
REGIST <- matrix(0, hhpt, 1)   #登録からの経過月数の対数を格納
PLAY <- matrix(0, hhpt, 1)   #月あたりのプレイ時間を格納


##データを発生させる
for(i in 1:hh){
  r <- (i-1)*pt+1
  ##メンバーの衣装を決定
  #穂乃果ちゃんの衣装
  p <- c(runif(cl))
  hono <- t(rmultinom(1, 1, p))
  hono.b <- rbind(hono, hono, hono)
  CLOTH.ho[r:(r+2), ] <- hono.b
  
  #ことりちゃんの衣装
  koto <- t(rmultinom(1, 1, p))
  koto.b <- rbind(koto, koto, koto)
  CLOTH.ko[r:(r+2), ] <- koto.b
  
  #海未ちゃんの衣装
  umi <- t(rmultinom(1, 1, p))
  umi.b <- rbind(umi, umi, umi)
  CLOTH.um[r:(r+2), ] <- umi.b
  
  #真姫ちゃんの衣装
  maki <- t(rmultinom(1, 1, p))
  maki.b <- rbind(maki, maki, maki)
  CLOTH.ma[r:(r+2), ] <- maki.b
  
  #かよちんの衣装
  kayo <- t(rmultinom(1, 1, p))
  kayo.b <- rbind(kayo, kayo, kayo)
  CLOTH.ha[r:(r+2), ] <- kayo.b
  
  #凛ちゃんの衣装
  rin <- t(rmultinom(1, 1, p))
  rin.b <- rbind(rin, rin, rin)
  CLOTH.ri[r:(r+2), ] <- rin.b
  
  #のんたんの衣装
  nozo <- t(rmultinom(1, 1, p))
  nozo.b <- rbind(nozo, nozo, nozo)
  CLOTH.no[r:(r+2), ] <- nozo.b
  
  #えりちの衣装
  eri <- t(rmultinom(1, 1, p))
  eri.b <- rbind(eri, eri, eri)
  CLOTH.er[r:(r+2), ] <- eri.b
  
  #にこにーの衣装
  nico <- t(rmultinom(1, 1, p))
  nico.b <- rbind(nico, nico, nico)
  CLOTH.ni[r:(r+2), ] <- nico.b
  
  ##メンバーごとのshare of voiceを決定
  sv.m1 <- c(betamus, runif(1, 0, 1)) + c(rnorm(9, 0, runif(1, 0, 3.5)), runif(1, 0, 2))
  sv.m2 <- ifelse(sv.m1 < 0, runif(1, 0, 1), sv.m1) 
  sv.m3 <- sv.m2 / sum(sv.m2)
  svr <- rbind(sv.m3, sv.m3, sv.m3)
  SoV[r:(r+2), ] <- svr
  
  ##メンバーごとのセンター設定回数を決定
  lambda <- c(betamus, runif(1, 0, 1.0)) + c(runif(9, 0, 6), runif(1, 0, 2))
  c_cnt <- rpois(10, lambda)
  c_cntr <- rbind(c_cnt, c_cnt, c_cnt)
  CENTER[r:(r+2), ] <- c_cntr
  
  ##登録からの経過月数の対数
  rt <- rnorm(1, 12, 10)
  rt.c <- ifelse(rt > 36, 36, rt)
  rt.c <- ifelse(rt.c < 2, 2, rt.c)
  REGIST[r:(r+2), ] <- log(rt.c)
  
  ##月あたりのプレイ時間の対数
  play <- rnorm(1, 5, 6)
  play.c <- ifelse(play < 0, abs(play)+1, play)
  PLAY[r:(r+2), ] <- log(play.c)

  ##好きなメンバーを決定
  #効用関数を定義
  Honoka <- betaho + CLOTH.ho[r:(r+2), 1:cl] %*% betacs + SoV[r:(r+2), 1] * betasv + CENTER[r:(r+2), 1] * betace +
            REGIST[r:(r+2)] * betatm[1] + PLAY[r:(r+2)] * betapt[1] 
  
  Kotori <- betaho + CLOTH.ko[r:(r+2), 1:cl] %*% betacs + SoV[r:(r+2), 2] * betasv + CENTER[r:(r+2), 2] * betace +
            REGIST[r:(r+2)] * betatm[2] + PLAY[r:(r+2)] * betapt[2] 
  
  Umi <- betaho + CLOTH.um[r:(r+2), 1:cl] %*% betacs + SoV[r:(r+2), 3] * betasv + CENTER[r:(r+2), 3] * betace +
         REGIST[r:(r+2)] * betatm[3] + PLAY[r:(r+2)] * betapt[3] 
  
  Maki <- betaho + CLOTH.ma[r:(r+2), 1:cl] %*% betacs + SoV[r:(r+2), 4] * betasv + CENTER[r:(r+2), 4] * betace +
          REGIST[r:(r+2)] * betatm[4] + PLAY[r:(r+2)] * betapt[4] 
  
  Hanayo <- betaho + CLOTH.ha[r:(r+2), 1:cl] %*% betacs + SoV[r:(r+2), 5] * betasv + CENTER[r:(r+2), 5] * betace +
            REGIST[r:(r+2)] * betatm[5] + PLAY[r:(r+2)] * betapt[5] 
  
  Rin <- betaho + CLOTH.ri[r:(r+2), 1:cl] %*% betacs + SoV[r:(r+2), 6] * betasv + CENTER[r:(r+2), 6] * betace +
         REGIST[r:(r+2)] * betatm[6] + PLAY[r:(r+2)] * betapt[6] 
  
  Nozomi <- betaho + CLOTH.no[r:(r+2), 1:cl] %*% betacs + SoV[r:(r+2), 7] * betasv + CENTER[r:(r+2), 7] * betace +
            REGIST[r:(r+2)] * betatm[7] + PLAY[r:(r+2)] * betapt[7]  
  
  Eri <- betaho + CLOTH.er[r:(r+2), 1:cl] %*% betacs + SoV[r:(r+2), 8] * betasv + CENTER[r:(r+2), 8] * betace +
         REGIST[r:(r+2)] * betatm[8] + PLAY[r:(r+2)] * betapt[8] 
  
  Nico <- betaho + CLOTH.ni[r:(r+2), 1:cl] %*% betacs + SoV[r:(r+2), 9] * betasv + CENTER[r:(r+2), 9] * betace +
          REGIST[r:(r+2)] * betatm[9] + PLAY[r:(r+2)] * betapt[9] 
  
  Mob <- SoV[r:(r+2), 10] * betasv + CENTER[r:(r+2), 10] * betace
  
  #効用関数の行列を作成
  u <- c(exp(Honoka), exp(Kotori), exp(Umi), exp(Maki), exp(Hanayo), exp(Rin), exp(Nozomi), 
         exp(Eri), exp(Nico), exp(Mob))
  U[r:(r+2), ] <- u
}

##選択結果を発生させる
name <- c("Honoka", "Kotori", "Umi", "Maki", "Hanayo", "Rin", "Nozomi", "Eri", "Nico", "Mob")
colnames(U) <- name

##1位選択を発生
#確率を計算
U1 <- U[RANK==1, ]   #1位の効用を取得
Pr1 <- U1 / rowSums(U1)   #確率を計算

#多項乱数で1位選択を発生させる
First <- t(apply(Pr1, 1, function(x) rmultinom(1, 1, x)))
colnames(First) <- name
(first_sum <- colSums(First))   #メンバーごとの1位に選ばれた回数

##2位選択を発生
#確率を計算
U2 <- abs(First-1) * U[RANK==2, ]   #1位で選択したメンバーを選択肢から除く
Pr2 <- U2 / rowSums(U2)

#多項乱数で2位選択を発生させる
Second <- t(apply(Pr2, 1, function(x) rmultinom(1, 1, x)))
colnames(Second) <- name
(second_sum <- colSums(Second))   #メンバーごとの2位に選ばれた回数

##3位選択を発生
#確率を計算
U3 <- abs(First+Second-1) * U[RANK==3, ]   #1位と2位で選択したメンバーを選択肢から除く
Pr3 <- U3 / rowSums(U3)

#多項乱数で3位選択を発生させる
Third <- t(apply(Pr3, 1, function(x) rmultinom(1, 1, x)))
colnames(Third) <- name
(third_sum <- colSums(Third))   #メンバーごとの3位に選ばれた回数


####発生させた順位と説明変数の要約集計####
#順位の要約
(rank_sum <- rbind(first_sum, second_sum, third_sum))   #メンバーごとに選ばれた順位の集計
colSums(rank_sum)   #1～3位まで選ばれた累計数
round(Pr1, 3)   #メンバーごとの1位に選ばれる確率
round(Pr2, 3)   #メンバーごとの2位に選ばれる確率
round(Pr3, 3)   #メンバーごとの3位に選ばれる確率
round(colMeans(Pr1), 3)   #メンバーごとの1位に選ばれる平均確率
round(apply(Pr2, 2, function(x) mean(x[x!=0])), 3)   #メンバーごとの2位に選ばれる平均確率
round(apply(Pr3, 2, function(x) mean(x[x!=0])), 3)   #メンバーごとの3位に選ばれる平均確率

#説明変数の要約
round(colMeans(SoV[RANK==1, ]), 2)   #メンバーごとのshare of voice
round(colMeans(CENTER[RANK==1, ]), 2)   #メンバーごとのセンター設定回数
round(mean(exp(PLAY[RANK==1, ])), 3)   #プレイ時間平均
round(mean(exp(REGIST[RANK==1, ])), 3)   #登録からの経過時間平均


####ランクロジットモデルで推定####
##対数尤度を定義
fr <- function(b, RANK, First, Second, Third, CLOTH.ho, CLOTH.ko, CLOTH.um, CLOTH.ma, CLOTH.ha, 
               CLOTH.ri, CLOTH.no, CLOTH.er, CLOTH.ni, SoV, CENTER, REGIST, PLAY, member, cl){
  #パラメータの設定
  betaho <- b[1]   #穂乃果ちゃんの回帰係数
  betako <- b[2]   #ことりちゃんの回帰係数
  betaum <- b[3]   #海未ちゃんの回帰係数
  betama <- b[4]   #真姫ちゃんの回帰係数
  betaha <- b[5]   #かよちんの回帰係数
  betari <- b[6]   #凛ちゃんの回帰係数
  betano <- b[7]   #のんたんの回帰係数
  betaer <- b[8]   #えりちかの回帰係数
  betani <- b[9]   #にこにーの回帰係数
  betacs <- b[10:(10+cl-1)]   #衣装の回帰係数
  betasv <- b[(10+cl)]   #share of voiceの回帰係数
  betace <- b[(10+cl+1)]   #メンバーごとのセンター設定回数
  betatm <- b[(10+cl+2):(10+cl+2+member-2)]   #登録からの経過時間の回帰係数
  betapt <- b[(10+cl+2+member-1):(10+cl+2+2*member-2)]   #月あたりの平均プレイ時間の回帰係数
  
  
  #メンバーごとの効用関数を定義
  Honoka <- betaho + CLOTH.ho[, 1:cl] %*% betacs + SoV[, 1] * betasv + CENTER[, 1] * betace +
            REGIST * betatm[1] + PLAY * betapt[1] 
  
  Kotori <- betaho + CLOTH.ko[, 1:cl] %*% betacs + SoV[, 2] * betasv + CENTER[, 2] * betace +
            REGIST * betatm[2] + PLAY * betapt[2] 
  Umi <- betaho + CLOTH.um[, 1:cl] %*% betacs + SoV[, 3] * betasv + CENTER[, 3] * betace +
         REGIST * betatm[3] + PLAY * betapt[3] 
  
  Maki <- betaho + CLOTH.ma[, 1:cl] %*% betacs + SoV[, 4] * betasv + CENTER[, 4] * betace +
          REGIST * betatm[4] + PLAY * betapt[4] 
  
  Hanayo <- betaho + CLOTH.ha[, 1:cl] %*% betacs + SoV[, 5] * betasv + CENTER[, 5] * betace +
            REGIST * betatm[5] + PLAY * betapt[5] 
  
  Rin <- betaho + CLOTH.ri[, 1:cl] %*% betacs + SoV[, 6] * betasv + CENTER[, 6] * betace +
         REGIST * betatm[6] + PLAY * betapt[6] 
  
  Nozomi <- betaho + CLOTH.no[, 1:cl] %*% betacs + SoV[, 7] * betasv + CENTER[, 7] * betace +
            REGIST * betatm[7] + PLAY * betapt[7]  
  
  Eri <- betaho + CLOTH.er[, 1:cl] %*% betacs + SoV[, 8] * betasv + CENTER[, 8] * betace +
         REGIST * betatm[8] + PLAY * betapt[8] 
  
  Nico <- betaho + CLOTH.ni[, 1:cl] %*% betacs + SoV[, 9] * betasv + CENTER[, 9] * betace +
          REGIST * betatm[9] + PLAY * betapt[9] 
  
  Mob <- SoV[, 10] * betasv + CENTER[, 10] * betace
  
  U
  #効用関数の行列を定義
  U <- cbind(exp(Honoka), exp(Kotori), exp(Umi), exp(Maki), exp(Hanayo), exp(Rin), exp(Nozomi),
             exp(Eri), exp(Nico), exp(mob))
  
  #順位ごとに対数尤度を定義
  #1位選択
  U1 <- U[RANK==1, ]   #1位選択の効用を定義
  LLF <- log((U1 / rowSums(U1))^First)   #対数尤度を計算
  
  #2位選択
  U2 <- abs(First-1) * U[RANK==2, ]   #1位で選択したメンバーを選択肢から除く
  LLS <- log((U2 / rowSums(U2))^Second)   #対数尤度を計算
  
  #3位選択
  U3 <- abs(First+Second-1) * U[RANK==3, ]   #1位と2位で選択したメンバーを選択肢から除く
  LLT <- log((U3 / rowSums(U3))^Third)   #対数尤度を計算
  
  #対数尤度の和
  LL <- sum(LLF+LLS+LLT)
  return(LL)
}
