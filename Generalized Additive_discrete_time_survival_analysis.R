#####��ʉ����@���U���Ԑ������f��#####
library(MASS)
library(VGAM)
library(mgcv)
library(fda)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(324)
##�f�[�^�̐ݒ�
t <- 24*4   #�����P�ʂ�3�N
n <- 500   #�T���v����
N <- t*n   #���T���v����
col <- 17   #�ϐ���

#id�Atime�A�����ԍ���ݒ�
id <- rep(1:n, rep(t, n))
time <- rep(1:t, n)
no <- 1:N
ID <- data.frame(id, time, no)

##�p�����[�^�̐ݒ�
#�g�����h����
T <- 10000
for(t1 in 1:T){
  tb <- -1.2   #�����l
  trend <- numeric()
  s <- seq(0.85, 0.2, length=t)
  for(i in 1:t){
    r <- rnorm(5, tb, 0.08)
    sort <- sort(r)
    bi <- rbinom(1, 1, s[i])
    bb <- ifelse(bi == 1, sort[4], sort[2])
    tb <- bb
    trend <- c(trend, bb)
  }
  if(max(exp(trend)/(1+exp(trend))) > 0.35 && min(exp(trend)/(1+exp(trend))) < 0.1) break
  print(t1)
}  
plot(exp(trend)/(1+exp(trend)), type="l", lwd=1, xlab="half_month", ylab="p")


##��A�����̐ݒ�
cont <- 8   #�A���ϐ�
cnt <- 2   #�J�E���g�f�[�^��
bin <- 4   #��l�ϐ�

#�A���ϐ��̔���
X_cont <- matrix(runif(cont*N, 0, 1), N, cont)   

#�J�E���g�f�[�^�̔���
X_cnt <- matrix(rpois(cnt*N, 1.5), N, cnt)   

#��l�ϐ��̔���
X_bin <- matrix(0, N, bin)
pr_b <- runif(bin, 0.3, 0.7)
for(i in 1:bin){
  bin <- rbinom(N, 1, pr_b[i])
  X_bin[, i] <- bin
}

#�O��̍w������̌o�ߎ���
X_l <- rep(0, N)
X_q <- rep(0, N)

#�����̐ݒ�
month_l <- (rep(1:t, n)/2)/12
month_q <- month_l^2
month <- data.frame(month_l, month_q)

X <- data.frame(cont=X_cont, cnt=X_cnt, bin=X_bin, hist_l=X_l, hist_q=X_q)

##�p�����[�^�̐ݒ�
beta <- c(runif((col-3), -0.30, 0.40), 2.64, -3.57)
beta0 <- -0.82
beta.t <- beta
beta0.t <- beta0

##�O��̍w������̊��Ԃ̐����ϐ����L�^���Ȃ���A�V�~�����[�V�����f�[�^�𔭐�������
y <- matrix(NA, N, 1)   #�����ϐ����i�[
for(i in 1:n){
  for(j in 1:t){
    r <- (i-1)*t+j
    #�w���f�[�^�̔���
    xb <- trend[j] + beta0 + as.matrix(X[r, ]) %*% beta   #���`�֐�
    y[r, ] <- rbinom(1, 1, exp(xb) / (1 + exp(xb)))   #�w���𔭐�������
    
    #�O�񂩂�̍w���������L�^
    if(id[r]==n && time[r]==t) break   #�ŏI�s��break
    
    if(y[r, 1]==0 && sum(y[is.na(y)==0 & id==i, ])==0 && time[r]!=t) {
      X$hist_l[r+1] <- 0; X$hist_q[r+1] <- 0
    } 
    if(y[r, 1]==1 && time[r]!=t) {
      X$hist_l[r+1] <- 0.5/12; X$hist_q[r+1] <- 0.5^2/12
    } 
    if(y[r, 1]==0 && sum(y[is.na(y)==0 & id==i, ])!=0 && time[r]!=t) {
      h <- subset(1:t, y[id==i, 1]==1); 
      X$hist_l[r+1] <- (((j+1)-h[length(h)])/2)/12; X$hist_q[r+1] <- X$hist_l[r+1]^2
    }
    if(time[r]==t) {
      X$hist_l[r+1] <- 0; X$hist_q[r+1] <- 0
    }
  }
  print(i)
}

##�����������f�[�^�̊m�F�ƏW�v
round(YX <- data.frame(id, time, y, X), 3)   #�f�[�^�̊m�F
mean(y)   #�w����
summary(YX[, 3:length(YX)])   #�f�[�^�̗v��
by(YX[, 3:length(YX)], YX$id, summary)   #�l���Ƃ̗v��
by(YX[, 3:length(YX)], YX$time, summary)   #���Ԃ��Ƃ̗v��
by(YX$y, YX$time, mean)   #���Ԃ��Ƃ̗v��


####���U���Ԑ������f���𐄒�####
##���_�~�[�̓������ΐ��ޓx
#�ΐ��ޓx�̐ݒ�
fr2 <- function(b, month, x, t, y){
  #�p�����[�^�̐ݒ�
  alpha <- b[1]
  betam <- b[2:(t+1)]
  beta <- b[(t+2):length(b)]
  
  #�ޓx���`���č��v����
  Xb <- alpha + as.matrix(month) %*% betam + as.matrix(X) %*% beta 
  p <- exp(Xb) / (1 + exp(Xb))
  LLS <- y*log(p) + (1-y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

##���_�~�[�̓����Ă��Ȃ��ΐ��ޓx
#�ΐ��ޓx�̐ݒ�
fr1 <- function(b, x, y){
  #�p�����[�^�̐ݒ�
  alpha <- b[1]
  beta <- b[(2):length(b)]
  
  #�ޓx���`���č��v����
  Xb <- alpha + as.matrix(X) %*% beta 
  p <- exp(Xb) / (1 + exp(Xb))
  LLS <- y*log(p) + (1-y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

####time�̃_�~�[�ϐ���p���đΐ��ޓx���ő剻����####
#�f�[�^�Ə����p�����[�^�̐ݒ�
#time�̃_�~�[�ϐ�(2��������)���쐬
tm <- list()
for(i in 1:n){
  tm[[i]] <- matrix(rep(c(1, 1, 0), c(1, 1, 96)), t, t/2)
}
TM <- do.call(rbind, tm)   #���X�g���s��ɒ���

para <- ncol(X)+ncol(TM)+1   #�p�����[�^��

##�����l�����j���[�g���@�Ō��肷��
bf0 <- runif((1+ncol(X)))

#�����l�����j���[�g���@�Ō��肷��
res1 <- optim(bf0, fr1, gr=NULL, X, y,
             method="BFGS", hessian=FALSE, control=list(fnscale=-1))
bf <- res1$par   #�p�����[�^�̐���l

##time�̐���l�𐄒肷��
#�����l��ݒ�
b0 <- c(bf[1], runif(t/2, -1, 1), bf[2:length(bf)])

##���j���[�g���@�ŉ���
res <- optim(b0, fr2, gr=NULL, TM, X, t/2, y,
             method="BFGS", hessian=TRUE, control=list(fnscale=-1))

#���肳�ꂽ�p�����[�^
beta <- res$par
round(beta[1], 3); round(beta[(t/2+2):length(beta)], 3)
round(beta0.t, 3); round(beta.t, 3)

(LL <- res$value)   #�ő剻���ꂽ�ΐ��ޓx
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(N)*length(beta))   #BIC

####time�̕����������s���ăg�����h�𐄒�####
##3���X�v���C�����̐ݒ�
##�ߓ_�����߂�

#�p�����[�^�̐���l
beta.trend <- beta[2:(t/2+1)]
plot(1:(t/2), beta.trend, xlab="month", ylab="trend", pch=1)

tp <- seq(1, t/2, length=t/4)   #�ߓ_
bs <- create.bspline.basis(rangeval=c(1, t/2), breaks=tp)   ##�ߓ_����Ǐ��I�ȑ����������
bp1 <- fdPar(bs, Lfdobj=2)
bp1$lambda <- 5   #�������p�����[�^�̐ݒ�

##�����������s���ăO���t���쐬
sm2 <- smooth.basis(1:(t/2), beta.trend, bp1)

#�O���t���쐬
#����l�̃g�����h
par(mfrow=c(1,2)) 
plot(1:(t/2), beta.trend, xlab="month", ylab="trend")
lines(sm2, lwd=2, col="blue")

#�^�̃g�����h
index <- rep(c(1, 0), t/2)
plot(1:(t/2), trend[index==1], type="l", xlab="month", ylab="true_trend", lwd=2)
par(mfrow=c(1,1))
