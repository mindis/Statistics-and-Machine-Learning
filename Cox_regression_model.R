####�J��Ԃ��̂Ȃ��R�b�N�X��A���f��(�Z�~�p�����g���b�N�������f��)####
library(MASS)
library(survival)
library(reshape2)
library(plyr)

####�f�[�^�̔���####
#set.seed(9843)
t <- 100   #�ϑ�����
n <- 1000   #�T���v����
col <- 10   #�ϐ���

##�x�[�X���C���n�U�[�h�̐ݒ�
tb <- -4.7   #�����l
trend <- numeric()
s <- seq(0.5, 0.9, length=t)
for(i in 1:t){
  r <- rnorm(5, tb, 0.07)
  sort <- sort(r)
  bi <- rbinom(1, 1, s[i])
  bb <- ifelse(bi == 1, sort[4], sort[2])
  tb <- bb
  trend <- c(trend, bb)
}
plot(1:length(trend), trend, type="l", xlab="time", lwd=2)   #�g�����h�̃v���b�g
round(exp(trend)/(1+exp(trend)), 3) 

##�����ϐ��̔����Ɖ�A�����̐ݒ�
##�ÓI�A���ϐ��ƐÓI��l�ϐ��̔���
#�ÓI�ϐ��̐�
s.cont <- 5
s.bin <- 2

#�ÓI�A���ϐ��̔���
X.scont <- matrix(rnorm(n*s.cont, 0, 1), nrow=n, ncol=s.cont)

#���I��l�ϐ��̔���
X.sbin <- matrix(0, n, s.bin)
pr_b <- runif(s.bin, 0.3, 0.7)
for(i in 1:s.bin){
  bin <- rbinom(n, 1, pr_b[i])
  X.sbin[, i] <- bin
}

##���I�A���ϐ��Ɠ��I��l�ϐ��̔���
#���I�ϐ�
d.cont <- 1
d.bin <- 2

##���I�A���ϐ��̔���
type1 <- 10   #�A���ϐ�(�O���ϐ�)�̃^�C�v��
X.d <- matrix(0, t, 10)
for(i in 1:type1){
  X.d[, i] <- matrix(rnorm(t, 0, 1), t, 1)
}

#���I�A���ϐ��̌l���ƂɊ��蓖�Ă�
r <- runif(type1, 0.2, 0.7)
X.dcont <- matrix(0, n, t)
X.dd <- t(X.d)
for(i in 1:n){
  m <- which.max(t(rmultinom(1, 1, r)))
  X.dcont[i, ] <- X.dd[m, ]
  print(i)
}

##���I��l�ϐ��̔���
type2 <- 10   #��l�ϐ�(��l�ϐ�)�̃^�C�v��
X.db <- matrix(0, t, 10)

for(i in 1:type2){
  X.db[, i] <- matrix(rbinom(t, 1, runif(1, 0.3, 0.7)), t, 1)
}

#���I��l�ϐ����l���ƂɊ��蓖�Ă�
r <- runif(type, 0.2, 0.7)
X.dbin1 <- matrix(0, n, t)
X.dd <- t(X.db)

for(i in 1:n){
  m <- which.max(t(rmultinom(1, 1, r)))
  X.dbin1[i, ] <- X.dd[m, ]
  print(i)
}

##���I�ȓ��I��l�ϐ��̔���
X.dbin2 <- matrix(0, n, t)

for(i in 1:n){
  r <- runif(1, -1.5, 1.5)
  for(j in 1:t){
    p <- exp(r)/(1+exp(r))
    X.dbin2[i, j] <- rbinom(1, 1, p)
    r <- r + rnorm(1, 0, 0.025)
  }
}

##�f�[�^���p�l���f�[�^�`���ɕύX
#ID�Ǝ��Ԃ���ю��ʔԍ�
ID <- rep(1:n, rep(t, n))
time <- rep(1:t, n)
No. <- 1:(n*t)

#�x�[�X���C���n�U�[�h���p�l���`���ɕϊ�
TREND <- rep(trend, n)

#�ÓI�ϐ����p�l���`���ɕϊ�
X.s <- matrix(0, nrow=n*t, ncol=s.cont+s.bin)

for(i in 1:n){
  x <- cbind(X.scont, X.sbin)[i, ]
  X.s[ID==i, ] <- matrix(x, nrow=t, ncol=s.cont+s.bin, byrow=T)
}

#���I�ϐ����p�l���`���ɕϊ�
xdc <- matrix(t(X.dcont), nrow=n*t, 1, byrow=T)
xdb1 <- matrix(t(X.dbin1), nrow=n*t, 1, byrow=T)
xdb2 <- matrix(t(X.dbin2), nrow=n*t, 1, byrow=T)

X.d <- cbind(xdc, xdb1, xdb2)   #�f�[�^�̌���

#���ׂẴf�[�^������
round(X <- data.frame(No., ID, time, S=X.s, D=X.d), 2)

##��A�W���̐ݒ�
betasc <- runif(s.cont, -0.15, 0.15) #�ÓI�A���ϐ��̉�A�W��
betasb <- runif(s.bin, -0.2, 0.2)   #�ÓI��l�ϐ��̉�A�W��
betadc <- runif(d.cont, -0.15, 0.15)   #���I�A���ϐ��̉�A�W��
betadb <- runif(d.bin, -0.2, 0.2)   #���I��l�ϐ��̉�A�W��
beta <- c(betasc, betasb, betadc, betadb)   #��A�W��������

##�����N�֐����v�Z
s1 <- which.max(colnames(X)=="S.1")
logit <- TREND + as.matrix(X[, s1:ncol(X)]) %*% beta
Pr <- exp(logit)/(1+exp(logit))

##�C�x���g���Ԃ𔭐�������(�C�x���g������������ł��؂�)
Y <- c()
for(i in 1:n){
  yr <- c()
  for(j in 1:t){
    yr <- c(yr, rbinom(1, 1, Pr[X$ID==i & time==j, ]))
    if(max(yr[1:j])>0) break
  }
  yh <- if(max(yr)>0) which.max(yr)+round(runif(1, 0, 1), 1) else {0}
  Y <- c(Y, yh)
  print(i)
}
hist(Y, breaks=20, col="grey")   #���ʂ̕p�x������
table(Y)

##Cox��A�p�̃f�[�^�Z�b�g���쐬
yt <- rep(Y, rep(t, n))
YXt <- data.frame(X, yt)

#Cox��A�̃C�x���g���Ԃ��Ƃɓ���W�����ƂɃO���[�s���O����
YXl <- list()
for(i in 1:n){
  if(Y[i]==0) next   #�C�x���g���������Ă��Ȃ��Ȃ玟�̃T���v����
  index_y <- trunc(Y[i])   #�C�x���g���Ԃ𒊏o 
  index_x <- subset(1:length(Y), Y[i] <= Y | Y == 0)   #�C�x���g���������Ă��Ȃ��T���v���𒊏o
  YXz <- data.frame(YXt[YXt$ID %in% index_x & YXt$time %in% index_y, ], T=i)   #�C���f�b�N�X�ōi�����ϐ����f�[�^�𒊏o
  YXl[[i]] <- YXz
  print(i)
}
YX <- do.call(rbind, YXl)
index_z <- unique(YX$T)

Z <- list()
for(i in 1:length(index_z)){
  h <- subset(1:length(YX$yt[YX$T==index_z[i]]), YX$yt[YX$T==index_z[i]]>0)
  m <- min(subset(YX$yt[YX$T==index_z[i]], YX$yt[YX$T==index_z[i]]>0))
  z <- as.numeric(YX$yt[YX$T==index_z[i]]==m)
  Z[[i]] <- z
  print(i)
}

##�f�[�^�Z�b�g���܂Ƃ߂�
Zn <- unlist(Z)   #���X�g���x�N�g����
YXz <- data.frame(YX[2:length(YX)], Zn, Z=ifelse(YX$yt>0, 1, 0))   #�f�[�^������
rownames(YXz) <- 1:nrow(YXz)

#�����ϐ����������o��
val <- subset(1:ncol(YXz), colnames(YXz)==c("S.1", "D.3"))
XX <- YXz[, val[1]:val[2]]

####Cox��A���f���𐄒肷��####
##�^�C�f�[�^�����
#�C�x���g���Ԃ��Ƃ̏o�����𐔂���
(tac <- table(Y[Y!=0]))   

#�C�x���g���Ԃ̃��j�[�N����肵�A�����ɕ��בւ�
yu <- unique(YXz$yt[YXz$yt!=0])   
sortlist <- order(yu)   
(y.uni <- yu[sortlist])

sum(as.numeric(rownames(tac))-y.uni)   #���l����v���Ă��邩�`�F�b�N

#�C�x���g���Ԃ̃��j�[�N���ƂɃ^�C�f�[�^���L�^
index_d <- list()
index_m <- list()
for(i in 1:length(y.uni)){
  ZT <- YXz[YXz$Zn==1 & YXz$yt==y.uni[i], "T"]
  ind <- subset(1:nrow(YXz), YXz$T==ZT[1])
  mi <- subset(1:nrow(YXz), YXz$Zn==1 & YXz$yt==y.uni[i])
  index_d[[i]] <- ind
  index_m[[i]] <- mi[1:tac[i]]
  print(i)
}

##�ΐ������ޓx���`
fr <- function(b, D, M, X, uni){
  beta <- b[1:ncol(X)]
  link.f <- exp(as.matrix(X) %*% beta)
  LLs <- c()
  
  #�^�C�f�[�^���܂ޑΐ������ޓx(Efron�@)
  for(i in 1:length(uni)){
    LLd <- log(prod(sum(link.f[D[[i]]]) - (1:length(M[[i]])-1)/length(M[[i]])*sum(link.f[M[[i]]])))
    LLm <- colSums(X[M[[i]], ]) %*% beta
    LLs <- c(LLs, LLm-LLd)
  }
  LL <- sum(LLs)
  return(LL)
}

##�����ޓx���ő剻���āACox��A�𐄒�
D <- index_d   #�������Ԃ�t�̃��X�N�W���C�x���g
M <- index_m   #�C�x���g����t�ɔ��������T���v��
b0 <- runif(10, -0.5, 0.5)   #�p�����[�^�̏����l

#�����ޓx�����j���[�g���@�ōő剻
fit <- optim(b0, fr, gr=NULL, D=D, M=M, X=XX, uni=y.uni,
             method="BFGS", hessian=T, control=list(fnscale=-1))

fit$value   #�ő剻���ꂽ�ΐ��ޓx
round(betan <- fit$par, 2)    #���肳�ꂽ�p�����[�^
round(exp(betan), 2)   #���肳�ꂽ�p�����[�^(�n�U�[�h��)
round(exp(beta), 2)   #�^�̃p�����[�^(�I�b�Y��)

betan/sqrt(-diag(solve(fit$hessian)))   #t�l
(AIC <- -2*fit$value + 2*length(fit$par))   #AIC


####�x�[�X���C���n�U�[�h�֐��Ɛ����֐��̐���####

