#####�}���R�t���f���ɂ��Web�T�C�g��V�s������#####
library(MASS)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(93417)
N <- 2000   #�T���v����
S <- 5   #�T�C�g�y�[�W��

##�T�C�g��V�m���Ɨ��E�m���̐ݒ�
#�}���R�t���ڍs��̐ݒ�
p1 <- c(0.45, 0.2, 0.1, 0.05, 0.2)
p2 <- c(0.3, 0.3, 0.1, 0.2, 0.1)
p3 <- c(0.3, 0.2, 0.2, 0.1, 0.2)
p4 <- c(0.3, 0.3, 0.1, 0.15, 0.15)
p5 <- c(0.2, 0.3, 0.1, 0.2, 0.2)
p <- rbind(p1, p2, p3, p4, p5)   #�f�[�^�̌���
pt <- p

#�����f�B���O�y�[�W�̊m���̐ݒ�
pf <- c(0.5, 0.3, 0.1, 0.05, 0.05)
pft <- pf

#���E�m���̊֐��̐ݒ�
surv <- function(beta, X){
  logit <- beta[1] + as.matrix(X) %*% beta[-1]
  pr <- exp(logit)/(1+exp(logit))
  return(pr)
}

#���E���f���̉�A�W���̐ݒ�
alpha <- -2.0
beta.time <- 0.10
beta.page <- c(-0.6, -0.4, 0.15, 0.30)
beta <- c(alpha, beta.time, beta.page)
betat <- beta

x <- t(rmultinom(30, 1, runif(5)))
X <- cbind(1:30, x[, -3])
surv(beta, X)


##���U���ԃ}���R�t�������f���ŃV�~�����[�V�����f�[�^�𔭐�
#�f�[�^�̊i�[�p�I�u�W�F�N�g�̐ݒ�
X <- matrix(0, 200000, S+1)
Z <- c()
PZ <- c()

for(i in 1:N){
#�����f�B���O�y�[�W������
  print(i)
  r <- sum(rowSums(X)!=0)
  l <- t(rmultinom(1, 1, pf))
  x <- t(c(1, l[-3])) 
  X[r+1, ] <- t(c(1, l))
  
  #2��ڈȍ~�T�C�g��V
  M <- matrix(0, 100, S)
  
  for(t in 2:1000){
    #���E�������ǂ���������
    pz <- surv(beta=beta, X=x)   #���E�m��
    z <- rbinom(1, 1, pz)   #���E�������ǂ����̌���
    PZ <- c(PZ, pz)
    Z <- c(Z, z)
    if(z==1) break   #z=1�ŗ��E
    
    #�ǂ̃y�[�W�ɐi�񂾂�������
    l <- (1-z) * (rmultinom(1, 1, p[which.max(l), ]))   #�}���R�t����
    M[1, ] <- t(l)
    x <- t(c(t, l[-3]))
    X[r+t, ] <- t(c(t, l))
  }
}

##�s�v�ȃf�[�^���폜����
index <- rowMeans(X)!=0   #�]���Ă���f�[�^���������
X <- X[index, ]
summary(X)   #�f�[�^�̗v��

##ID�̐ݒ�
id <- c()
for(i in 1:N){
  r <- rep(i, subset(1:length(Z), Z==1)[i]-length(id))
  id <- c(id, r)  
}

##���ׂẴf�[�^������
ZX <- data.frame(id, Pr=PZ, Z, t=X[, 1], X=X[, -1])
round(ZX, 2)


####�Ŗޖ@�ŗ��U���Ԑ����}���R�t���ڃ��f���𐄒�####
#�}���R�t�������藧�̂ŁA�����f�B���O�y�[�W�A�y�[�W�J�ځA���E�͂��ꂼ��Ɨ��Ɖ���ł���B
#���������āA3�̃p�����[�^�͓Ɨ��ɋ��߂�B

##�y�[�W���ꂼ��̃����f�B���O�m�������߂�
Pl <- colMeans(X[X[, 1]==1, 2:ncol(X)])
round(Pl, 3)


##�y�[�W���ڊm���s������߂�
#�����f�B���O�y�[�W�ŗ��E���Ă��郆�[�U�[����菜��
index1 <- subset(1:N, as.numeric(table(id))!=1)
R <- c()
for(i in 1:length(index1)){
  r <- subset(1:length(id), id==index1[i])
  R <- c(R, r)
}
Xm <- X[R, ]   #�y�[�W���ڂ�����f�[�^�̂ݎ��o��

#�}���R�t���ڍs������߂�
Pr <- matrix(0, 0, S) 
for(s in 1:S){
  index1 <- subset(1:nrow(Xm), Z[R]!=1 & Xm[, s+1]==1)   #Z=1�̃f�[�^����菜���At-1�̃f�[�^�����o��
  Xm1 <- Xm[index1, 2:ncol(Xm)]   #t�̗����菜��
  Nm <- sum(Xm1[, s])   #t-1��s�̃f�[�^�����v����
  Mm <- colSums(Xm[index1+1, 2:ncol(Xm)])   #t�̃f�[�^�����v����
  pr <- Mm / Nm   #�m�����v�Z
  Pr <- rbind(Pr, pr)
}
rownames(Pr) <- c("P1", "P2", "P3", "P4", "P5")
round(Pr, 3)


##���U���Ԑ������f���ŗ��E���𐄒�
#���U���Ԑ������f���̑ΐ��ޓx��ݒ�
loglike <- function(beta, y, X){
  logit <- beta[1] + as.matrix(X) %*% beta[-1]    #���W�b�g�̌v�Z
  p <- exp(logit)/(1+exp(logit))   #�m���̌v�Z
  
  #�ΐ��ޓx�̌v�Z
  LLs <- y*log(p) + (1-y)*log(1-p)
  LL <- sum(LLs)
  return(LL)
}

#�ΐ��ޓx���ő剻����
b0 <- c(-1.0, 0.2, runif(S-1, -1, 1))
res <- optim(b0, loglike, y=Z, X=X[, -4], method="BFGS", hessian=TRUE, control=list(fnscale=-1))

#���茋��
beta <- res$par
round(beta, 3)   #���肳�ꂽ��A�W��
round(betat, 3)   #�^�̉�A�W��


##���ׂĂ̐��茋�ʂ�\��
#�������茋�ʁA�E���^�̒l
round(Pl, 2); round(pft, 3)   #�y�[�W�ʂ̃����f�B���O�m��
round(Pr, 2); round(pt, 2)   #�}���R�t���ڍs��
round(beta, 2); round(betat, 2)   #�T�C�g���E���f���̃p�����[�^


####���E���̃V�~�����[�V����####
n <- 500   #�V�~�����[�V��������T���v����
t <- 100   #�V�~�����[�V�����������

#�f�[�^�̊i�[�p�I�u�W�F�N�g�̐ݒ�
Xs <- matrix(0, n*t, S+1)
Zs <- c()
PZs <- c()
for(i in 1:n){
  #�����f�B���O�y�[�W������
  print(i)
  r <- sum(rowSums(Xs)!=0)
  l <- t(rmultinom(1, 1, Pl))
  x <- t(c(1, l[-3])) 
  Xs[r+1, ] <- t(c(1, l))
  
  #2��ڈȍ~�T�C�g��V
  M <- matrix(0, 100, S)
  
  for(t in 2:t){
    #���E�������ǂ���������
    p <- surv(beta=beta, X=x)   #���E�m��
    z <- rbinom(1, 1, p)   #���E�������ǂ����̌���
    PZs <- c(PZs, p)
    Zs <- c(Z, z)
    
    #�ǂ̃y�[�W�ɐi�񂾂�������
    l <- (rmultinom(1, 1, Pr[which.max(l), ]))   #�}���R�t����
    M[1, ] <- t(l)
    x <- t(c(t, l[-3]))
    Xs[r+t, ] <- t(c(t, l))
  }
  #�ŏI���Ԃ̗��E��
  p <- surv(beta=beta, X=x)   #���E�m��
  z <- rbinom(1, 1, p)   #���E�������ǂ����̌���
  PZs <- c(PZs, p)
  Zs <- c(Z, z)
}

#�y�[�W���ڐ��݂̂��l���������E��
s <- matrix(0, t, S-1)
Xt <- data.frame(1:t, s)
X.surv <- surv(beta=beta, X=Xt)

##�f�[�^������
ids <- rep(1:n, rep(t, n))   #id��ݒ�
ZXs <- data.frame(id=ids, p=PZs, t=Xs[, 1], S=Xs[, -1])   #�f�[�^������
round(ZXs, 3)


##���E��������
plot(1:t, X.surv, type="l", ylab="���E��", xlab="�y�[�W���ڐ�", lwd=3, col=2, ylim=c(0.1, 1))
for(i in 1:5){
  lines(1:t, ZXs[ZXs$id==i, ]$p, type="l", lty=i)
}
lines(1:t, X.surv, type="l", ylab="���E��", xlab="�y�[�W���ڐ�", lwd=3, col=2, ylim=c(0.1, 1))