#####�ړI�ϐ��Ɍ������܂񂾃��W�b�g���f��#####
library(MASS)
library(bayesm)
library(R2WinBUGS)
library(rstan)
library(reshape2)
library(plyr)
library(lattice)
library(ggplot2)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
col <- 15   #�p�����[�^��
N <- 4000   #�T���v����

##�����ϐ��̔���
#�A���ϐ��̔���
cont <- 7   #�A���ϐ��̃p�����[�^��
X.cont <- matrix(rnorm(N*cont, -1, 1), N, cont)

#��l�ϐ��̔���
bin <- 3   #��l�ϐ��̃p�����[�^��
X.bin <- matrix(0, N, bin)
for(i in 1:bin){
  r <- runif(1, 0.2, 0.8)
  X.bin[, i] <- rbinom(N, 1, r)
}

#���l�ϐ��̔���
multi <- 5   #���l�ϐ��̃p�����[�^��
m <- runif(5)
X.ma <- t(rmultinom(N, 1, m))
zm <- which.min(colSums(X.ma))
X.multi <- X.ma[, -zm]

#�f�[�^�̌���
round(X <- data.frame(cont=X.cont, bin=X.bin, multi=X.multi), 2)

##��A�W���̐ݒ�
alpha0 <- 0.6
beta.cont <- runif(cont, 0, 0.65)
beta.bin <- runif(bin, -0.7, 0.9)
beta.multi <- runif(multi-1, -0.8, 1.2)
betaT <- c(alpha0, beta.cont, beta.bin, beta.multi)

##�����ϐ��̔���
#�m���̌v�Z
logit <- alpha0 + as.matrix(X) %*% betaT[-1]   #���W�b�g
P <- exp(logit)/(1+exp(logit))   #�m���̌v�Z
hist(P, col="grey", main="�m���̕��z")

#�x���k�[�C�����ŉ����ϐ��𔭐�
Y <- rbinom(N, 1, P)
round(cbind(Y, P), 2)   #�����ϐ��Ɗm���̔�r


##���t�f�[�^�������_���Ȍ����ɏ]���č폜����
#�����_���Ȍ���������
alpha.na <- 0.8
beta1.na <- runif(cont, 0, 0.7)
beta2.na <- runif(bin, -0.8, 1.4)

#���W�b�g�Ɗm���̌v�Z
logit.na <- alpha.na + as.matrix(X[, 1:cont]) %*% beta1.na + as.matrix(X.bin) %*% beta2.na
P.na <- exp(logit.na)/(1+exp(logit.na))

#�x���k�[�C�����ŉ����ϐ�������������
z.na <- rbinom(N, 1, P.na)

#�����x�N�g���̍쐬
Y.na <- ifelse(z.na==1, NA, Y)
YZ <- data.frame(Y, z.na, Y.na)

####EM�A���S���Y���Ō����̂���f�[�^�𐄒�####
##���W�X�e�B�b�N��A���f���̑ΐ��ޓx���`
loglike <- function(b, X, Y){
  #�p�����[�^�̐ݒ�
  alpha <- b[1]
  beta <- b[2:(col)]
  
  #�ޓx���`���č��v����
  logit <- alpha + as.matrix(X) %*% as.vector(beta) 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}


b <- beta
Z <- z
##���S�f�[�^�ł̃��W�X�e�B�b�N��A���f���̑ΐ��ޓx
fr <- function(b, X.comp, X.na, Y.comp, Z, k){
  beta0 <- b[1]
  beta <- b[2:(k+1)]
  
  #���S�f�[�^�̑ΐ��ޓx
  logit.comp <- beta0 + as.matrix(X.comp) %*% beta 
  P.comp <- exp(logit.comp) / (1 + exp(logit.comp))
  LLS.comp <- Y.comp*log(P.comp) + (1-Y.comp)*log(1-P.comp)  
  LL.comp <- sum(LLS.comp)
  
  #�s���S�f�[�^�̑ΐ��ޓx
  logit.na <- beta0 + as.matrix(X.na) %*% beta
  P.na <- exp(logit.na) / (1 + exp(logit.na))
  LLs.na <- Z[, 1]*log(P.na) + Z[, 2]*log(1-P.na)
  LL.na <- sum(LLs.na)
  LL <- LL.comp + LL.na
  return(LL)
}

##�ϑ��f�[�^�ł̖ޓx�Ɛ��ݕϐ�z�̌v�Z
obsll <- function(b, X.comp, Y.comp, X.na, r, N.na, k){
  beta0 <- b[1]
  beta <- b[2:(k+1)]
  
  #���S�f�[�^�̑ΐ��ޓx
  logit.comp <- beta0 + as.matrix(X.comp) %*% beta 
  P.comp <- exp(logit.comp) / (1 + exp(logit.comp))
  LLS.comp <- Y.comp*log(P.comp) + (1-Y.comp)*log(1-P.comp)  
  LL.comp <- sum(LLS.comp)
  
  #�s���S�f�[�^�̖ޓx���v�Z
  logit.na <- beta0 + as.matrix(X.na) %*% beta
  P.na <- exp(logit.na) / (1 + exp(logit.na))
  
  #�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̌v�Z
  #������
  R <- matrix(r, N.na, 2, byrow=T)
  
  #���ݕϐ�z�̌v�Z
  LLr <- R * cbind(P.na, 1-P.na)
  z0 <- matrix(apply(LLr, 1, sum), N.na, 2)   #z�̕���
  z1 <- LLr/z0   #z�̌v�Z
  
  #�ϑ��f�[�^�̑ΐ��ޓx
  LL.na <- sum(log(apply(matrix(r, N.na, 2, byrow=T) * cbind(P.na, 1-P.na), 1, sum)))
  LLobz <- LL.na
  rval <- list(LLobz=LLobz, z1=z1)
  return(rval)
}

##�f�[�^�̐ݒ�
#�f�[�^�����S�f�[�^�ƌ����f�[�^�Ƀ\�[�g
index.na <- subset(1:nrow(YZ), is.na(YZ$Y.na)==TRUE)

#���S�f�[�^
YZ.comp <- YZ[-index.na, ]
Y.comp <- Y[-index.na]
X.comp <- X[-index.na, ]

#�����̂���f�[�^
YZ.na <- YZ[index.na, ]
Y.na <- Y[index.na]
X.na <- X[index.na, ]

#���������ϐ�
YZ <- rbind(YZ.comp, YZ.na)


##EM�A���S���Y���̐ݒ�Ə����l�̐ݒ�
iter <- 0
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̏����l��ݒ�
tol <- 0.1

#�p�����[�^�̏����l�̐ݒ�
#�ΐ��ޓx���ő剻����
b0 <- c(rep(0, ncol(X.comp)+1))   #�����p�����[�^�̐ݒ�
res <- optim(b0, loglike, gr=NULL, X=X.comp, Y=Y.comp, method="BFGS", hessian=TRUE, control=list(fnscale=-1))
beta <- res$par + runif(length(res$par), 0.25, 0.25)
r <- c(0.5, 0.5)   #�������̏����l

#�ϑ��f�[�^�̖ޓx�Ɛ��ݕϐ�z�̏����l
obsllz <- obsll(b=beta, X.comp=X.comp, Y.comp=Y.comp, X.na=X.na, r=r, N.na=nrow(X.na), k=ncol(X.comp))
LL1 <- obsllz$LLobz
z <- obsllz$z1
z

####EM�A���S���Y���ɂ��s���S�f�[�^���W�X�e�B�b�N��A���f��####
##���S�f�[�^�ł̉�A�W���𐄒�(M�X�e�b�v)
while(dl >= tol){   #dl��tol�ȏ�Ȃ�J��Ԃ�
  res <- optim(beta, fr, X.comp=X.comp, Y.comp=Y.comp, X.na=X.na, Z=z, k=ncol(X.comp), method="BFGS", 
               hessian=FALSE, control=list(fnscale=-1))
  
  beta <- res$par
  r <- apply(z, 2, sum) / nrow(X.na)   #�������̌v�Z
  
  ##E�X�e�b�v
  obsllz <- obsll(b=beta, X.comp=X.comp, Y.comp=Y.comp, X.na=X.na, r=r, N.na=nrow(X.na), k=ncol(X.comp))
  LL <- obsllz$LLobz
  z <- obsllz$z1
  
  iter <- iter+1
  dl <- abs(LL- LL1)
  LL1 <- LL
  print(LL)
}

round(beta, 3)
round(betaT, 3)

round(cbind(YZ.na, z), 3)