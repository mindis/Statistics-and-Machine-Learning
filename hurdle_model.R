#####�[���ߏ�|�A�\����A���f��(�n�[�h�����f��)#####
library(MASS)
library(pscl)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
N <- 5000   #�T���v����
N.zero <- 2000   #�������^�̃[���̃T���v���� 
N.pois <- 3000   #���������̃T���v����
k <- 15   #�����ϐ���

##�����ϐ��̔���
#�A���ϐ��̔���
cont <- 10   #�A���ϐ��̐����ϐ���
X.cont <- matrix(rnorm(N*cont, 0, 1), N, cont)  

#��l�ϐ��̔���
bin <- 5   #��l�ϐ��̐����ϐ���
X.bin <- matrix(0, N, bin)
for(i in 1:bin){
  Pb <- runif(1, 0.2, 0.8)
  X.bin[, i] <- rbinom(N, 1, Pb)
}

#�f�[�^�̌���
X <- data.frame(cont=X.cont, bin=X.bin)

##��A�W���̐ݒ�
#�|�A�\����A�̉�A�W��
betap0 <- 0.68
betap.c <- runif(cont, 0.1, 0.5)
betap.b <- runif(bin, -0.3, 0.65)
betap <- c(betap.c, betap.b)

##�|�A�\���������牞���ϐ��𔭐�
lambda <- exp(betap0 + as.matrix(X[1:N.pois, ]) %*% c(betap))
lambda.all <- exp(betap0 + as.matrix(X) %*% c(betap))
y <- rpois(N.pois, lambda)

#�|�A�\����������̉����ϐ��ƃ[���̉����ϐ�������
Y <- c(y, rep(0, N.zero))
zz <- c(rep(1, N.pois), rep(0, N.zero))

####EM�A���S���Y���Ńn�[�h�����f���𐄒�####
##���S�f�[�^�ł̃|�A�\����A���f���̑ΐ��ޓx
fr <- function(b, X, Y, k, zpt){
  beta0 <- b[1]
  beta <- b[2:(k+1)]
  
  #�ޓx���`���Ęa�����
  #�|�A�\����A���f���̕��ύ\��
  lambda <- exp(beta0 + as.matrix(X) %*% beta)
  
  LLpois <- Y*log(lambda)-lambda - lfactorial(Y)  #�|�A�\�����f���̑ΐ��ޓx
  LLzero <- dpois(Y, 0)   #�[���̖ޓx
  LLzeros <- log(ifelse(LLzero==0, 10^(-300), LLzero))   #0���������ޓx�ɒu������
  LL <- sum(zpt * cbind(LLpois, LLzeros))   #���݊m���ŏd�݂������ΐ��ޓx�̘a�����
  return(LL)
}


##�ϑ��f�[�^�ł̖ޓx�Ɛ��ݕϐ�z�̌v�Z
obsll <- function(x, X, Y, r, N, k){
  beta0 <- x[1]
  beta <- x[2:(k+1)]
  
  #�ޓx���`���Ęa�����
  lambda <- exp(beta0 + as.matrix(X) %*% beta)
  
  #�ޓx�Ƒΐ��ޓx���v�Z
  LLpois <- Y*log(lambda)-lambda - lfactorial(Y)   #�|�A�\�����f���̑ΐ��ޓx
  LLzero <- dpois(Y, 0)   #�[���̖ޓx
  LLzeros <- log(ifelse(LLzero==0, 10^(-300), LLzero))   #0���������ޓx�ɒu������
  
  LLe <- exp(cbind(LLpois, LLzeros))   #�ΐ��ޓx��ޓx�ɖ߂�
  LLe2 <- ifelse(LLe < 10^(-300), 10^(-300), LLe)   #�ޓx��0�̉ӏ��͏������ޓx�ɒu������
  
  #�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̌v�Z
  #������
  R <- matrix(r, N, 2, byrow=T)
  
  #���ݕϐ�z�̌v�Z
  LLr <- R * LLe2
  z0 <- matrix(apply(LLr, 1, sum), N, 2)   #z�̕���
  z1 <- LLr / z0   #z�̌v�Z
  
  #�ϑ��f�[�^�̑ΐ��ޓx
  LLobz <- sum(log(apply(matrix(r, N, 2, byrow=T) * LLe2, 1, sum)))   #�ϑ��f�[�^�ł̑ΐ��ޓx
  rval <- list(LLobz=LLobz, z1=z1)
  return(rval)
}

##EM�A���S���Y���̏����l�̐ݒ�
iter <- 0

#�p�����[�^�̏����l�̐ݒ�
fs <- glm(Y ~ as.matrix(X), family=poisson)
beta <- as.numeric(fs$coef) + runif(length(fs$coef), -0.25, 0.25)
r <- c(0.5, 0.5)   #�������̏����l

#�ϑ��f�[�^�̖ޓx�Ɛ��ݕϐ�z�̏����l
obsllz <- obsll(x=beta, X=X, Y=Y, r=r, N=N, k=k)
LL1 <- obsllz$LLobz
z <- obsllz$z1

dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̏����l��ݒ�
tol <- 0.01

####EM�A���S���Y��####
##���S�f�[�^�ł̃n�[�h�����f���̉�A�W���𐄒�(M�X�e�b�v)
while(dl >= tol){   #dl��tol�ȏ�Ȃ�J��Ԃ�
  res <- optim(beta, fr, X=X, Y=Y, k=k, zpt=z, method="BFGS", 
               hessian=FALSE, control=list(fnscale=-1))
  
  beta <- res$par
  r <- apply(z, 2, sum) / N   #�������̌v�Z
  
  ##E�X�e�b�v
  obsllz <- obsll(x=beta, X=X, Y=Y, r=r, N=N, k=k)
  LL <- obsllz$LLobz
  z <- obsllz$z1
  
  iter <- iter+1
  dl <- abs(LL- LL1)
  LL1 <- LL
  print(LL)
}

####���茋�ʂƓ��v��####
round(beta, 3)   #���肳�ꂽbeta
round(c(betap0, betap), 3)   #�^��beta

round(as.numeric(r), 3)   #���肳�ꂽ������
round(c(N.pois/N, N.zero/N), 3)   #�^�̍�����

#�^�̌��ʂƂ̔�r
YZ <- round(data.frame(Zt=zz, Y, pois=lambda.all, z1=z[, 1], z0=z[, 2]), 3)

#�K���x
(LL <- obsllz$LLobz)   #�ϑ��f�[�^�̑ΐ��ޓx
-2*(LL) + 2*(length(beta)+length(r))   #AIC