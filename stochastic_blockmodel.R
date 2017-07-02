#####���N���X�^�����O(�m���I�u���b�N���f��)#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(gtools)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

#set.seed(318)

####�f�[�^�̔���####
#�f�[�^�̐ݒ�
N <- 1000   #���[�U�[��
K <- 150   #�A�C�e����
sg_n <- 4   #���[�U�[�̃Z�O�����g��
sg_k <- 3   #�A�C�e���̃Z�O�����g��

##�p�����[�^�ƃZ�O�����g�����肵�ϑ��f�[�^�𔭐�������
#���[�U�[�Z�O�����g�𔭐�
alpha.s1 <- rep(5, sg_n)
pi.s1 <- rdirichlet(1, alpha.s1)
z.s1 <- t(rmultinom(N, 1, pi.s1))
zno1 <- z.s1 %*% 1:sg_n

#�A�C�e���Z�O�����g�̔���
alpha.s2 <- rep(5, sg_k)
pi.s2 <- rdirichlet(1, alpha.s2)
z.s2 <- t(rmultinom(K, 1, pi.s2))
zno2 <- z.s2 %*% 1:sg_k

#�ϑ��ϐ��̃p�����[�^�̐ݒ�
#���[�U�[�Z�O�����g�~�A�C�e���Z�O�����g�̃x�[�^���O���z�̃p�����[�^�𔭐�
theta.s <- rbeta(sg_n*sg_k, 1.0, 2.0)
round(theta.m <- matrix(theta.s, nrow=sg_n, ncol=sg_k), 3)

#�x���k�[�C���z����ϑ��s��𔭐�������
X <- matrix(0, nrow=N, ncol=K)
for(i in 1:N){
  print(i)
  for(j in 1:K){
  X[i, j] <- rbinom(1, 1, theta.m[zno1[i], zno2[j]])
  }
}


####���Ӊ��M�u�X�T���v�����O�ŋ��N���X�^�����O�𐄒�####
##MCMC�A���S���Y���̐ݒ�
R <- 20000
keep <- 4
sbeta <- 1.5

##�n�C�p�[�p�����[�^�̐ݒ�
#�f�B�N�������z�̃p�����[�^
alpha1 <- rep(5, sg_n)
alpha2 <- rep(5, sg_k)

#�x�[�^���z�̃p�����[�^
alpha_b <- 1.0
theta_b <- 1.5

####�����������z���f���Ő��ݕϐ��̊����̏����l��ݒ�####
##���ݕϐ��̊����̏����l��ݒ�
##�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z���v�Z���邽�߂̊֐�
LLobz <- function(theta, r, Y, ID, n, k, freq, v){
  #�ޓx�Ƒΐ��ޓx���v�Z
  LLind <- matrix(0, nrow=n, ncol=k)
  for(i in 1:k){
    Li <- apply(cbind(Y, freq), 1, function(x) dmultinom(x[1:v], x[v+1], theta[i, ]))
    LLind[, i] <- Li 
  }
  
  #�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̌v�Z
  #������
  R <- matrix(r, nrow=n, ncol=k, byrow=T)
  
  #�l�ʂ̐��݊m���̌v�Z
  LLd <- matrix(0, nrow=n, ncol=k)
  LLl <- log(LLind)
  for(i in 1:n){
    if(NROW(ID[ID[, 2]==i, ]==1)) { 
      LLd[i, ] <- LLl[i, ]  
      } else {
      LLd[i, ] <- apply(LLl[ID[, 2]==i, ], 2, sum)
    }
  }

  #��������h�����߂ɑΐ��ޓx��-744�ȉ��̏ꍇ�͑ΐ��ޓx�𐓏グ����
  LL.min <- apply(LLd, 1, min)
  index.loss <- subset(1:nrow(LLd), (LL.min + 743) < 0)
  lplus <- -matrix((LL.min[index.loss] + 743), nrow=length(index.loss), ncol=k)
  LLd[index.loss, ] <- LLd[index.loss, ] + lplus
  
  #���݊m��z�̌v�Z
  LLho <- R * exp(LLd)
  z <- LLho / matrix(rowSums(LLho), nrow=n, ncol=k)

  #�ϑ��f�[�^�̑ΐ��ޓx���v�Z
  LLosum <- sum(log(apply(matrix(r, nrow=n, ncol=k, byrow=T) * exp(LLd), 1, sum)))
  rval <- list(LLobz=LLosum, z=z, LL=LLd)
  return(rval)
}


####EM�A���S���Y���Ń��[�U�[�̐��ݕϐ��̏����l�����肷��####
##EM�A���S���Y���̐ݒ�
#�X�V�X�e�[�^�X
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 1
iter <- 0
ID <- cbind(1:N, 1:N)

##�����l�̐ݒ�
alpha <- rep(5, K)   #���O���z�̃p�����[�^
r <- c(0.25, 0.25, 0.25, 0.25)   #�������̏����l

#�p�����[�^�̏����l
theta.f <- matrix(0, nrow=sg_n, ncol=K)
for(i in 1:sg_n){
  minmax <- colSums(X)
  pf <- runif(K, min(minmax), max(minmax))
  theta.f[i, ] <- pf/sum(pf)
}

#�ΐ��ޓx�̏�����
L <- LLobz(theta=theta.f, r=r, Y=X, ID=ID, n=N, k=sg_n, freq=rowSums(X), v=K)
LL1 <- L$LLob
z <- L$z

##EM�A���S���Y��
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  #E�X�e�b�v�̌v�Z
  z <- L$z   #���ݕϐ�z�̏o��
  zpt <- matrix(0, nrow=N, ncol=sg_n)
  for(i in 1:N){
    zpt[ID[, 2]==i, ] <- matrix(z[i, ], nrow=length(ID[ID[, 2]==i, 2]), ncol=sg_n, byrow=T)
  }
  
  #M�X�e�b�v�̌v�Z�ƍœK��
  #theta�̐���
  theta <- matrix(0, nrow=sg_n, ncol=K)
  for(j in 1:sg_n){
    #���S�f�[�^�̑ΐ��ޓx����theta�̐���ʂ��v�Z
    theta.seg <- (colSums(zpt[, j]*X)) / (sum(zpt[, j]*X))
    theta[j, ] <- as.matrix(theta.seg)
  }
  #�������𐄒�
  r <- apply(z, 2, sum) / N
  
  #�ϑ��f�[�^�̑ΐ��ޓx���v�Z
  L <- LLobz(theta=theta, r=r, Y=X, ID=ID, n=N, k=sg_n, freq=rowSums(X), v=K)
  LL <- L$LLob   #�ϑ��f�[�^�̑ΐ��ޓx
  iter <- iter+1   
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

##�����Z�O�����g������
Z1 <- as.numeric(t(apply(z, 1, function(x) rmultinom(1, 1, x))) %*% 1:sg_n)


####EM�A���S���Y���ŃA�C�e���̐��ݕϐ��̏����l�����肷��####
##EM�A���S���Y��������ɓ����܂ŃT���v�����O���J��Ԃ�
for(rp in 1:1000){
  #�X�V�X�e�[�^�X
  dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
  tol <- 1
  iter <- 0
  ID <- cbind(1:K, 1:K)
  
  ##�����l�̐ݒ�
  alpha <- rep(5, K)   #���O���z�̃p�����[�^
  r <- c(0.25, 0.25, 0.25)   #�������̏����l
  
  ##�A�C�e������ɂ���ƃ��[�U�[��(�ϐ�)�������̂Ń��[�U�[�������_���T���v�����O����
  N_samp <- 200
  index.u <- sample(1:nrow(X), N_samp)
  X_item <- t(X[index.u, ])
  
  #�p�����[�^�̏����l
  theta.f <- matrix(0, nrow=sg_k, ncol=N_samp)
  for(i in 1:sg_k){
    minmax <- colSums(X_item)
    pf <- runif(N_samp, min(minmax), max(minmax))
    theta.f[i, ] <- pf/sum(pf)
  }
  
  #�ΐ��ޓx�̏�����
  L <- LLobz(theta=theta.f, r=r, Y=X_item, ID=ID, n=K, k=sg_k, freq=rowSums(X_item), v=N_samp)
  LL1 <- L$LLob
  z <- L$z
  
  
  ##EM�A���S���Y��
  LL.vec <- c()
  while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
    #E�X�e�b�v�̌v�Z
    z <- L$z   #���ݕϐ�z�̏o��
    zpt <- matrix(0, nrow=N_samp, ncol=sg_k)
    for(i in 1:K){
      zpt[ID[, 2]==i, ] <- matrix(z[i, ], nrow=length(ID[ID[, 2]==i, 2]), ncol=sg_k, byrow=T)
    }
    
    #M�X�e�b�v�̌v�Z�ƍœK��
    #theta�̐���
    theta <- matrix(0, nrow=sg_k, ncol=N_samp)
    for(j in 1:sg_k){
      #���S�f�[�^�̑ΐ��ޓx����theta�̐���ʂ��v�Z
      theta.seg <- (colSums(zpt[, j]*X_item)) / (sum(zpt[, j]*X_item))
      theta[j, ] <- as.matrix(theta.seg)
    }
    #�������𐄒�
    r <- apply(z, 2, sum) / K
    
    #�ϑ��f�[�^�̑ΐ��ޓx���v�Z
    L <- LLobz(theta=theta, r=r, Y=X_item, ID=ID, n=K, k=sg_k, freq=rowSums(X_item), v=N_samp)
    LL <- L$LLob   #�ϑ��f�[�^�̑ΐ��ޓx
    iter <- iter+1   
    dl <- LL-LL1
    LL1 <- LL
    LL.vec <- c(LL.vec, LL)
    print(LL)
    if(is.nan(LL)==TRUE) break
  }
  print(rp)
  LL.vec <- LL.vec[is.nan(LL.vec)==FALSE]
  if(max(LL.vec)==LL.vec[length(LL.vec)]) break
}

##�����Z�O�����g������
Z2 <- as.numeric(t(apply(z, 1, function(x) rmultinom(1, 1, x))) %*% 1:sg_k)

##���v�ʂ̏����l���v�Z
#���[�U�[�Z�O�����g�̓��v��
M1k <- as.numeric(table(Z1))

Mkl <- matrix(0, nrow=sg_n, sg_k)
Xkl <- array(0, dim=c(N, K, sg_n*sg_k))
Zkl <- array(0, dim=c(N, K, sg_n*sg_k))
for(i in 1:sg_n){
  for(j in 1:sg_k){
    r <- (i-1) + j
    Zkl[, , r] <- as.matrix(ifelse(Z1==i, 1, 0), nrow=N, ncol=1) %*% ifelse(Z2==j, 1, 0)
    Xkl[, , r] <- X * Zkl[, , r]
    Mkl[i, j] <- sum(Xkl[, , r])
  }
}

####���Ӊ��M�u�X�T���v�����O�Ő��ݕϐ�z���T���v�����O####
##���[�U�[�̐��ݕϐ����T���v�����O
