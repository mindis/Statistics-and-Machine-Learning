#####�x�C�W�A�����݈��q���f��#####
library(MASS)
library(lda)
library(RMeCab)
library(bayesm)
library(extraDistr)
library(matrixStats)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(654978)

####�f�[�^�̔���####
#�f�[�^��ݒ�
hh <- 2000   #�T���v����
colums <- 100   #�ϐ���
k <- 10   #���ݕϐ���

#�o�C�i���s��𔭐�
Z0 <- matrix(0, nrow=hh, ncol=k)
for(j in 1:k){
  p <- runif(1, 0.25, 0.5)
  Z0[, j] <- rbinom(hh, 1, p)
}
r0 <- colMeans(Z0)   #������

#���q�s��𔭐�
X0 <- matrix(rnorm(hh*k, 0, 1), nrow=k, ncol=colums)

#�����ϐ��𔭐�
Cov0 <- 0.5
y <- Z0 %*% X0 + matrix(rnorm(hh*colums, 0, Cov0), nrow=hh, ncol=colums)


####�}���R�t�A�������e�J�����@�Ő��ݓ������f���𐄒�####
##�A���S���Y���̐ݒ�
R <- 10000
keep <- 4
iter <- 0

##���O���z�̐ݒ�
#���K���z�̕��U�̎��O���z�̃p�����[�^
Deltabar <- matrix(0, nrow=k, ncol=colums)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta <- 0.01 * diag(rep(1, k))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu <- colums   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(rep(1, colums)) #�t�E�B�V���[�g���z�̃p�����[�^
s0 <- 0.01
v0 <- 0.01

##�����l�̐ݒ�
#�����l�̌��𐶐�
iter <- 5000
Z_array <- array(0, dim=c(hh, k, iter))
storage.mode <- "integer"
sq_error <- rep(0, iter)

for(i in 1:iter){
  print(i)
  for(j in 1:k){
    p <- runif(1, 0.25, 0.6)
    Z_array[, j, i] <- rbinom(hh, 1, p)
  }
  X <- solve(t(Z_array[, , i]) %*% Z_array[, , i]) %*% t(Z_array[, , i]) %*% y
  sq_error[i] <- sum((y - Z_array[, , i] %*% X)^2)
}
min(sq_error)


#�x�X�g�ȏ����l��I��
bt <- which.min(sq_error)
Z <- Z_array[, , bt]
X <- solve(t(Z) %*% Z) %*% t(Z) %*% y
Cov <- 0.1
r <- colMeans(Z)
rm(Z_array)

##�T���v�����O���ʂ̊i�[�p�z��
Zi <- array(0, dim=c(hh, k, R/keep))
FACTOR <- array(0, dim=c(k, colums, R/keep))
SIGMA <- rep(0, R/keep)
storage.mode(Zi) <- "integer"
gc(); gc()


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  z1_rate <- matrix(0, nrow=hh, ncol=k)
  
  ##���ݕϐ�Z���T���v�����O
  for(j in 1:k){
    
    #�p�^�[���ʂɑΐ��ޓx���v�Z
    z0 <- z1 <- Z
    z0[, j] <- 0; z1[, j] <- 1
    LLi0 <- rowSums(dnorm(y, z0 %*% X, Cov, log=TRUE))
    LLi1 <- rowSums(dnorm(y, z1 %*% X, Cov, log=TRUE))
    
    #logsumexp�̖ޓx���v�Z
    LLi <- cbind(LLi0, LLi1)
    
    LLi_max <- matrix(apply(LLi, 1, max), nrow=hh, 2)
    r_matrix <- matrix(c(1-r[j], r[j]), nrow=hh, ncol=2, byrow=T)  

    #�����m���̃p�����[�^��ݒ�
    expl <- r_matrix * exp(LLi - LLi_max)
    expl_log <- log(expl)
    expl_max <- log(max(expl))
    z1_rate[, j] <- exp(expl_log[, 2] - (log(rowSums(exp(expl_log - expl_max))) + expl_max))   #�Z�O�����g�����m��
    
    #�x���k�[�C���z�����ݕϐ��𐶐�
    Z[, j] <- as.integer(z1_rate[, j] > runif(hh))
  }
 
  #���������X�V
  r <- colMeans(Z)
  
  ##���ϗʉ�A���f���ň��q�s��X���X�V
  out <- rmultireg(Y=y, X=Z, Bbar=Deltabar, A=ADelta, nu=nu, V=V)
  X <- out$B
  
  ##�t�K���}���z����W���΍����X�V
  er <- as.numeric(y - Z %*% X)
  s <- s0 + t(er) %*% er
  v <- v0 + hh*colums
  Cov <- sqrt(1/(rgamma(1, v/2, s/2)))   #�t�K���}���z����sigma^2���T���v�����O

  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    Zi[, , mkeep] <- Z
    FACTOR[, , mkeep] <- X
    SIGMA[mkeep] <- Cov
      
    #�T���v�����O���ʂ̕\��
    print(rp)
    print(round(rbind(r, r0), 3))
    print(round(c(Cov, Cov0), 3))
    print(cbind(Z[1:20, ], Z0[1:20, ]))
  }
}

####�T���v�����O���ʂ̊m�F�ƓK���x�̊m�F####
burnin <- 250   #�o�[���C������(1000�T���v���܂�)
RS <- R/keep 

##�T���v�����O���ʂ̉���
matplot(t(FACTOR[, 1, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(FACTOR[, 2, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(FACTOR[, 3, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
plot(1:RS, SIGMA, type="l", xlab="�T���v�����O��")

##����v��l
cbind(apply(Zi[, , burnin:RS], c(1, 2), mean), Z0)
round(cbind(t(apply(FACTOR[, , burnin:RS], c(1, 2), mean)), t(X0)), 2)
round(c(mean(SIGMA[burnin:RS]), Cov0), 3)


pnorm(100)
pnorm(-100)