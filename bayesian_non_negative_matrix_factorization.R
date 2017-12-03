#####�x�C�W�A���񕉒l�s����q����#####
library(MASS)
library(matrixStats)
library(FAdist)
library(NMF)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)
 
####�f�[�^�̔���####
#�f�[�^�̐ݒ�
hh <- 3000   #���[�U�[��
category <- 200   #�J�e�S���[��
k <- 10   #���ݕϐ���

##�񕉒l�s����q�����̉���ɏ]���f�[�^�𐶐�
#�K���}���z���p�����[�^��ݒ�
alpha01 <- 0.2; beta01 <- 1.0
alpha02 <- 0.15; beta02 <- 0.8
W0 <- matrix(rgamma(hh*k, alpha01, beta01), nrow=hh, ncol=k)
H0 <- matrix(rgamma(category*k, alpha02, beta02), nrow=k, ncol=category)
WH <- W0 %*% H0

#�|�A�\�����z���f�[�^�𐶐�
Data <- matrix(0, nrow=hh, ncol=category)
for(j in 1:category){
  Data[, j] <- rpois(hh, WH[, j])
}
colSums(Data)
rowSums(Data)
LL <- sum(dpois(Data, W0 %*% H0, log=TRUE))
sparse_data <- as(Data, "CsparseMatrix")   #�X�p�[�X�s��̍쐬


####�}���R�t�A�������e�J�����@��NMF�𐄒�####
##�A���S���Y���̐ݒ�
R <- 10000
keep <- 4
iter <- 0

##���O���z�̐ݒ�
alpha1 <- 0.01; beta1 <- 0.01
alpha2 <- 0.01; beta2 <- 0.01

##�����l�̐ݒ�
W <- matrix(rgamma(hh*k, 0.1, 0.1), nrow=hh, ncol=k)
H <- matrix(rgamma(category*k, 0.1, 0.1), nrow=k, ncol=category)

##�T���v�����O���ʂ̕ۑ��p�z��
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, category, R/keep))
lambda <- array(0, dim=c(hh, category, k))


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�⏕�ϐ�lambda���X�V
  lambda <- array(0, dim=c(hh, category, k))
  WH <- W %*% H
  for(j in 1:k){
    lambda[, , j] <- W[, j] %*% t(H[j, ]) / WH
  }
  
  ##�K���}���z���W���T���v�����O
  for(j in 1:k){
    w1 <- alpha1 + rowSums(lambda[, , j] * Data)
    w2 <- beta1 + sum(H[j, ])
    W[, j] <- rgamma(hh, w1, w2)   
  }
  W <- W / matrix(colSums(W), nrow=hh, ncol=k, byrow=T) * hh/5   #�e��x�N�g���𐳋K��
  
  ##�⏕�ϐ�lambda���X�V
  lambda <- array(0, dim=c(hh, category, k))
  WH <- W %*% H
  for(j in 1:k){
    lambda[, , j] <- W[, j] %*% t(H[j, ]) / WH
  }

  ##�K���}���z���H���T���v�����O
  for(j in 1:k){
    h1 <- alpha2 + colSums(lambda[, , j] * Data)
    h2 <- beta2 + sum(W[, j])
    H[j, ] <- rgamma(category, h1, h2)  
  }

  ##�T���v�����O���ʂ̕ۑ��ƕ\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    W_array[, , mkeep] <- W[, 1:k]
    H_array[, , mkeep ] <- H[1:k, ]
    print(rp)
    print(c(sum(dpois(Data, W %*% H, log=T)), LL))
    print(round(cbind(W[1:10, 1:k], W0[1:10, 1:k]), 3))
    print(round(cbind(H[1:k, 1:7], H0[, 1:7]), 3))
  }
}

####�T���v�����O���ʂ̗v��Ɖ���####
burnin <- 2000/keep
RS <- R/keep

##�T���v�����O���ʂ̉���
#���x�N�g��W�̃p�����[�^�̉���
matplot(t(W_array[1:10, 1, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(W_array[1:10, 2, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(W_array[1:10, 3, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(W_array[1:10, 4, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(W_array[1:10, 5, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(W_array[11:20, 6, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(W_array[11:20, 7, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(W_array[11:20, 8, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(W_array[11:20, 9, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(W_array[11:20, 10, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")

#�����W��H�̃p�����[�^�̉���
matplot(t(H_array[1, 1:10, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(H_array[2, 1:10, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(H_array[3, 1:10, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(H_array[4, 1:10, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(H_array[5, 1:10, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(H_array[6, 11:20, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(H_array[7, 11:20, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(H_array[8, 11:20, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(H_array[9, 11:20, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(H_array[10, 11:20, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")


##�T���v�����O���ʂ̗v��
round(cbind(apply(W_array[, , burnin:RS], c(1, 2), mean), W0), 3)
round(apply(W_array[, , burnin:RS], c(1, 2), sd), 3)
round(cbind(t(apply(H_array[, , burnin:RS], c(1, 2), mean)), t(H0)), 3)
round(t(apply(H_array[, , burnin:RS], c(1, 2), sd)), 3)

##���㕽�ς���\���l���o��
W_mu <- apply(W_array[, , burnin:RS], c(1, 2), mean)
H_mu <- apply(H_array[, , burnin:RS], c(1, 2), mean)
WH_mu <- W_mu %*% H_mu

#�\�����ʂƃf�[�^���r
Predict <- matrix(0, nrow=hh, ncol=category*2)
index <- rep(0:1, category)
Predict[, index==0] <- round(WH_mu, 3)
Predict[, index==1] <- Data
