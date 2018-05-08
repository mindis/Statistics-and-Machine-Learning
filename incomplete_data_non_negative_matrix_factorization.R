#####�����̂���񕉒l�s����q����#####
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
hh <- 5000   #���[�U�[��
item <- 1500  #�J�e�S���[��
k <- 10   #���ݕϐ���

##ID�̐ݒ�
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)


##�񕉒l�s����q�����̉���ɏ]���f�[�^�𐶐�
#�K���}���z���p�����[�^��ݒ�
alpha01 <- 0.25; beta01 <- 1.0
alpha02 <- 0.15; beta02 <- 0.85
W0 <- matrix(rgamma(hh*k, alpha01, beta01), nrow=hh, ncol=k)
H0 <- matrix(rgamma(item*k, alpha02, beta02), nrow=k, ncol=item)
WH <- W0 %*% H0

#�����L���̃x�[�^���z�̃p�����[�^��ݒ�
beta1 <- rbeta(hh, 9.5, 10.0)   #���[�U-�w���m��
beta2 <- rbeta(item, 7.5, 8.0)   #�A�C�e���w���m��

#�|�A�\�����z���f�[�^�𐶐�
Data0 <- matrix(0, nrow=hh, ncol=item)
for(j in 1:item){
  Data0[, j] <- rpois(hh, WH[, j])
}

##����������w���f�[�^�𐶐�
#�����s��𐶐�
Z0 <- matrix(0, nrow=hh, ncol=item)
for(j in 1:item){
  deficit <- rbinom(hh, 1, beta1 * beta2[j])
  Z0[, j] <- deficit   #��������
}

#�����C���f�b�N�X
Z <- Z0 * Data0 > 0
z_vec <- as.numeric(t(Z))
index_z <- which(z_vec==1)
N <- length(index_z)

#�����̂���w���x�N�g��
user_id <- user_id0[index_z]
item_id <- item_id0[index_z]
y <- as.numeric(t(Data0))[index_z]


##�x�X�g�ȃp�����[�^�ɑ΂���ΐ��ޓx
LLc1 <- sum(dpois(Data0, W0 %*% H0, log=TRUE))   #���S�f�[�^�ɑ΂���ΐ��ޓx
LLc2 <- sum(as.numeric(t(dpois(Data0, W0 %*% H0, log=TRUE)))[index_z])   #�s���S�f�[�^�ɑ΂���ΐ��ޓx
sparse_data <- as(Data0, "CsparseMatrix")   #�X�p�[�X�s��̍쐬


####�}���R�t�A�������e�J�����@��NMF�𐄒�####
##�A���S���Y���̐ݒ�
R <- 5000
keep <- 4
iter <- 0
disp <- 10

##���O���z�̐ݒ�
alpha1 <- 0.1; beta1 <- 1
alpha2 <- 0.1; beta2 <- 1

##�����l�̐ݒ�
W <- matrix(rgamma(hh*k, 0.1, 0.25), nrow=hh, ncol=k)
H <- matrix(rgamma(item*k, 0.1, 0.25), nrow=k, ncol=item)

##�T���v�����O���ʂ̕ۑ��p�z��
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, item, R/keep))
lambda <- array(0, dim=c(hh, item, k))

##���[�U�[����уA�C�e���̃C���f�b�N�X���쐬
user_list <- user_vec <- list()
item_list <- item_vec <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
  user_vec[[i]] <- rep(1, length(user_list[[i]]))
}
for(j in 1:item){
  item_list[[j]] <- which(item_id==j)
  item_vec[[j]] <- rep(1, length(item_list[[j]]))
}


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�⏕�ϐ�lambda���X�V
  lambda <- matrix(0, nrow=N, ncol=k)
  WH <- as.numeric(t(W %*% H))[index_z]
  for(j in 1:k){
    lambda[, j] <- as.numeric(t(W[, j] %*% t(H[j, ])))[index_z] / WH
  }
  
  ##�K���}���z��胆�[�U�[�����s��W���T���v�����O
  #���[�U�[���Ƃ̃K���}���z�̃p�����[�^��ݒ�
  W1 <- W2 <- matrix(0, nrow=hh, ncol=k)
  W1_T <- t(lambda * y)   #�v�f���Ƃ̊��Ғl
  for(i in 1:hh){
    W1[i, ] <- alpha1 + W1_T[, user_list[[i]], drop=FALSE] %*% user_vec[[i]]
    W2[i, ] <- beta1 + H[, item_id[user_list[[i]]], drop=FALSE] %*% user_vec[[i]]
  }
  
  #�K���}���z��胆�[�U�[�����s��W���T���v�����O
  W <- matrix(rgamma(hh*k, W1, W2), nrow=hh, ncol=k)
  W <- W / matrix(colSums(W), nrow=hh, ncol=k, byrow=T) * hh/5   #�e��x�N�g���𐳋K��
  
  
  ##�⏕�ϐ�lambda���X�V
  lambda <- matrix(0, nrow=N, ncol=k)
  WH <- as.numeric(t(W %*% H))[index_z]
  for(j in 1:k){
    lambda[, j] <- as.numeric(t(W[, j] %*% t(H[j, ])))[index_z] / WH
  }
  
  ##�K���}���z���A�C�e�������s��H���T���v�����O
  #�A�C�e�����Ƃ̃K���}���z�̃p�����[�^��ݒ�
  H1 <- H2 <- matrix(0, nrow=item, ncol=k)
  H1_T <- t(lambda * y)   #�v�f���Ƃ̊��Ғl
  for(i in 1:item){
    H1[i, ] <- alpha1 + H1_T[, item_list[[i]], drop=FALSE] %*% item_vec[[i]]
    H2[i, ] <- beta1 + t(W[user_id[item_list[[i]]], , drop=FALSE]) %*% item_vec[[i]]
  }
  #�K���}���z��胆�[�U�[�����s��W���T���v�����O
  H <- t(matrix(rgamma(item*k, H1, H2), nrow=item, ncol=k))
  
  
  ##�T���v�����O���ʂ̕ۑ��ƕ\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    W_array[, , mkeep] <- W[, 1:k]
    H_array[, , mkeep] <- H[1:k, ]
  }
  
  if(rp%%disp==0){
    print(rp)
    print(c(sum(dpois(y, as.numeric(t(W %*% H))[index_z], log=TRUE)), LLc2))
    print(round(W[1:10, ], 3))
  }
}

####�T���v�����O���ʂ̊m�F####
##�T���v�����O���ʂ��v���b�g
matplot(t(W_array[100, , ]), type="l")
matplot(t(H_array[, 1, ]), type="l")

##���㕽�ς��v�Z


##�e�X�g�f�[�^�ɑ΂���ΐ��ޓx
sum(dpois(as.numeric(t(Data0))[-index_z], as.numeric(t(W %*% H))[-index_z], log=TRUE))
sum(dpois(as.numeric(t(Data0))[-index_z], as.numeric(t(W0 %*% H0))[-index_z], log=TRUE))
round(cbind(as.numeric(t(Data0))[-index_z], as.numeric(t(W %*% H))[-index_z], as.numeric(t(W0 %*% H0))[-index_z]), 3)
round(cbind(y, as.numeric(t(W %*% H))[index_z], as.numeric(t(W0 %*% H0))[index_z]), 3)

