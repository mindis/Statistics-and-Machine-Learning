#####�[���ߏ�񕉒l�s����q����#####
options(warn=2)
library(MASS)
library(NMF)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(data.table)
library(ggplot2)

#set.seed(319547)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 3000   #���[�U�[��
item <- 500   #�A�C�e����
k <- 10   #���ݕϐ���

##�[���ߏ�NMF�̉���Ɋ�Â��f�[�^�𐶐�
#�K���}���z���p�����[�^��ݒ�
alpha01 <- 0.25; beta01 <- 1.0
alpha02 <- 0.25; beta02 <- 0.8
W <- WT <- matrix(rgamma(hh*k, alpha01, beta01), nrow=hh, ncol=k)
H <- HT <- matrix(rgamma(item*k, alpha02, beta02), nrow=k, ncol=item)
WH0 <- W %*% H

#���ݍw���s��z�𐶐�
tau <- rbeta(hh, 3, 4.5)   #�x�[�^���z�̃p�����[�^
Z <- matrix(0, nrow=hh, ncol=item)
for(i in 1:hh){
  Z[i, ] <- rbinom(item, 1, tau[i])
}
mean(tau)

#�|�A�\�����z����w���s��𐶐�
WHT <- WH <- WH0 * Z
Data <- matrix(0, nrow=hh, ncol=item)
for(j in 1:item){
  Data[, j] <- rpois(hh, WH[, j])
}


####�}���R�t�A�������e�J�����@�Ń[���ߏ�NMF�𐄒�####
##�|�A�\�����z�̑ΐ��ޓx�֐�
pois <- function(Data, lambda, const){
  LLi <- Data*log(lambda) - lambda - const
  return(LLi)
}

##�A���S���Y���̐ݒ�
R <- 5000
keep <- 2
disp <- 10
iter <- 0

##���O���z�̐ݒ�
alpha1 <- 0.01; beta1 <- 0.01
alpha2 <- 0.01; beta2 <- 0.01

##�p�����[�^�̐^�l
W <- WT
H <- HT
r <- rowMeans(Z)

##�����l�̐ݒ�
W <- matrix(rgamma(hh*k, 0.1, 0.1), nrow=hh, ncol=k)
H <- matrix(rgamma(item*k, 0.1, 0.1), nrow=k, ncol=item)
r <- rowMeans(Data > 0)

##�p�����[�^�̊i�[�p�z��
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, item, R/keep))
Z_data <- matrix(0, nrow=hh, ncol=item)


##�f�[�^�̐ݒ�
z <- as.numeric(Data > 0)
const <- lfactorial(Data)   #�|�A�\�����z�̑ΐ��ޓx�̒萔
const_vec <- as.numeric(const)
data_vec <- as.numeric(Data)   #�f�[�^���x�N�g����
u_id <- rep(1:hh, item)
t_id <- rep(1:item, rep(hh, item))

#�C���f�b�N�X���쐬
index_zeros <- which(data_vec==0)   #���ݕϐ��̐���Ώۂ̃x�N�g��


##�ΐ��ޓx�̊�l
LLbest <- sum(dpois(Data, WT %*% HT, log=TRUE) * Z)
LLst <- sum(dpois(Data, mean(Data), log=TRUE))


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  ##���ݕϐ�Z���T���v�����O
  #���ݕϐ�z�̊����m���̃p�����[�^
  r_vec <- r[u_id][index_zeros]   #�������̃x�N�g��
  Li_zeros <- exp(-as.numeric(W %*% H)[index_zeros])   #�f�[�^���[���̎��̖ޓx
  z_rate <- r_vec*Li_zeros / (r_vec*Li_zeros + (1-r_vec))   #���ݕϐ��̊����m��
  
  #�񍀕��z������ݕϐ�z���T���v�����O
  z[index_zeros] <- rbinom(length(index_zeros), 1, z_rate)
  Zi <- matrix(z, nrow=hh, ncol=item)
  r <- rowMeans(Zi)   #���������X�V
  
  
  ##�⏕�ϐ�lambda���X�V
  lambda <- array(0, dim=c(hh, item, k))
  WH <-  W %*% H
  for(j in 1:k){
    lambda[, , j] <- (W[, j] %*% t(H[j, ]) * Zi / WH) 
  }

  ##�K���}���z���W���T���v�����O
  weights <- colMeans(Zi)
  for(j in 1:k){
    w1 <- alpha1 + rowSums(lambda[, , j] * Data)
    w2 <- beta1 + sum(H[j, ] * weights)
    W[, j] <- rgamma(hh, w1, w2)   
  }
  W <- W / matrix(colSums(W), nrow=hh, ncol=k, byrow=T) * hh/5   #�e��x�N�g���𐳋K��
  
  ##�⏕�ϐ�lambda���X�V
  lambda <- array(0, dim=c(hh, item, k))
  WH <-  W %*% H
  for(j in 1:k){
    lambda[, , j] <- (W[, j] %*% t(H[j, ]) * Zi / WH)
  }
  
  ##�K���}���z���H���T���v�����O
  weights <- rowMeans(Zi)
  for(j in 1:k){
    h1 <- alpha2 + colSums(lambda[, , j] * Data)
    h2 <- beta2 + sum(W[, j] * weights)
    H[j, ] <- rgamma(item, h1, h2)  
  }
  
  ##�T���v�����O���ʂ̕ۑ��ƕ\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    W_array[, , mkeep] <- W
    H_array[, , mkeep ] <- H
    if(rp > 1000){
      Z_data <- Z_data + Zi
    }
    
    if(rp%%disp==0){
      print(rp)
      print(c(sum(dpois(Data, W %*% H, log=T) * Zi), LLbest, LLst))
      print(round(c(mean(Z), mean(Zi)), 3))
      print(round(cbind(W[1:10, 1:k], WT[1:10, 1:k]), 3))
      print(round(cbind(H[1:k, 1:7], HT[, 1:7]), 3))
    }
  }
}
