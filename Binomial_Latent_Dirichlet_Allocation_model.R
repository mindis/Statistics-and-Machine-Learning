#####Binomial Latent Dirichlet Allocation model#####
options(warn=0)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
r <- 10   #�]���X�R�A��
k1 <- 7   #���[�U�[�g�s�b�N��
k2 <- 8   #�A�C�e���g�s�b�N��
K <- matrix(1:(k1*k2), nrow=k1, ncol=k2, byrow=T)   #�g�s�b�N�̔z��
hh <- 3000   #���r���A�[��
item <- 1000   #�A�C�e����

##ID�ƌ����x�N�g���̍쐬
#ID�����ݒ�
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)

#�����x�N�g�����쐬
for(rp in 1:100){
  m_vec <- rep(0, hh*item)
  for(i in 1:item){
    prob <- rbeta(1, 8.0, 50.0)
    m_vec[item_id0==i] <- rbinom(hh, 1, prob)
  }
  m_index <- which(m_vec==1)
  
  #���S��ID��ݒ�
  user_id <- user_id0[m_index]
  item_id <- item_id0[m_index]
  d <- length(user_id)   #�����r���[��
  
  #���ׂẴp�^�[��������������break
  if(length(unique(user_id))==hh & length(unique(item_id))==item) break
}

##�p�����[�^�̐ݒ�
#�f�B���N�����z�̎��O���z��ݒ�
alpha11 <- rep(0.2, k1)
alpha12 <- rep(0.2, k2)
alpha2 <- rep(0.15, r)

#�f�B���N�����z����p�����[�^�𐶐�
theta1 <- thetat1 <- extraDistr::rdirichlet(hh, alpha11)
theta2 <- thetat2 <- extraDistr::rdirichlet(item, alpha12)
phi <- phit <- extraDistr::rdirichlet(k1*k2, alpha2)

#���f���Ɋ�Â��f�[�^�𐶐�
y <- matrix(0, nrow=d, ncol=r)
y_vec <- rep(0, d)
Z1 <- matrix(0, nrow=d, ncol=k1)
Z2 <- matrix(0, nrow=d, ncol=k2)

for(i in 1:d){
  if(i%%10000==0){
    print(i)
  }
  #�f�[�^�C���f�b�N�X�𒊏o
  u_index <- user_id[i]
  i_index <- item_id[i]
  
  #�g�s�b�N�𐶐�
  z1 <- as.numeric(rmnom(1, 1, theta1[u_index, ]))
  z2 <- as.numeric(rmnom(1, 1, theta2[i_index, ]))
  
  #�g�s�b�N�Ɋ�Â��]���X�R�A�𐶐�
  k_index <- K[which.max(z1), which.max(z2)]
  y[i, ] <- rmnom(1, 1, phi[k_index, ])
  y_vec[i] <- which.max(y[i, ])
  
  #�f�[�^���i�[
  Z1[i, ] <- z1
  Z2[i, ] <- z2
}

####�}���R�t�A�������e�J�����@��Bi LDA�𐄒�####
##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, w, k){
  #���S�W�����v�Z
  Bur <- theta[w, ] * t(phi)[wd, ]   #�ޓx
  Br <- Bur / rowSums(Bur)   #���S��
  r <- colSums(Br) / sum(Br)   #������
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##�A���S���Y���̐ݒ�
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##���O���z�̐ݒ�
alpha1 <- 0.25
alpha2 <- 0.25
beta <- 0.25

##�p�����[�^�̐^�l
theta1 <- thetat1
theta2 <- thetat2
phi <- phit

##�p�����[�^�̏����l��ݒ�
#�g�s�b�N���z�̏����l
theta1 <- extraDistr::rdirichlet(hh, rep(1.0, k1))
theta2 <- extraDistr::rdirichlet(item, rep(1.0, k2))
phi <- extraDistr::rdirichlet(k1*k2, rep(1.0, r)) 

##�p�����[�^�̊i�[�p�z��
THETA1 <- array(0, dim=c(hh, k1, R/keep))
THETA2 <- array(0, dim=c(item, k2, R/keep))
PHI <- array(0, dim=c(k1*k2, r, R/keep))
SEG1 <- matrix(0, nrow=d, ncol=k1)
SEG2 <- matrix(0, nrow=d, ncol=k2)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"


##�f�[�^�ƃC���f�b�N�X�̐ݒ�
#�C���f�b�N�X�̐ݒ�
user_index <- user_vec <- list()
item_index <- item_vec <- list()
y_list <- y_ones <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
  user_vec[[i]] <- rep(1, length(user_index[[i]]))
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
  item_vec[[j]] <- rep(1, length(item_index[[j]]))
}
for(j in 1:r){
  y_list[[j]] <- which(y_vec==j)
  y_ones[[j]] <- rep(1, length(y_list[[j]]))
}

#�f�[�^�̐ݒ�
r_vec <- rep(1, r)
k1_vec <- rep(1, k2)
k2_vec <- rep(1, k1)
index_k1 <- rep(1:k1, rep(k2, k1))
index_k2 <- rep(1:k2, k1)

##�ΐ��ޓx�̊�l�ƍŗǒl
#��̑ΐ��ޓx
par <- colSums(y) / d
LLst <- sum(y %*% log(par))

#�ŗǂ̑ΐ��ޓx
phi11 <- t(phit)[y_vec, ] * thetat2[item_id, index_k2]
par_u0 <- matrix(0, nrow=d, ncol=k1)
for(j in 1:k1){
  par_u0[, j] <- phi11[, K[j, ]] %*% k1_vec
}
par_u1 <- thetat1[user_id, ] * par_u0   #���[�U�[�g�s�b�N�̊��Җޓx
LLbest <- sum(log(rowSums(par_u1)))


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##���[�U�[�g�s�b�N���T���v�����O
  #���[�U�[�g�s�b�N�̏����t���m��
  phi11 <- t(phi)[y_vec, ] * theta2[item_id, index_k2]
  par_u0 <- matrix(0, nrow=d, ncol=k1)
  for(j in 1:k1){
    par_u0[, j] <- phi11[, K[j, ]] %*% k1_vec
  }
  par_u1 <- theta1[user_id, ] * par_u0   #���[�U�[�g�s�b�N�̊��Җޓx
  
  #���ݕϐ��̊����m������g�s�b�N���T���v�����O
  z1_rate <- par_u1 / as.numeric(par_u1 %*% k2_vec)
  Zi1 <- rmnom(d, 1, z1_rate)
  Zi1_T <- t(Zi1)
  
  
  ##�A�C�e���g�s�b�N���T���v�����O
  #�A�C�e���g�s�b�N�̏����t���m��
  phi12 <- t(phi)[y_vec, ] * theta1[user_id, index_k1]
  par_i0 <- matrix(0, nrow=d, ncol=k2)
  for(j in 1:k2){
    par_i0[, j] <- phi12[, K[, j]] %*% k2_vec
  }
  par_i1 <- theta2[item_id, ] * par_i0   #���[�U�[�g�s�b�N�̊��Җޓx
  
  #���ݕϐ��̊����m������g�s�b�N���T���v�����O
  z2_rate <- par_i1 / as.numeric(par_i1 %*% k1_vec)
  Zi2 <- rmnom(d, 1, z2_rate)
  Zi2_T <- t(Zi2)
  
  #���[�U�[�ƃA�C�e���g�s�b�N�𓝍�
  Zi <- Zi1[, index_k1] * Zi2[, index_k2]
  Zi_T <- t(Zi)
  
  
  ##�p�����[�^���T���v�����O
  #���[�U�[�̃g�s�b�N���z���T���v�����O
  wusum0 <- matrix(0, nrow=d, ncol=k1)
  for(i in 1:hh){
    wusum0[i, ] <- Zi1_T[, user_index[[i]]] %*% user_vec[[i]]
  }
  wusum <- wusum0 + alpha1   #�f�B���N�����z�̃p�����[�^
  theta1 <- extraDistr::rdirichlet(hh, wusum)   #�f�B���N�����z����theta11���T���v�����O
  
  #�A�C�e���̃g�s�b�N���z���T���v�����O
  wisum0 <- matrix(0, nrow=item, ncol=k2)
  for(i in 1:item){
    wisum0[i, ] <- Zi2_T[, item_index[[i]]] %*% item_vec[[i]]
  }
  wisum <- wisum0 + alpha2   #�f�B���N�����z�̃p�����[�^
  theta2 <- extraDistr::rdirichlet(item, wisum)   #�f�B���N�����z����theta11���T���v�����O
  
  
  ##�]���X�R�A���z���T���v�����O
  vsum0 <- matrix(0, nrow=k1*k2, ncol=r)
  for(j in 1:r){
    vsum0[, j] <- Zi_T[, y_list[[j]]] %*% y_ones[[j]]
  }
  vsum <- vsum0 + beta   #�f�B���N�����z�̃p�����[�^
  phi <- extraDistr::rdirichlet(k1*k2, vsum)   #�f�B���N�����z����eta���T���v�����O
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA1[, , mkeep] <- theta1
    THETA2[, , mkeep] <- theta2
    PHI[, , mkeep] <- phi
  }  
  
  #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
  if(rp%%keep==0 & rp >= burnin){
    SEG1 <- SEG1 + Zi1
    SEG2 <- SEG2 + Zi2
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    LL <- sum(log(rowSums(par_u1)))
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(c(LL, LLbest, LLst))
    print(round(cbind(theta1[1:10, ], thetat1[1:10, ]), 3))
  }
}
