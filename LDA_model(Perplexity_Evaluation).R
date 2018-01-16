#####Latent Dirichlet Allocation���f��(Perplexity)#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(21437)

####�f�[�^�̔���####
#set.seed(423943)
#�f�[�^�̐ݒ�
k <- 15   #�g�s�b�N��
d <- 3000   #������
v <- 500   #��b��
w <- rpois(d, rgamma(d, 70, 0.40))   #1����������̒P�ꐔ
f <- sum(w)


#�p�����[�^�̐ݒ�
alpha0 <- rep(0.25, k)   #�����̃f�B���N�����O���z�̃p�����[�^
alpha1 <- rep(0.1, v)   #�P��̃f�B���N�����O���z�̃p�����[�^

#�f�B���N�������̔���
thetat <- theta <- rdirichlet(d, alpha0)   #�����̃g�s�b�N���z���f�B���N���������甭��
phit <- phi <- rdirichlet(k, alpha1)   #�P��̃g�s�b�N���z���f�B���N���������甭��

#�������z�̗�������f�[�^�𔭐�
WX <- matrix(0, nrow=d, ncol=v)
Z <- list()

for(i in 1:d){
  z <- rmnom(w[i], 1, theta[i, ])   #�����̃g�s�b�N���z�𔭐�
  z_vec <- z %*% c(1:k)   #�g�s�b�N�������x�N�g����
  
  wn <- rmnom(w[i], 1, phi[z_vec, ])   #�����̃g�s�b�N�J���P��𐶐�
  wdn <- colSums(wn)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  WX[i, ] <- wdn
  Z[[i]] <- z
  print(i)
}

####�g�s�b�N���f������̂��߂̃f�[�^�Ɗ֐��̏���####
##���ꂼ��̕������̒P��̏o������ѕ⏕���̏o�����x�N�g���ɕ��ׂ�
##�f�[�^����pID���쐬
n1 <- 2000   #�w�K�p�T���v����
n2 <- 1000   #���ؗp�T���v����
WX1 <- WX[1:n1, ]   #�w�K�p�f�[�^
w1 <- w[1:n1]
WX2 <- WX[(n1+1):d, ]   #���ؗp�f�[�^
w2 <- w[(n1+1):d]
ID_list <- list()
wd_list <- list()

#�������Ƃɕ���ID����ђP��ID���쐬
for(i in 1:nrow(WX1)){
  print(i)
  
  #�P���ID�x�N�g�����쐬
  ID_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
  num2 <- which(num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
}

#���X�g���x�N�g���ɕϊ�
ID_d <- unlist(ID_list)
wd <- unlist(wd_list)

##�C���f�b�N�X���쐬
doc_list <- list()
word_list <- list()
for(i in 1:length(unique(ID_d))) {doc_list[[i]] <- which(ID_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- which(wd==i)}
gc(); gc()


####�}���R�t�A�������e�J�����@�őΉ��g�s�b�N���f���𐄒�####
##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #���S�W���̊i�[�p
  for(j in 1:k){
    #���S�W�����v�Z
    Bi <- rep(theta[, j], w) * phi[j, wd]   #�ޓx
    Bur[, j] <- Bi   
  }
  
  Br <- Bur / rowSums(Bur)   #���S���̌v�Z
  r <- colSums(Br) / sum(Br)   #�������̌v�Z
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##�A���S���Y���̐ݒ�
R <- 2500   #�T���v�����O��
keep <- 2   #2���1��̊����ŃT���v�����O���ʂ��i�[
iter <- 0
burnin <- 500/keep
disp <- 10

##���O���z�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z
alpha01 <- 1.0
beta0 <- 0.5


##�p�����[�^�̏����l
theta <- rdirichlet(n1, rep(1, k))   #�����g�s�b�N�̃p�����[�^�̏����l
phi <- rdirichlet(k, colSums(WX1)/sum(WX1)*v)   #�P��g�s�b�N�̃p�����[�^�̏����l

##�p�����[�^�̊i�[�p�z��
f1 <- length(ID_d)
THETA <- array(0, dim=c(n1, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
SEG <- matrix(0, nrow=f1, ncol=k)
storage.mode(SEG) <- "integer"
gc(); gc()

##MCMC����p�z��
wsum0 <- matrix(0, nrow=n1, ncol=k)
vf0 <- matrix(0, nrow=k, ncol=v)


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�P��g�s�b�N���T���v�����O
  #�P�ꂲ�ƂɃg�s�b�N�̏o�������v�Z
  word_rate <- burden_fr(theta, phi, wd, w1, k)$Br
  
  #�������z����P��g�s�b�N���T���v�����O
  Zi <- rmnom(f1, 1, word_rate)   
  z_vec <- Zi %*% 1:k
  
  ##�P��g�s�b�N�̃p�����[�^���X�V
  #�f�B�N�������z����theta���T���v�����O
  for(i in 1:n1){
    wsum0[i, ] <- colSums(Zi[doc_list[[i]], ])
  }
  wsum <- wsum0 + alpha01 
  theta <- extraDistr::rdirichlet(n1, wsum)
  
  #�f�B�N�������z����phi���T���v�����O
  for(j in 1:v){
    vf0[, j] <- colSums(Zi[word_list[[j]], ])
  }
  vf <- vf0 + beta0
  phi <- extraDistr::rdirichlet(k, vf)
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    
    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(rp%%keep==0 & rp >= burnin){
      SEG <- SEG + Zi
    }
    
    #�T���v�����O���ʂ��m�F
    if(rp%%disp==0){
      print(rp)
      #print(round(cbind(theta[1:10, ], thetat[1:10, ]), 3))
      print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
    }
  }
}


####�T���v�����O���ʂ̉����Ɨv��####
burnin <- 1000/keep   #�o�[���C������
RS <- R/keep

##�T���v�����O���ʂ̉���
#�����̃g�s�b�N���z�̃T���v�����O����
matplot(t(THETA[1, , ]), type="l", ylab="�p�����[�^", main="����1�̃g�s�b�N���z�̃T���v�����O����")
matplot(t(THETA[2, , ]), type="l", ylab="�p�����[�^", main="����2�̃g�s�b�N���z�̃T���v�����O����")
matplot(t(THETA[3, , ]), type="l", ylab="�p�����[�^", main="����3�̃g�s�b�N���z�̃T���v�����O����")
matplot(t(THETA[4, , ]), type="l", ylab="�p�����[�^", main="����4�̃g�s�b�N���z�̃T���v�����O����")

#�P��̏o���m���̃T���v�����O����
matplot(t(PHI[1, 1:10, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N1�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[2, 11:20, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N2�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[3, 21:30, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N3�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[4, 31:40, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N4�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[5, 41:50, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N5�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[6, 51:60, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N6�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[7, 61:70, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N7�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[8, 71:80, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N8�̒P��̏o�����̃T���v�����O����")

##�T���v�����O���ʂ̗v�񐄒��
#�g�s�b�N���z�̎��㐄���
topic_mu <- apply(THETA[, , burnin:(R/keep)], c(1, 2), mean)   #�g�s�b�N���z�̎��㕽��
round(cbind(topic_mu, thetat[1:n1, ]), 3)
round(topic_sd <- apply(THETA[, , burnin:(R/keep)], c(1, 2), sd), 3)   #�g�s�b�N���z�̎���W���΍�

#�P��o���m���̎��㐄���
word_mu <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #�P��̏o�����̎��㕽��
round(rbind(word_mu, phit)[, 1:50], 3)


####�e�X�g�f�[�^��p����Perplexity��]��####
##�C���f�b�N�X�̍쐬
ID2_list <- list()
wd2_list <- list()

#�������Ƃɕ���ID����ђP��ID���쐬
for(i in 1:nrow(WX2)){
  print(i)
  
  #�P���ID�x�N�g�����쐬
  ID2_list[[i]] <- rep(i, w2[i])
  num1 <- (WX2[i, ] > 0) * (1:v)
  num2 <- which(num1 > 0)
  W1 <- WX2[i, (WX2[i, ] > 0)]
  number <- rep(num2, W1)
  wd2_list[[i]] <- number
}

#���X�g���x�N�g���ɕϊ�
ID2_d <- unlist(ID2_list)
wd2 <- unlist(wd2_list)

##�C���f�b�N�X���쐬
doc2_list <- list()
word2_list <- list()
for(i in 1:length(unique(ID2_d))) {doc2_list[[i]] <- which(ID2_d==i)}
for(i in 1:length(unique(wd2))) {word2_list[[i]] <- which(wd2==i)}
gc(); gc()

#�A���S���Y���̐ݒ�
R <- 2500   #�T���v�����O��
keep <- 2   #2���1��̊����ŃT���v�����O���ʂ��i�[
iter <- 0
burnin <- 500/keep
disp <- 10
f2 <- length(wd2)

#�����l��ݒ�
theta <- extraDistr::rdirichlet(n2, rep(1, k))

#�p�����[�^�̊i�[�p�z��
THETA2 <- array(0, dim=c(n2, k, R/keep))
SEG2 <- matrix(0, nrow=f2, k)

for(rp in 1:R){
  
  ##�P��g�s�b�N���T���v�����O
  #�P�ꂲ�ƂɃg�s�b�N�̏o�������v�Z
  word_rate <- burden_fr(theta, word_mu, wd2, w2, k)$Br
  
  #�������z����P��g�s�b�N���T���v�����O
  Zi <- rmnom(f2, 1, word_rate)   
  z_vec <- Zi %*% 1:k
  
  ##�P��g�s�b�N�̃p�����[�^���X�V
  #�f�B�N�������z����theta���T���v�����O
  wsum0 <- matrix(0, nrow=n2, ncol=k)
  for(i in 1:n2){
    wsum0[i, ] <- colSums(Zi[doc2_list[[i]], ])
  }
  wsum <- wsum0 + alpha01 
  theta <- extraDistr::rdirichlet(n2, wsum)

  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA2[, , mkeep] <- theta
    
    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(rp%%keep==0 & rp >= burnin){
      SEG2 <- SEG2 + Zi
    }
    
    #�T���v�����O���ʂ��m�F
    if(rp%%disp==0){
      print(rp)
      print(round(cbind(theta[1:10, ], thetat[(n1+1):(n1+10), ]), 3))
    }
  }
}

##Perplexity���v�Z
#���j�O�������f����Perplexity���v�Z
par <- colSums(WX1) / sum(WX1)
LL1 <- sum(WX2 %*% log(par))
Perx1 <- exp(-LL1 / f2)   #Perplexity

#�g�s�b�N���f����Perplexity���v�Z
theta_mu <- apply(THETA2[, , burnin:(R/keep)], c(1, 2), mean)   #theta�̎��㕽��
LL2 <- sum(log(rowSums(burden_fr(theta_mu, word_mu, wd2, w2, k)$Bur)))   #�ΐ��ޓx
Perx2 <- exp(-LL2 / f2)   #Perplexity
cbind(Perx1, Perx2)
