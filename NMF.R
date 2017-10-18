#####�I���̞B�����̂��߂̃}���`�N���X���ރ��f��####
library(MASS)
library(matrixStats)
library(flexmix)
library(glmnet)
library(mlogit)
library(nnet)
library(FAdist)
library(NMF)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

####�f�[�^�̔���####
hh <- 3000
seg <- 2

##�����ϐ��̔���
#���Ƃ̐����ϐ��s��𔭐�
topic <- 10   #�g�s�b�N��
k1 <- 300   #�����ϐ���
freq <- rpois(hh, 200)   #�|�A�\�����z����p�x�𔭐�

#�f�B���N�����z����o���m���𔭐�
#�p�����[�^�̐ݒ�
alpha0 <- runif(topic, 0.1, 1.0)   #�����̃f�B���N�����O���z�̃p�����[�^
theta0 <- rdirichlet(hh, alpha0)   #�����̃g�s�b�N���z���f�B���N���������甭��

alpha1 <- matrix(0, nrow=topic, ncol=k1)
phi0 <- matrix(0, nrow=topic, ncol=k1)
for(i in 1:topic){
  alpha1[i, ] <- rgamma(k1, 0.4, 0.1)   #�P��̃f�B���N�����O���z�̃p�����[�^
  phi0[i, ] <- rdirichlet(1, alpha1[i, ])   #�P��̃g�s�b�N���z���f�B���N���������甭��
}

#�������z�̗�������f�[�^�𔭐�
X0 <- matrix(0, hh, k1)
Topic <- list()

for(i in 1:hh){
  z <- t(rmultinom(freq[i], 1, theta0[i, ]))   #�����̃g�s�b�N���z�𔭐�
  
  zn <- z %*% c(1:topic)   #0,1�𐔒l�ɒu��������
  zdn <- cbind(zn, z)   #apply�֐��Ŏg����悤�ɍs��ɂ��Ă���
  
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi0[x[1], ])))   #�����̃g�s�b�N����P��𐶐�
  wdn <- colSums(wn)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  X0[i, ] <- wdn  
  Topic[[i]] <- zdn[, 1]
  print(i)
}

##�񕉒l�s����q�����ŕp�x�s������k
k <- 10
X0_trance <- t(X0)
res1 <- nmf(X0, k, "brunet")
res2 <- nmf(X0_trance, k, "brunet")

H <- res2@fit@H
W <- res2@fit@W

test <- round(cbind(t(solve(t(W) %*% W) %*% t(W) %*% X0_trance), t(H)), 3)


#�^�̃g�s�b�N�̏o���m���Ɛ��肳�ꂽ�g�s�b�N�m�����r
t_rate <- matrix(0, hh, topic) 
for(i in 1:hh){
  rate0 <- table(Topic[[i]])/freq[i]
  rate <- rep(0, topic)
  index <- subset(1:topic, 1:topic %in% names(rate0))
  rate[index] <- rate0
  t_rate[i, ] <- rate
}

#�œK�ȃg�s�b�N����10
opt.topic <- 10
Topic_rate <- round(cbind(t_rate, t(res2[[opt.topic-1]]@fit@H)/matrix(rowSums(t(res2[[opt.topic-1]]@fit@H)), nrow=hh, ncol=topic)), 3)

#�g�s�b�N�̏o���m��������ϐ��Ƃ���
X <- (t(res2[[opt.topic-1]]@fit@H)/matrix(rowSums(t(res2[[opt.topic-1]]@fit@H)), nrow=hh, ncol=opt.topic))[, -opt.topic]
XM <- cbind(1, X)