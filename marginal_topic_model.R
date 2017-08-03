#####���Ӊ��g�s�b�N���f��#####
library(MASS)
library(lda)
library(RMeCab)
library(gtools)
library(reshape2)
library(plyr)
library(ggplot2)


####�f�[�^�̔���####
#set.seed(423943)
#�f�[�^�̐ݒ�
k <- 5   #�g�s�b�N��
d <- 200   #������
v <- 100   #��b��
w <- 200   #1����������̒P�ꐔ 

#�p�����[�^�̐ݒ�
alpha0 <- runif(k, 0.1, 1.25)   #�����̃f�B���N�����O���z�̃p�����[�^
alpha1 <- rep(0.5, v)   #�P��̃f�B���N�����O���z�̃p�����[�^

#�f�B���N�������̔���
theta0 <- rdirichlet(d, alpha0)   #�����̃g�s�b�N���z���f�B���N���������甭��
phi0 <- rdirichlet(k, alpha1)   #�P��̃g�s�b�N���z���f�B���N���������甭��

#�������z�̗�������f�[�^�𔭐�
WX <- matrix(0, d, v)
ZS <- list()

for(i in 1:d){
  z <- t(rmultinom(w, 1, theta0[i, ]))   #�����̃g�s�b�N���z�𔭐�
  
  zn <- z %*% c(1:k)   #0,1�𐔒l�ɒu��������
  zdn <- cbind(zn, z)   #apply�֐��Ŏg����悤�ɍs��ɂ��Ă���
  
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi0[x[1], ])))   #�����̃g�s�b�N����P��𐶐�
  wdn <- colSums(wn)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  WX[i, ] <- wdn  
  ZS[[i]] <- zdn[, 1]
  print(i)
}
barplot(colSums(z), names.arg=c("seg1", "seg2", "seg3", "seg4", "seg5"))
round(colSums(WX)/sum(WX), 3)   #�P��̏o���p�x


####���Ӊ��M�u�X�T���v�����O�Ńg�s�b�N���f���𐄒�####