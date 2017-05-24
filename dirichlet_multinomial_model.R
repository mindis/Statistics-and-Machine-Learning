#####�f�B�N�����������z���f��#####
library(MASS)
library(gtools)
library(reshape2)
library(plyr)
library(ggplot2)

####�f�[�^�̔���####
#set.seed(3479)
#�f�[�^�̐ݒ�
N1 <- 1500   #�N���X1�̃T���v����
N2 <- 1500   #�N���X2�̃T���v����
N <- N1 + N2   #���T���v����
col <- 50   #�ϐ���
n <- round(rpois(N, 200), 0)   #�T���v��������̕p�x

##�f�B�N�����������z����f�[�^�𔭐�������
##�N���X1�̃f�[�^����
#�f�B�N�������z����̗�������
alpha1 <- runif(col, 0.1, 1.5)   #�p�����[�^
theta1 <- rdirichlet(N1, alpha1)   #�f�B���N�������𔭐�
round(theta1[, 1:15], 3)   #�f�[�^���m�F

#�������z����̗�������
Z1 <- cbind(n[1:N1], theta1)
W1 <- t(apply(Z1, 1, function(x) rmultinom(1, x[1], x[-1])))
W1[, 1:20]

##�N���X2�̃f�[�^����
#�f�B�N�������z����̗�������
alpha2 <- runif(col, 0.5, 2.0)   #�p�����[�^
theta2 <- rdirichlet(N2, alpha2)   #�f�B���N�������𔭐�
round(theta2[, 1:15], 3)   #�f�[�^���m�F

#�������z����̗�������
Z2 <- cbind(n[(N1+1):N], theta2)
W2 <- t(apply(Z2, 1, function(x) rmultinom(1, x[1], x[-1])))
W2[, 1:20]