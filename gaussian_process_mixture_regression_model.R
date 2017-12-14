#####�K�E�X�ߒ�������A���f��#####
library(GPfit)
library(kernlab)
library(MASS)
library(Matrix)
library(matrixStats)
library(extraDistr)
library(reshape2)
library(dplyr)

#set.seed(85347)

####�f�[�^�̔���####
##�f�[�^�̔���
week <- 7   #����������
n <- week*104    #�T���v����
Data <- cbind(1, 1:n/10, matrix(diag(week), nrow=n, ncol=week, byrow=T)[, -1])   #�f�[�^�t���[��

##�����ϐ����O�����s��ɕϊ�
Kern1 <- kernelMatrix(rbfdot(sigma = 0.0001), Data[1:(n/2), ])
Kern2 <- kernelMatrix(rbfdot(sigma = 0.0001), Data[1:(n/2), ])

y1 <- mvrnorm(1, rep(35, n/2), Kern1)
y2 <- mvrnorm(1, rep(35, n/2), Kern2)
plot(1:n, c(y1, y2), type="l")

Kern1


