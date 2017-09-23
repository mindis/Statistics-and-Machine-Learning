#####�ϗʌ��ʃ��W�b�g���f���ɂ���ʉ��\�����_#####
library(MASS)
library(nlme)
library(glmm)
library(survival)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(57389)

####�f�[�^�̔���####
n1 <- 500   #���[�U�[��
n2 <- 150   #�ΏۃQ�[����
N <- n1*n2   #���T���v����
genre <- 6   #�Q�[���W��������

##ID�̐ݒ�
u.id <- rep(1:n1, rep(n2, n1))
g.id <- rep(1:n2, n1)
ID <- data.frame(no=1:N, u.id, g.id)

####�����ϐ��̔���####
##�Q�[���̓����ϐ��s����쐬
Disc <- rbinom(n2, 1, 0.6)
Genre0 <- t(rmultinom(n2, 1, runif(genre, 0.5, 2.0)))
Genre <- Genre0[, -which.min(colSums(Genre0))]
X1 <- cbind(1, Disc, Genre)   #�f�[�^�̌���

##�ϗʌ��ʂ̃f�U�C���s���ݒ�
Z1 <- matrix(0, nrow=N, ncol=n1*ncol(X1))
Z2 <- matrix(0, nrow=N, ncol=n2)
for(i in 1:n1){
  print(i)
  r <- ((i-1)*ncol(X1)+1):((i-1)*ncol(X1)+ncol(X1))
  Z1[ID$u.id==i, r] <- X1
}

for(i in 1:n2){
  print(i)
  Z2[ID$g.id==i, i] <- 1 
}

####�����ϐ��̔���####
##�p�����[�^�̐ݒ�
#���ύ\����ݒ�
theta0 <- runif(1, -1.0, -0.6)

#�l�n�D�̕ϗʌ��ʂ̐ݒ�
Cov01 <- diag(c(0.75, 0.35, runif(genre-1, 0.2, 0.4)))
alpha0 <- as.numeric(t(mvrnorm(n1, rep(0, ncol(X1)), Cov01)))

#�Q�[�����l�̕ϗʌ��ʂ̐ݒ�
Cov02 <- runif(1, 0.75, 0.8)
beta0 <- rnorm(n2, 0, Cov02)

##���W�b�g�Ɗm���̌v�Z
logit <- theta0 + Z1 %*% alpha0 + Z2 %*% beta0
Pr0 <- exp(logit)/(1+exp(logit))

##�����ϐ��̔���
Y <- rbinom(N, 1, Pr0)

#�����������ϐ��̉����ƏW�v
mean(Y); sum(Y); mean(Pr0); summary(Pr0)
hist(Pr0, col="grey", xlab="�m��", main="�w���m���̕��z")

