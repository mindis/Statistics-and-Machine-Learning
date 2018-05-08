#####�������q���̓��f��#####
options(warn=0)
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
library(ggplot2)

#set.seed(529867)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
k <- 10   #��ꐔ
hh <- 3000   #���[�U�[��
v <- 100   #���ڐ�

##�p�����[�^�̐ݒ�
beta1 <- 17.5
beta2 <- 25.0
mu_vec <- rep(0, k)
cov <- diag(k)
sigma <- 1

##�����ϐ��̐���
Z <- ZT <- matrix(rbinom(hh*k, 1, rbeta(hh*k, beta1, beta2)), nrow=hh, ncol=k)   #���ݕϐ��s��
X <- XT<- t(mvrnorm(v, mu_vec, cov))   #�����s��
Mu <- Z %*% X   #���ύ\���𐶐�
Y <- Mu + matrix(rnorm(hh*v, 0, sigma), nrow=hh, ncol=v)   #�����ϐ�


####�}���R�t�A�������e�J�����@�Ŗ������q���͂𐄒�####
##�A���S���Y���̐ݒ�
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10
g <- 5

##���O���z�̐ݒ�
sigma <- 1
mu_vec <- rep(0, k)
cov <- diag(k)
alpha <- 0.0001   #IBP�̎��O���z

##�����l�̐ݒ�
k0 <- 1
Z <- matrix(rbinom(hh, 1, 0.4), nrow=hh, ncol=)
X <- t(rnorm(v, 0, 1))
Mu <- Z %*% X

##�f�[�^�̐ݒ�
max_seg <- 20
vec <- 1:v
z_vec <- rep(1, hh)


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�x���k�[�C���z���������ݕϐ�z�𐶐�
  Z_new <- matrix(0, nrow=hh, ncol=k0+1)
  Z0 <- Z1 <- Z
  
  for(j in 1:k0){
    #�ΐ��ޓx��ݒ�
    Z1[, j] <- 1; Z0[, j] <- 0   #���ݕϐ���u������
    LLi0 <- as.numeric(dnorm(Y, Z0 %*% X, sigma, log=TRUE) %*% vec)
    LLi1 <- as.numeric(dnorm(Y, Z1 %*% X, sigma, log=TRUE) %*% vec)
    LLi <- cbind(LLi0, LLi1)
    Li <- exp(LLi - rowMaxs(LLi)) + 10^-100
    
    #���ݕϐ��̊����m������z�𐶐�
    gamma0 <- sum(Z[, j]) / hh   #���O���z
    gamma_par <- cbind((1-gamma0) * Li[, 1], gamma0 * Li[, 2])
    gamma <- gamma_par[, 2] / rowSums(gamma_par)   #���ݕϐ��̊����m��
    Z_new[, j] <- rbinom(hh, 1, gamma)   #�x���k�[�C���z������ݕϐ��𐶐�
  }
  Z <- Z_new
  
  ##IBP����V�������ݕϐ�z�𐶐�
  #�ΐ��ޓx��ݒ�
  X_new <- t(rnorm(v, 0, sigma))
  lambda <- alpha/hh
  LLi0 <- as.numeric(dnorm(Y, Z[, 1:k0] %*% X, sigma, log=TRUE) %*% vec)
  LLi1 <- as.numeric(dnorm(Y, Z[, 1:k0] %*% X + z_vec %*% X_new, sigma, log=TRUE) %*% vec) + dpois(k0+1, lambda, log=TRUE)
  LLi <- cbind(LLi0, LLi1)
  Li <- exp(LLi - rowMaxs(LLi)) 

  #���ݕϐ��̊����m������z�𐶐�
  gamma_par <- cbind((1-alpha/hh) * Li[, 1], alpha/hh * Li[, 2])
  gamma <- gamma_par[, 2] / rowSums(gamma_par)   #���ݕϐ��̊����m��
  Z[, (k0+1)] <- rbinom(hh, 1, gamma)   #�x���k�[�C���z������ݕϐ��𐶐�
  
  #�V�������ݕϐ�����萔�ȉ��Ȃ���ݕϐ����폜
  if(sum(Z[, (k0+1)]) <= g){
    Z <- Z[, -(k0+1)]
  }
  k0 <- ncol(Z)   #��ꐔ���X�V
  
  lambda <- alpha/hh
  dpois(10, 10-lambda, log=TRUE)
  
  
  ##�����s����X�V
  #���ϗʉ�A���f���̕��σx�N�g��
  ZZ <- t(Z) %*% Z + diag(k0)
  inv_ZZ <- solve(ZZ)
  beta_mu <- inv_ZZ %*% t(Z) %*% Y
  
  #���ϗʐ��K���z�������s����T���v�����O
  X <- matrix(0, nrow=k0, ncol=v)
  for(j in 1:v){
    X[, j] <- mvrnorm(1, beta_mu[, j], inv_ZZ)
  }
  
  LL <- sum(dnorm(Y, Z %*% X, sigma, log=TRUE))
  print(LL)
  print(k0)
  print(colSums(Z))
}


sum(dnorm(Y, ZT %*% XT, sigma, log=TRUE))


