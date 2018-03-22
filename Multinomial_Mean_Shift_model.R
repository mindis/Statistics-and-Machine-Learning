#####Mean Shift Model#####
options(warn=2)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(SMC)
library(RcppSMC)
library(SMC)
library(KFAS)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(93441)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
d <- 1000   #������
v <- 300   #��b��
a <- rpois(d, rgamma(d, 15, 0.7))   #��؂蕶�͐�
f1 <- sum(a)   #�����͐�
w <- rtpois(f1, 8, 5, Inf)   #��؂蕶�͂��Ƃ̒P�ꐔ
f2 <- sum(w)   #����b��

##ID�̐ݒ�
w_id1 <- rep(1:d, a)


##�p�����[�^�̐ݒ�
#���O���z�̃p�����[�^
alpha01 <- rep(0.2, v)   #�f�B���N�����z�̃p�����[�^
alpha11 <- 20   #�x�[�^���z�̃p�����[�^
beta11 <- 100

#�p�����[�^�𐶐�
gamma <- rbeta(d, alpha11, beta11)   #�؊����m���̃p�����[�^

#���f���Ɋ�Â��f�[�^�𐶐�
WX <- matrix(0, nrow=f1, ncol=v)
storage.mode(WX) <- "integer"
Z1_list <- list()
Z2_list <- list()
theta1_list <- list()
theta2_list <- list()

for(i in 1:d){
  
  #�؊����ϐ��𐶐�
  z <- rbinom(a[i], 1, gamma[i])
  z[1] <- 1   #1���ڂ�z=1
  z_no <- cumsum(z)   #�؊��������ɕϊ�
  
  #�f�B���N�����z����p�����[�^�𐶐�
  theta <- extraDistr::rdirichlet(sum(z), alpha01)
  theta_all <- theta[z_no, ]
  
  #�������z����P��𐶐�
  words <- rmnom(a[i], w[w_id1==i], theta_all)
  
  #�f�[�^���i�[
  WX[w_id1==i, ] <- words
  Z1_list[[i]] <- z
  Z2_list[[i]] <- z_no
  theta1_list[[i]] <- theta
  theta2_list[[i]] <- theta_all
}

#���X�g��ϊ�
Z <- unlist(Z1_list)
Zn <- unlist(Z2_list)
theta <- do.call(rbind, theta1_list)
theta_all <- do.call(rbind, theta2_list)
rm(theta1_list); rm(theta2_list)


####���q�t�B���^��Multinomial Mean Shift model�𐄒�####


##���q�t�B���^�̐ݒ�
s <- 3000   #���q��
alpha1 <- 0.1
beta1 <- rep(1/v, v)

#�p�����[�^�𐶐�
par <- WX[1, ]
beta0 <- (par + alpha1) / sum(par + alpha1)

#�ޓx���X�V
Li0 <- as.numeric(WX[2, ] %*% t(log(beta0)))
Li1 <- as.numeric(WX[2, ] %*% t(log(beta1)))
Li <- cbind(Li0, Li1)

#�p�����[�^�����T���v�����O
resample_par <- exp(Li0 + Li1 - max(Li0 + Li1))
resample_par / sum(resample_par)



#���ݕϐ�z�̊����m�����X�V
z_par <- exp(Li - rowMaxs(Li))
z_rate <- z_par / rowSums(z_par)
z_rate
rbinom(s, z_rate[, 2])


WX[2, ] %*% log((WX[1, ] + alpha1) / sum(WX[1, ] + alpha1))
WX[2, ] %*% log(rep(alpha1, v) / (v * alpha1))

