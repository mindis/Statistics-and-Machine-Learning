######�����K�E�X�ߒ���A���f��#####
options(warn=2)
library(MASS)
library(kernlab)
library(GPfit)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(5698)

####�f�[�^�̔���####
#�f�[�^�̐ݒ�
seg <- 8   #������
d <- 5000   #�T���v����
k <- 100   #���͕ϐ���
w <- rpois(d, rgamma(d, 15.0, 0.5))
w[w < 3] <- ceiling(runif(sum(w < 3), 3, 10))

#�Z�O�����g�̐���
alpha01 <- as.numeric(extraDistr::rdirichlet(1, rep(20.0, seg)))
Z <- rmnom(d, 1, alpha01)
z <- as.numeric(Z %*% 1:seg)

index_z <- list()
for(j in 1:seg){
  index_z[[j]] <- which(z==j)
} 


##�f�[�^�̐���
#�������z������͕ϐ��𐶐�
for(j in 1:1000){
  alpha0 <- rep(0.15, k)
  theta <- thetat <- extraDistr::rdirichlet(seg, alpha0)   #�������z�̃p�����[�^
  Data <- rmnom(d, w, theta[z, ])
  if(min(colSums(Data)) >= 10) break
}
sparse_data <- as(Data, "CsparseMatrix")


#�J�[�l���֐��̐���
kern_data <- list()
n0 <- rep(0, seg)
for(j in 1:seg){
  index <- index_z[[j]]
  kern_data[[j]] <- Data[index, ] %*% t(Data[index, ])
  n0[j] <- nrow(kern_data[[j]])
}

#�K�E�X�ߒ�����Z�O�����g���Ƃɉ����ϐ��𐶐�
y <- rep(0, d)
sigma <- rep(0, seg)
for(j in 1:seg){
  K <- kern_data[[j]]
  n <- n0[j]
  sigma[j] <- runif(1, 0.5, 5.0)
  y[index_z[[j]]] <- mvrnorm(1, rep(0, n0[j]), K + diag(sigma[j], n))   #�K�E�X�ߒ����牞���ϐ��𐶐�
}


####MCMC-EM�A���S���Y���ō����K�E�X�ߒ���A���f���𐄒�####
##�A���S���Y���̐ݒ�
R <- 5000   #�T���v�����O��
keep <- 2   #2���1��̊����ŃT���v�����O���ʂ��i�[
disp <- 10
iter <- 0
burnin <- 1000/keep

#���O���z�̐ݒ�
alpha01 <- rep(10, seg)
alpha02 <- rep(1, k)
r <- matrix(1/seg, nrow=d, ncol=seg, byrow=T)

#�����l�̐ݒ�
theta <- extraDistr::rdirichlet(seg, rep(5, k))
tau1 <- rep(1, seg)
tau2 <- rep(0.5, seg)

#�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(seg, k, R/keep))
SEG <- matrix(0, nrow=d, ncol=seg)
TAU1 <- matrix(0, nrow=R/keep, ncol=seg)
TAU2 <- matrix(0, nrow=R/keep, ncol=seg)

#�f�[�^�̐ݒ�
const <- lfactorial(w) - rowSums(lfactorial(Data))


####�p�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�����������z���f�������͕ϐ��̊������T���v�����O
  #�������z�̑ΐ��ޓx�֐�
  LLi <- as.matrix(const + sparse_data %*% t(log(theta)))

  #���ݕϐ��̊���Z���T���v�����O
  LL_max <- rowMaxs(LLi)
  LH <- exp(LLi - LL_max)   #�ޓx�ɕϊ�
  z_rate <- r * LH / rowSums(r * LH)   #���ݕϐ��̊����m��
  Zi <- rmnom(d, 1, z_rate)   #���ݕϐ����T���v�����O
  z_vec <- as.numeric(Zi %*% 1:seg)
  
  ##�����������z�̃p�����[�^���X�V
  #�������̍X�V
  rsum <- colSums(Zi) + alpha01
  r <- matrix(extraDistr::rdirichlet(1, rsum), nrow=d, ncol=seg, byrow=T)

  #�������z�̃p�����[�^���X�V
  wsum0 <- t(Data) %*% Zi
  wsum <- t(wsum0 + alpha02)
  theta <- extraDistr::rdirichlet(seg, wsum)   #�f�B���N�����z����theta���T���v�����O
  
  
  ##�Z�O�����g�����Ɋ�Â��K�E�X�ߒ���A���f����EM�A���S���Y���Ő���
  index_zi <- list()
  n <- rep(0, seg)
  LLm <- rep(0, j)
  beta_list <- list()
  
  for(j in 1:seg){
    #���ݕϐ�z�̃C���f�b�N�X���쐬
    index_zi[[j]] <- which(Zi[, j]==1)
    index <- index_zi[[j]]
  
    #�f�[�^�ƃJ�[�l���֐���ݒ�
    y_vec <- y[index]
    data <- Data[index, ]
    K <- data %*% t(data)
    KK <- K %*% K
    n[j] <- length(index_zi[[j]])
    n0 <- n[j]
    tau_s1 <- tau1[j]
    tau_s2 <- tau2[j]
  
    #E�X�e�b�v�ŉ�A�W���𐄒�
    beta <- solve(KK + diag(tau_s1/tau_s2, n0)) %*% t(K) %*% y_vec
    delta <- diag(tau_s1, n0) + tau_s2*KK
    beta_list[[j]] <- beta
    
    #M�X�e�b�v�Ńn�C�p�[�p�����[�^�𐄒�
    tau1_inv <- (sum(abs(beta^2)) + sum(diag(solve(delta)))) / n0
    tau1[j] <- 1 / tau1_inv
    tau2_inv <- (sum((y_vec - K %*% beta)^2) + sum(diag(KK %*% solve(delta)))) / n0
    tau2[j] <- 1 / tau2_inv
    
    #���Ӗޓx�̍X�V
    diag_tau2 <- diag(tau2_inv, n0)
    tau1_KK <- tau1_inv * KK
    Lm <- -1/2*abs(diag_tau2 + tau1_KK) - 1/2*as.numeric((t(y_vec) %*% solve(diag_tau2 + tau1_KK) %*% y_vec))
    LLm[j] <- sum(Lm)
  }
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    TAU1[mkeep, ] <- tau1
    TAU2[mkeep, ] <- tau2
    
    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(rp%%keep==0 & rp >= burnin){
      SEG <- SEG + Zi
    }
    
    #�T���v�����O���ʂ��m�F
    if(rp%%disp==0){
      print(rp)
      print(LLm)
      print(round(rbind(r[1, ], colMeans(Z)), 3))
      print(round(cbind(theta[, 1:10], thetat[, 1:10]), 3))
    }
  }
}

####���茋�ʂ̊m�F�Ɨv��####
burnin <- 1000/keep
RS <- R/keep

##�T���v�����O���ʂ̉���
#�������z�̃p�����[�^�̉���
matplot(t(THETA[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA[3, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA[7, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

#�������p�����[�^�̉���
matplot(TAU1/TAU2, type="l", xlab="�T���v�����O��", ylab="�p�����[�^")


##���茋�ʂ̗v�񐄒��
#�������p�����[�^�̗v�񓝌v��
tau_mu <- colMeans(TAU1[burnin:RS, ])/colMeans(TAU2[burnin:RS, ])

##���ݕϐ��̊������Ƃɉ�A�W���𐄒�
#���ݕϐ��̊���
Zi <- SEG / rowSums(SEG)

#�d�ݕt���f�[�^��ݒ�
data_list <- list()
K_list <- list()
y_list <- list()
beta_list <- list()

for(j in 1:seg){
  data0 <- Zi[, j] * Data; y_vec0 <- Zi[, j] * y   
  data_list[[j]] <- data <- data0[rowSums(abs(data0)) > 0, ]
  y_list[[j]] <- y_vec <- y_vec0[abs(y_vec0) > 0]
  n <- nrow(data)
  
  #��A�W���𐄒�
  K_list[[j]] <- K <- data %*% t(data)   #�J�[�l���֐�
  KK <- K %*% K
  beta_list[[j]] <- solve(KK + diag(tau_mu, n)) %*% t(K) %*% y_vec
}

#�\�����ʂ𐄒�
matplot(cbind(K_list[[1]] %*% beta_list[[1]], y_list[[1]]), type="l", xlab="�T���v���ԍ�", ylab="�p�����[�^")
matplot(cbind(K_list[[3]] %*% beta_list[[3]], y_list[[3]]), type="l", xlab="�T���v���ԍ�", ylab="�p�����[�^")
matplot(cbind(K_list[[5]] %*% beta_list[[5]], y_list[[5]]), type="l", xlab="�T���v���ԍ�", ylab="�p�����[�^")
matplot(cbind(K_list[[7]] %*% beta_list[[7]], y_list[[7]]), type="l", xlab="�T���v���ԍ�", ylab="�p�����[�^")
