#####�g�s�b�N���f���ɂ��m���I�u���b�N���f��#####
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

#set.seed(2578)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
d <- 10000   #���[�U�[��
k1 <- 10   #���[�U�[�̃Z�O�����g��
k2 <- 12   #�g�s�b�N��
v <- 1000   #�A�C�e����
w <- rpois(d, rgamma(d, 15, 0.5))   #�w���A�C�e����
w[w < 5] <- ceiling(runif(sum(w < 5), 5, 15))
f <- sum(w)

##�p�����[�^�̐ݒ�
#�f�B�N�������z�̃p�����[�^��ݒ�
alpha01 <- rep(0.25, k2)   #�g�s�b�N���z�̃p�����[�^
alpha11 <- rep(0.1, v)   #�A�C�e���w���m���̃p�����[�^

#�f�B�N�������z�̃p�����[�^�𐶐�
omega <- omegat <- rep(1/k1, k1)   #������
theta <- thetat <- extraDistr::rdirichlet(k1, alpha01)   #�g�s�b�N���z�̃p�����[�^
phi <- phit <- extraDistr::rdirichlet(k2, alpha11)   #�P�ꕪ�z�̃p�����[�^


##�������z���A�C�e���w���s��𐶐�
WX <- matrix(0, nrow=d, ncol=v)
Z1_list <- list()
Z2_list <- list()

for(i in 1:d){
  #���[�U�[�̃Z�O�����g�𐶐�
  z1 <- rmnom(1, 1, omega)
  z1_vec <- as.numeric(z1 %*% 1:k1)
  
  #�����������[�U�[�Z�O�����g����g�s�b�N�𐶐�
  z2 <- rmnom(w[i], 1, theta[z1_vec, ])
  z2_vec <- as.numeric(z2 %*% 1:k2)
  
  #�g�s�b�N����A�C�e���w���𐶐�
  wd <- rmnom(w[i], 1, phi[z2_vec, ])
  WX[i, ] <- colSums(wd)
  
  #�p�����[�^���i�[
  Z1_list[[i]] <- z1
  Z2_list[[i]] <- z2
}


#���X�g�`����ϊ�
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)
v <- length(which(colSums(WX) > 0))
index <- which(colSums(WX) > 0)
WX <- WX[, index]
phit <- phit[, index]
const <- lfactorial(w) - rowSums(lfactorial(WX))   #�������z�̖��x�֐��̑ΐ��ޓx�̒萔

##�f�[�^����pID���쐬
ID_list <- list()
wd_list <- list()

#���l���Ƃɋ��lID����ђP��ID���쐬
for(i in 1:nrow(WX)){
  
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
for(i in 1:d) {doc_list[[i]] <- which(ID_d==i)}
for(i in 1:v) {word_list[[i]] <- which(wd==i)}
gc(); gc()


####�}���R�t�A�������e�J�����@�Ŋm���I�g�s�b�N�u���b�N���f���𐄒�####
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

##���ݕϐ�z���v�Z����֐�
LLobz <- function(WX, theta, r, const, d, k){
  
  #�������z�̑ΐ��ޓx
  log_theta <- log(t(theta))
  LLi <- const + WX %*% log_theta
  
  #logsumexp�̖ޓx
  LLi_max <- matrix(apply(LLi, 1, max), nrow=d, ncol=k)
  r_matrix <- matrix(r, nrow=d, ncol=k, byrow=T)
  
  #�����m���̃p�����[�^��ݒ�
  expl <- r_matrix * exp(LLi - LLi_max)
  expl_log <- log(expl)
  expl_max <- matrix(log(max(expl[1, ])), nrow=d, ncol=k)
  z <- exp(expl_log - (log(rowSums(exp(expl_log - expl_max))) + expl_max))   #�Z�O�����g�����m��
}

##�A���S���Y���̐ݒ�
R <- 10000
keep <- 2
burnin <- 1000/keep
iter <- 0
disp <- 10
LLt <- sum(log(rowSums(burden_fr(thetat[as.numeric(Z1 %*% 1:k1), ], phit, wd, w, k2)$Bur)))

##���O���z�̐ݒ�
alpha1 <- 1
beta1 <- 1

##�����l�̐ݒ�
r <- rep(1/k1, k1)   #�������̏����l
Zi1 <- rmnom(d, 1, rep(1/k1, k1))   #�Z�O�����g�����̏����l
theta <- extraDistr::rdirichlet(k1, rep(100, k2))   #�g�s�b�N���z�̏����l
phi <- extraDistr::rdirichlet(k2, colSums(WX)/sum(WX)*100)   #�P�ꕪ�z�̏����l

##�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(k1, k2, R/keep))
PHI <- array(0, dim=c(k2, v, R/keep))
SEG1 <- matrix(0, nrow=d, ncol=k1)
SEG2 <- matrix(0, nrow=f, ncol=k2)
Zi1 <- Z1


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�A�C�e�����ƂɃg�s�b�N���T���v�����O
  #�g�s�b�N���z�̊����m�����v�Z
  zi1 <- as.numeric(Zi1 %*% 1:k1)
  thetan <- theta[zi1, ]
  out <- burden_fr(thetan, phi, wd, w, k2)
  word_rate <- out$Br
  
  #�������z����g�s�b�N���T���v�����O
  Zi2 <- rmnom(f, 1, word_rate)
  z2_vec <- as.numeric(Zi2 %*% 1:k2)
  
  
  ##���[�U�[���ƂɃZ�O�����g���T���v�����O
  #�Z�O�����g�����m�����v�Z
  ZX <- matrix(0, nrow=d, ncol=k2)
  for(i in 1:d){
   ZX[i, ] <- colSums(Zi2[doc_list[[i]], , drop=FALSE])
  }
  LLi0 <- ZX %*% t(log(theta))
  LLho <- exp(LLi0 - apply(LLi0, 1, max))
  z1_rate <- r*LLho / rowSums(r*LLho)
  
  #�������z����Z�O�����g���T���v�����O
  Zi1 <- rmnom(d, 1, z1_rate)
  z1_vec <- as.numeric(Zi1 %*% 1:k1)
  
  #���������X�V
  #r <- matrix(colMeans(Zi1), nrow=d, ncol=k1, byrow=T)
  r <- matrix(1/k1, nrow=d, ncol=k1)
  
  ##�p�����[�^���X�V
  #�g�s�b�N���ztheta���T���v�����O
  wsum0 <- matrix(0, nrow=k1, ncol=k2)
  for(j in 1:k1){
    wsum0[j, ] <- colSums(ZX[Zi1[, j]==1, , drop=FALSE])
  }
  wsum <- wsum0 + alpha1
  theta <- extraDistr::rdirichlet(k1, wsum)
  
  #�P�ꕪ�zphi���T���v�����O
  vf0 <- matrix(0, nrow=k2, ncol=v)
  for(j in 1:v){
    vf0[, j] <- colSums(Zi2[word_list[[j]], , drop=FALSE])
  }
  vf <- vf0 + beta1
  phi <- extraDistr::rdirichlet(k2, vf)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    
    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(rp%%keep==0 & rp >= burnin){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
    }
    
    #�T���v�����O���ʂ��m�F
    if(rp%%disp==0){
      print(rp)
      print(c(sum(log(rowSums(out$Bur))), LLt))
      print(round(rbind(r=colMeans(Zi1), omega), 3))
      #print(round(cbind(theta, thetat), 3))
      print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
    }
  }
}

matplot(t(THETA[1, , ]), type="l")
matplot(t(THETA[2, , ]), type="l")
matplot(t(THETA[3, , ]), type="l")
matplot(t(PHI[, 1, ]), type="l")
matplot(t(PHI[, 2, ]), type="l")
matplot(t(PHI[, 3, ]), type="l")