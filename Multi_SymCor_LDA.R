#####Multi SymCor LDA#####
library(MASS)
library(lda)
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

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
L <- 2   #�f�[�^�Z�b�g��
d <- 3000   #���͐�
k1 <- 10
k2 <- 15
v1 <- 1000; v2 <- 1000   #��b��
v11 <- 600; v12 <- v1-v11
v21 <- 650; v22 <- v2-v21
w1 <- rpois(d, rgamma(d, 60, 0.5))   
w2 <- rpois(d, rgamma(d, 45, 0.75))
f1 <- sum(w1); f2 <- sum(w2)

##ID�̐ݒ�
d_id1 <- rep(1:d, w1)
d_id2 <- rep(1:d, w2)
no_id1 <- no_id2 <- c()
for(i in 1:d){
  no_id1 <- c(no_id1, 1:w1[i])
  no_id2 <- c(no_id2, 1:w2[i])
}

##�p�����[�^�̐ݒ�
#�f�B���N�����z�̎��O���z
alpha01 <- rep(0.1, k1)  
alpha02 <- rep(0.1, k2)
alpha11 <- c(rep(0.15, v11), rep(0.001, v12))
alpha12 <- c(rep(0.001, v11), rep(0.1, v12))
alpha21 <- c(rep(0.1, v21), rep(0.001, v22))
alpha22 <- c(rep(0.001, v21), rep(2.5, v22))


##�f�[�^�̐���
for(rp in 1:1000){
  print(rp)
  
  #�p�����[�^�𐶐�
  beta1 <- betat1 <- rbeta(d, 25, 40)   #����1�ŕ����Ԃŋ��ʃg�s�b�N���ǂ����̊m��
  beta2 <- betat2 <- rbeta(d, 100, 40)   #����2�ŕ����Ԃŋ��ʃg�s�b�N���ǂ����̊m��
  theta1 <- thetat1 <- extraDistr::rdirichlet(d, alpha01)   
  theta2 <- thetat2 <- extraDistr::rdirichlet(d, alpha02) 
  phi1 <- phit1 <- extraDistr::rdirichlet(k1, alpha11)
  phi2 <- phit2 <- extraDistr::rdirichlet(k2, alpha12)
  omega1 <- omegat1 <- extraDistr::rdirichlet(k2, alpha21)
  omega2 <- omegat2 <- as.numeric(extraDistr::rdirichlet(1, alpha22))
  
  
  ##���f���̉���Ɋ�Â��f�[�^�𐶐�
  #�f�[�^�̊i�[�p�z��
  flag1_list <- list(); flag2_list <- list()
  WX1 <- matrix(0, nrow=d, ncol=v1)
  WX2 <- matrix(0, nrow=d, ncol=v2)
  word11_list <- word12_list <- word21_list <- word22_list <- list()
  Z11_list <- Z12_list <- Z21_list <- list()
  
  for(i in 1:d){
    if(i%%100==0){
      print(i)
    }
    #���������f�[�^�̊i�[�p�z��
    z11 <- matrix(0, nrow=w1[i], ncol=k1)
    z12 <- matrix(0, nrow=w1[i], ncol=k2)
    z21 <- matrix(0, nrow=w2[i], ncol=k2)
    word11 <- word12 <- matrix(0, nrow=w1[i], ncol=v1)
    word21 <- word22 <- matrix(0, nrow=w2[i], ncol=v2)
    
    #�P�ꂲ�ƂɃg�s�b�N�̋��ʐ����ǂ����𐶐�
    flag1 <- rbinom(w1[i], 1, beta1[i])
    flag2 <- rbinom(w2[i], 1, beta2[i])
    index_flag1 <- which(flag1==1)
    index_flag2 <- which(flag2==1)
    
    #�����������ʐ��Ɋ�Â��g�s�b�N�𐶐�
    z11[-index_flag1, ] <- rmnom(sum(1-flag1), 1, theta1[i, ])
    z12[index_flag1, ] <- rmnom(sum(flag1), 1, theta2[i, ])
    z21[index_flag2, ] <- rmnom(sum(flag2), 1, theta2[i, ])
    z11_vec <- as.numeric(z11 %*% 1:k1)
    z12_vec <- as.numeric(z12 %*% 1:k2)
    z21_vec <- as.numeric(z21 %*% 1:k2)
    
    #���������g�s�b�N����P��𐶐�
    word11[-index_flag1, ]<- rmnom(sum(1-flag1), 1, phi1[z11_vec[-index_flag1], ])
    word12[index_flag1, ]<- rmnom(sum(flag1), 1, phi2[z12_vec[index_flag1], ])
    word21[index_flag2, ]<- rmnom(sum(flag2), 1, omega1[z21_vec[index_flag2], ])
    word22[-index_flag2, ]<- rmnom(sum(1-flag2), 1, omega2)
    
    #���������f�[�^���i�[
    flag1_list[[i]] <- flag1; flag2_list[[i]] <- flag2
    Z11_list[[i]] <- z11; Z12_list[[i]] <- z12; Z21_list[[i]] <- z21
    word11_list[[i]] <- as.numeric(word11 %*% 1:v1)
    word12_list[[i]] <- as.numeric(word12 %*% 1:v1)
    word21_list[[i]] <- as.numeric(word21 %*% 1:v2)
    word22_list[[i]] <- as.numeric(word22 %*% 1:v2)
    WX1[i, ] <- colSums(word11) + colSums(word12)
    WX2[i, ] <- colSums(word21) + colSums(word22) 
  }
  if(min(colSums(WX1)) > 0 & min(colSums(WX2)) > 0) break
}

#���X�g��ϊ�
flag1 <- unlist(flag1_list); flag2 <- unlist(flag2_list)
Z11 <- do.call(rbind, Z11_list); Z12 <- do.call(rbind, Z12_list)
Z21 <- do.call(rbind, Z21_list)
wd11 <- unlist(word11_list); wd12 <- unlist(word12_list)
wd21 <- unlist(word21_list); wd22 <- unlist(word22_list)
wd1 <- wd11 + wd12
wd2 <- wd21 + wd22
rm(Z11_list); rm(Z12_list); rm(Z21_list)
rm(word11_list); rm(word12_list); rm(word21_list); rm(word22_list)
gc(); gc()


####�}���R�t�A�������e�J�����@��Multi CorSym LDA�𐄒�####
##�g�s�b�N���f���̖ޓx�ƕ��S�����`����֐�
burden_fr <- function(theta, phi, wd, w, k){
  #���S�W�����v�Z
  Bur <- theta[w, ] * t(phi)[wd, ]   #�ޓx
  Br <- Bur / rowSums(Bur)   #���S��
  r <- colSums(Br) / sum(Br)   #������
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##�A���S���Y���̐ݒ�
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##�C���f�b�N�X�̐ݒ�
doc1_list <- doc2_list <- list()
doc1_vec <- doc2_vec <- list()
wd1_list <- wd2_list <- list()
wd1_vec <- wd2_vec <- list()

for(i in 1:d){
  doc1_list[[i]] <- which(d_id1==i)
  doc2_list[[i]] <- which(d_id2==i)
  doc1_vec[[i]] <- rep(1, length(doc1_list[[i]]))
  doc2_vec[[i]] <- rep(1, length(doc2_list[[i]]))
}
for(j in 1:v1){
  wd1_list[[j]] <- which(wd1==j)
  wd1_vec[[j]] <- rep(1, length(wd1_list[[j]]))
}
for(j in 1:v2){
  wd2_list[[j]] <- which(wd2==j)
  wd2_vec[[j]] <- rep(1, length(wd2_list[[j]]))
}
topic_vec1 <- rep(1, k1)
topic_vec2 <- rep(1, k2) 
 

##���O���z�̐ݒ�
#�g�s�b�N���f���̎��O���z
alpha01 <- 0.1
alpha02 <- 0.01

#�x�[�^���z�̎��O���z
beta01 <- 1
beta02 <- 1

##�p�����[�^�̐^�l
beta1 <- betat1; beta2 <- betat2
theta1 <- thetat1; theta2 <- thetat2
phi1 <- phit1; phi2 <- phit2
omega1 <- omegat1; omega2 <- omegat2

##�����l�̐ݒ�
beta1 <- beta2 <- rep(0.5, d)
theta1 <- extraDistr::rdirichlet(d, rep(1.0, k1))
theta2 <- extraDistr::rdirichlet(d, rep(1.0, k2))
phi1 <- extraDistr::rdirichlet(k1, c(rep(1.5, v11), rep(0.1, v12)))
phi2 <- extraDistr::rdirichlet(k2, c(rep(0.1, v11), rep(1.5, v12)))
omega1 <- extraDistr::rdirichlet(k2, c(rep(1.5, v21), rep(0.1, v22)))
omega2 <- as.numeric(extraDistr::rdirichlet(1, c(rep(0.1, v21), rep(1.5, v22))))


##�p�����[�^�̊i�[�p�z��
BETA1 <- matrix(0, nrow=R/keep, ncol=d)
BETA2 <- matrix(0, nrow=R/keep, ncol=d)
THETA1 <- array(0, dim=c(d, k1, R/keep))
THETA2 <- array(0, dim=c(d, k2, R/keep))
PHI1 <- array(0, dim=c(k1, v1, R/keep))
PHI2 <- array(0, dim=c(k2, v1, R/keep))
OMEGA1 <- array(0, dim=c(k2, v2, R/keep))
OMEGA2 <- matrix(0, nrow=R/keep, ncol=v2)
FLAG1 <- rep(0, f1)
FLAG2 <- rep(0, f2)
SEG11 <- matrix(0, nrow=f1, ncol=k1)
SEG12 <- matrix(0, nrow=f1, ncol=k2)
SEG21 <- matrix(0, nrow=f2, ncol=k2)

##��̑ΐ��ޓx���v�Z
#���j�O�������f���̑ΐ��ޓx
LLst <- sum(WX1 %*% log(colSums(WX1)/f1)) + sum(WX2 %*% log(colSums(WX2)/f2))

#�^�l�ł̑ΐ��ޓx
phit11 <- (phit1+10^-100) / rowSums(phit1+10^-100); phit12 <- (phit2+10^-100) / rowSums(phit2+10^-100)
omegat11 <- (omegat1+10^-100) / rowSums(omegat1+10^-100); omegat12 <- (omegat2+10^-100) / sum(omegat2+10^-100)
LLbest1 <- sum(log((1-flag1)*rowSums(thetat1[d_id1, ]*t(phit11)[wd1, ]) + flag1*rowSums(thetat2[d_id1, ]*t(phit12)[wd1, ])))
LLbest2 <- sum(log(flag2*rowSums(thetat2[d_id2, ]*t(omegat11)[wd2, ]) + (1-flag2)*omegat12[wd2]))      
LLbest <- LLbest1 + LLbest2


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##����1�̃g�s�b�N�̐����ߒ����T���v�����O
  #���Җޓx���`
  phi1_T <- t(phi1); phi2_T <- t(phi2)
  Lho11 <- theta1[d_id1, ] * phi1_T[wd1, ]   
  Lho12 <- theta2[d_id1, ] * phi2_T[wd1, ] 
  topic_par11 <- Lho11 %*% topic_vec1   #�Ǝ������̊��Җޓx
  topic_par12 <- Lho12 %*% topic_vec2   #���ʕ����̊��Җޓx
  
  #�x���k�[�C���z���琶���ߒ����T���v�����O
  beta_par1 <- (1-beta1)[d_id1] * topic_par11
  beta_par2 <- beta1[d_id1] * topic_par12
  flag_rate1<- beta_par2 / (beta_par1+beta_par2)   #���ݕϐ��̊����m��
  flag_vec1 <- rbinom(f1, 1, flag_rate1)
  
  
  ##����2�̃g�s�b�N�̐����ߒ����T���v�����O
  #���Җޓx���`
  omega1_T <- t(omega1)
  Lho21 <- theta2[d_id2, ] * omega1_T[wd2, ]
  topic_par21 <- Lho21 %*% topic_vec2   #���ʕ����̊��Җޓx
  topic_par22 <- omega2[wd2]   #�Ǝ������̊��Җޓx
  
  #�x���k�[�C���z���琶���ߒ����T���v�����O
  beta_par1 <- beta2[d_id2] * topic_par21
  beta_par2 <- (1-beta2)[d_id2] * topic_par22
  flag_rate2 <- beta_par1 / (beta_par1+beta_par2)   #���ݕϐ��̊����m��
  flag_vec2 <- rbinom(f2, 1, flag_rate2)
  
  
  ##���������T���v�����O
  #�x�[�^���z�̃p�����[�^
  freq1 <- tapply(flag_vec1, d_id1, sum)
  freq2 <- tapply(flag_vec2, d_id2, sum)
  
  #�x�[�^���z����p�����[�^���T���v�����O
  beta1 <- rbeta(d, freq1+beta01, w1-freq1+beta02)
  beta2 <- rbeta(d, freq2+beta01, w2-freq2+beta02)
  
  
  ##����1�̃g�s�b�N���T���v�����O
  #�f�[�^�̐ݒ�
  Zi11 <- matrix(0, nrow=f1, ncol=k1)
  Zi12 <- matrix(0, nrow=f1, ncol=k2)
  index1 <- which(flag_vec1==0)
  
  #�������z����Ǝ��g�s�b�N���T���v�����O
  z_rate11 <- (Lho11 / as.numeric(topic_par11))[index1, ]   #�g�s�b�N�̊����m��
  Zi11[index1, ] <- rmnom(length(index1), 1, z_rate11)   #�g�s�b�N�̃T���v�����O
  zi11_vec <- as.numeric(Zi11 %*% 1:k1)
  Zi11_T <- t(Zi11)
  
  #�������z���狤�ʃg�s�b�N���T���v�����O
  z_rate12 <- (Lho12 / as.numeric(topic_par12))[-index1, ]   #�g�s�b�N�̊����m��
  Zi12[-index1, ] <- rmnom(f1-length(index1), 1, z_rate12)   #�g�s�b�N�̃T���v�����O
  zi12_vec <- as.numeric(Zi12 %*% 1:k2)
  Zi12_T <- t(Zi12)
  
  
  ##����2�̃g�s�b�N���T���v�����O
  #�f�[�^�̐ݒ�
  Zi21 <- matrix(0, nrow=f2, ncol=k2)
  index2 <- which(flag_vec2==1)
  
  #�������z����Ǝ��g�s�b�N���T���v�����O
  z_rate21 <- (Lho21 / as.numeric(topic_par21))[index2, ]   #�g�s�b�N�̊����m��
  Zi21[index2, ] <- rmnom(length(index2), 1, z_rate21)   #�g�s�b�N�̃T���v�����O
  zi21_vec <- as.numeric(Zi21 %*% 1:k2)
  Zi21_T <- t(Zi21)
  
  
  ##�g�s�b�N���z�̃p�����[�^���T���v�����O
  #����1�̓Ǝ��g�s�b�N���z���T���v�����O
  wsum0 <- matrix(0, nrow=d, ncol=k1)
  for(i in 1:d){
    wsum0[i, ] <- Zi11_T[, doc1_list[[i]]] %*% doc1_vec[[i]]
  }
  wsum <- wsum0 + alpha01   #�f�B���N�����z�̃p�����[�^
  theta1 <- extraDistr::rdirichlet(d, wsum)   #�f�B���N�����z����p�����[�^���T���v�����O
  
  #���ʃg�s�b�N���z���T���v�����O
  wsum0 <- matrix(0, nrow=d, ncol=k2)
  for(i in 1:d){
    wsum0[i, ] <- Zi12_T[, doc1_list[[i]]] %*% doc1_vec[[i]] + Zi21_T[, doc2_list[[i]]] %*% doc2_vec[[i]]
  }
  wsum <- wsum0 + alpha01   #�f�B���N�����z�̃p�����[�^
  theta2 <- extraDistr::rdirichlet(d, wsum)   #�f�B���N�����z����p�����[�^���T���v�����O
  
  
  ##�P�ꕪ�z�̃p�����[�^���T���v�����O
  #����1�̒P�ꕪ�z���T���v�����O
  vsum01 <- matrix(0, nrow=k1, ncol=v1)
  vsum02 <- matrix(0, nrow=k2, ncol=v1)
  for(j in 1:v1){
    vsum01[, j] <- Zi11_T[, wd1_list[[j]], drop=FALSE] %*% wd1_vec[[j]]
    vsum02[, j] <- Zi12_T[, wd1_list[[j]], drop=FALSE] %*% wd1_vec[[j]]
  }
  vsum1 <- vsum01 + alpha02; vsum2 <- vsum02 + alpha02   #�f�B���N�����z�̃p�����[�^
  phi1 <- extraDistr::rdirichlet(k1, vsum1)   #�f�B���N�����z����Ǝ��g�s�b�N�̒P�ꕪ�z���T���v�����O
  phi2 <- extraDistr::rdirichlet(k2, vsum2)   #�f�B���N�����z���狤�ʃg�s�b�N�̒P�ꕪ�z���T���v�����O
  
  #����2�̒P�ꕪ�z���T���v�����O
  vsum01 <- matrix(0, nrow=k2, ncol=v2)
  vsum02 <- rep(0, v2)
  for(j in 1:v2){
    vsum01[, j] <- Zi21_T[, wd2_list[[j]], drop=FALSE] %*% wd2_vec[[j]]
    vsum02[j] <- (1-flag2[wd2_list[[j]], drop=FALSE]) %*% wd2_vec[[j]]
  }
  vsum1 <- vsum01 + alpha02; vsum2 <- vsum02 + alpha02   #�f�B���N�����z�̃p�����[�^
  omega1 <- extraDistr::rdirichlet(k2, vsum1)   #�f�B���N�����z���狤�ʃg�s�b�N�̒P�ꕪ�z���T���v�����O
  omega2 <- as.numeric(extraDistr::rdirichlet(1, vsum2))   #�f�B���N�����z�����ʌ�̒P�ꕪ�z���T���v�����O
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    BETA1[mkeep, ] <- beta1
    BETA2[mkeep, ] <- beta2
    THETA1[, , mkeep] <- theta1
    THETA2[, , mkeep] <- theta2
    PHI1[, , mkeep] <- phi1
    PHI2[, , mkeep] <- phi2
    OMEGA1[, , mkeep] <- omega1
    OMEGA2[mkeep, ] <- omega2
  } 
  
  #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
  if(rp%%keep==0 & rp >= burnin){
    FLAG1 <- FLAG1 + flag1; FLAG2 <- FLAG2 + flag2
    SEG11 <- SEG11 + Zi11; SEG12 <- SEG12 + Zi12
    SEG21 <- SEG21 + Zi21
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    LL1 <- sum(log((1-flag_rate1)*topic_par11 + flag_rate1*topic_par12))
    LL2 <- sum(log(flag_rate2*topic_par21 + (1-flag_rate2)*topic_par22))
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(c(LL1+LL2, LL1, LL2, LLbest, LLbest1, LLbest2, LLst))
    print(round(c(mean(beta1), mean(beta2), mean(betat1), mean(betat2)), 3))
    print(round(cbind(phi1[, c(1:5, (v11+1):(v11+5))], phit1[, c(1:5, (v11+1):(v11+5))]), 3))
    print(round(rbind(omega2[c(1:10, (v21+1):(v21+10))], omegat2[c(1:10, (v21+1):(v21+10))]), 3))
  }
}

####�T���v�����O���ʂ̗v��Ɖ���####
burnin <- 2000/keep
RS <- R/keep

##�T���v�����O���ʂ̉���
#�g�s�b�N���z�̉���
matplot(t(THETA1[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA1[10, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA1[100, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA1[1000, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA2[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA2[10, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA2[100, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA2[1000, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

#�P�ꕪ�z�̉���
matplot(t(PHI1[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(PHI1[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(PHI2[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(PHI2[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(OMEGA1[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(OMEGA1[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(OMEGA2, type="l", xlab="�T���v�����O��", ylab="�p�����[�^")


##�T���v�����O���ʂ̎��㕪�z
#�g�s�b�N�����ߒ��̎��㕽��
round(cbind(colMeans(BETA1[burnin:RS, ]), betat1), 3)
round(cbind(colMeans(BETA2[burnin:RS, ]), betat2), 3)

#�g�s�b�N���z�̎��㕽��
round(cbind(apply(THETA1[, , burnin:RS], c(1, 2), mean), thetat1), 3)
round(cbind(apply(THETA2[, , burnin:RS], c(1, 2), mean), thetat2), 3)

#�P�ꕪ�z�̎��㕽��
round(cbind(t(apply(PHI1[, , burnin:RS], c(1, 2), mean)), t(phit1)), 3)
round(cbind(t(apply(PHI2[, , burnin:RS], c(1, 2), mean)), t(phit2)), 3)
round(cbind(t(apply(OMEGA1[, , burnin:RS], c(1, 2), mean)), t(omegat1)), 3)
round(cbind(colMeans(OMEGA2[burnin:RS, ]), omegat2), 3)


##���ݕϐ��̃T���v�����O���ʂ̎��㕪�z
seg11_rate <- SEG11 / rowSums(SEG11); seg12_rate <- SEG12 / rowSums(SEG12)
seg21_rate <- SEG21 / rowSums(SEG21)
seg11_rate[is.nan(seg11_rate)] <- 0; seg12_rate[is.nan(seg12_rate)] <- 0
seg21_rate[is.nan(seg21_rate)] <- 0

#�g�s�b�N�������ʂ��r
round(cbind(rowSums(SEG11), seg11_rate, Z11), 3)
round(cbind(rowSums(SEG12), seg12_rate, Z12), 3)
round(cbind(rowSums(SEG21), seg21_rate, Z21), 3)