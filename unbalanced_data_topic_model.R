#####�s�ނ荇���f�[�^�̊m���I���݈Ӗ����(�g�s�b�N���f��)#####
library(MASS)
library(lda)
library(RMeCab)
detach("package:bayesm", unload=TRUE)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####�f�[�^�̔���####
#set.seed(423943)
#�f�[�^�̐ݒ�
k <- 8   #�g�s�b�N��
d <- 5000   #������
v <- 300   #��b��
w <- rpois(d, 250)   #1����������̒P�ꐔ

#�p�����[�^�̐ݒ�
alpha0 <- round(runif(k, 0.1, 1.25), 3)   #�����̃f�B���N�����O���z�̃p�����[�^
alpha1 <- rep(0.25, v)   #�P��̃f�B���N�����O���z�̃p�����[�^

#�f�B���N�������̔���
theta <- rdirichlet(d, alpha0)   #�����̃g�s�b�N���z���f�B���N���������甭��
phi <- rdirichlet(k, alpha1)   #�P��̃g�s�b�N���z���f�B���N���������甭��

#�������z�̗�������f�[�^�𔭐�
WX <- matrix(0, nrow=d, ncol=v)
Z <- list()
for(i in 1:d){
  z <- t(rmultinom(w[i], 1, theta[i, ]))   #�����̃g�s�b�N���z�𔭐�
  zn <- z %*% c(1:k)   #0,1�𐔒l�ɒu��������
  zdn <- cbind(zn, z)   #apply�֐��Ŏg����悤�ɍs��ɂ��Ă���
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi[x[1], ])))   #�����̃g�s�b�N����P��𐶐�
  wdn <- colSums(wn)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  WX[i, ] <- wdn  
  Z[[i]] <- zdn[, 1]
  print(i)
}

####EM�A���S���Y���Ńg�s�b�N���f���𐄒�####
####�g�s�b�N���f���̂��߂̃f�[�^�Ɗ֐��̏���####
##���ꂼ��̕������̒P��̏o�����x�N�g���ɕ��ׂ�
##�f�[�^����pID���쐬
ID_list <- list()
wd_list <- list()

#���l���Ƃɋ��lID����ђP��ID���쐬
for(i in 1:nrow(WX)){
  print(i)
  ID_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * c(1:v) 
  num2 <- subset(num1, num1 > 0)
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
for(i in 1:length(unique(ID_d))) {doc_list[[i]] <- subset(1:length(ID_d), ID_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- subset(1:length(wd), wd==i)}
gc(); gc()


##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #���S�W���̊i�[�p
  for(kk in 1:k){
    #���S�W�����v�Z
    Bi <- rep(theta[, kk], w) * phi[kk, c(wd)]   #�ޓx
    Bur[, kk] <- Bi   
  }
  Br <- Bur / rowSums(Bur)   #���S���̌v�Z
  r <- colSums(Br) / sum(Br)   #�������̌v�Z
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}


####EM�A���S���Y���̏����l��ݒ肷��####
##�����l�������_���ɐݒ�
#phi�̏����l
freq_v <- matrix(colSums(WX), nrow=k, ncol=v, byrow=T)   #�P��̏o����
rand_v <- matrix(trunc(rnorm(k*v, 0, (colSums(WX)/2))), nrow=k, ncol=v, byrow=T)   #�����_����
phi_r <- abs(freq_v + rand_v) / rowSums(abs(freq_v + rand_v))   #�g�s�b�N���Ƃ̏o�����������_���ɏ�����

#theta�̏����l
theta_r <- rdirichlet(d, runif(k, 0.2, 4))   #�f�B���N�����z���珉���l��ݒ�

###�p�����[�^�̍X�V
##���S���̌v�Z
bfr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, w=w, k=k)
Br <- bfr$Br   #���S��
r <- bfr$r   #������

##theta�̍X�V
tsum <- (data.frame(id=ID_d, Br=Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
theta_r <- tsum / matrix(w, nrow=d, ncol=k)   #�p�����[�^���v�Z

##phi�̍X�V
vf <- (data.frame(id=wd, Br=Br) %>%
         dplyr::group_by(id) %>%
         dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
phi_r <- t(vf) / matrix(colSums(vf), nrow=k, ncol=v)

#�ΐ��ޓx�̌v�Z
(LLS <- sum(log(rowSums(bfr$Bur))))


####EM�A���S���Y���Ńp�����[�^���X�V####
#�X�V�X�e�[�^�X
iter <- 1
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 0.1
LLo <- LLS   #�ΐ��ޓx�̏����l
LLw <- LLS

out <- cbind(ID_d, wd, round(bfr$Br, 3))

###�p�����[�^�̍X�V
##���S���̌v�Z
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  bfr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, w=w, k=k)
  Br <- bfr$Br   #���S��
  r <- bfr$r   #������
  
  ##theta�̍X�V
  tsum <- (data.frame(id=ID_d, Br=Br) %>%
             dplyr::group_by(id) %>%
             dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
  theta_r <- tsum / matrix(w, nrow=d, ncol=k)   #�p�����[�^���v�Z

  ##phi�̍X�V
  vf <- (data.frame(id=wd, Br=Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
  phi_r <- t(vf) / matrix(colSums(vf), nrow=k, ncol=v)
  
  #�ΐ��ޓx�̌v�Z
  LLS <- sum(log(rowSums(bfr$Bur)))
  
  iter <- iter+1
  dl <- LLS-LLo
  LLo <- LLS
  LLw <- c(LLw, LLo)
  print(LLo)
}

####���茋�ʂƓ��v��####
plot(1:length(LLw), LLw, type="l", xlab="iter", ylab="LL", main="�ΐ��ޓx�̕ω�", lwd=2)

(PHI <- data.frame(round(t(phi_r), 3), t=round(t(phi), 3)))   #phi�̐^�̒l�Ɛ��茋�ʂ̔�r
(THETA <- data.frame(w, round(theta_r, 3), t=round(theta, 3)))   #theta�̐^�̒l�Ɛ��茋�ʂ̔�r
r   #�������̐��茋��

round(colSums(THETA[, 1:k]) / sum(THETA[, 1:k]), 3)   #���肳�ꂽ�������̊e�g�s�b�N�̔䗦
round(colSums(THETA[, (k+1):(2*k)]) / sum(THETA[, (k+1):(2*k)]), 3)   #�^�̕������̊e�g�s�b�N�̔䗦

#AIC��BIC
tp <- dim(theta_r)[1]*dim(theta_r)[2]
pp <- dim(phi_r)[1]*dim(phi_r)[2]

(AIC <- -2*LLo + 2*(tp+pp)) 
(BIC <- -2*LLo + log(nrow(WX))*(tp+pp))

##���ʂ��O���t��
#theta�̃v���b�g(50�Ԗڂ̕����܂�)
barplot(theta_r[1:100, 1], ylim=c(0, 1), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(theta_r[, 1]), col=10, lty=5)
barplot(theta_r[1:100, 2], ylim=c(0, 1), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(theta_r[, 2]), col=10, lty=5)
barplot(theta_r[1:100, 3], ylim=c(0, 1), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(theta_r[, 2]), col=10, lty=5)
barplot(theta_r[1:100, 4], ylim=c(0, 1), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(theta_r[, 2]), col=10, lty=5)
barplot(theta_r[1:100, 5], ylim=c(0, 1), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(theta_r[, 2]), col=10, lty=5)

#phi�̃v���b�g(50�Ԗڂ̒P��܂�)
barplot(phi_r[1, 1:50], ylim=c(0, 0.05), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(phi_r[1, ]), col=10, lty=5)
barplot(phi_r[2, 1:50], ylim=c(0, 0.05), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(phi_r[2, ]), col=10, lty=5)
barplot(phi_r[3, 1:50], ylim=c(0, 0.05), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(phi_r[3, ]), col=10, lty=5)
barplot(phi_r[4, 1:50], ylim=c(0, 0.05), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(phi_r[4, ]), col=10, lty=5)
barplot(phi_r[5, 1:50], ylim=c(0, 0.05), col=c(1:ncol(phi_r)), density=50)
abline(h=mean(phi_r[5, ]), col=10, lty=5)