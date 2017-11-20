#####�����m���I���݈Ӗ����#####
library(MASS)
library(lda)
library(RMeCab)
detach("package:bayesm", unload=TRUE)
library(extraDistr)
library(matrixStats)
library(monomvn)
library(lars)
library(glmnet)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(58079)

####�f�[�^�̔���####
#set.seed(423943)
#�����f�[�^�̐ݒ�
k <- 15   #�g�s�b�N��
d <- 8000   #������
v <- 450   #��b��
w0 <- rpois(d, rgamma(d, 12.5, 0.175))   #1����������̒P�ꐔ
w <- ifelse(w0 < 20, 20, w0)
a <- 150   #�⏕�ϐ���
x0 <- rpois(d, rgamma(d, 12.5, 1.25))   #1�⏕����������̒P�ꐔ
x <- ifelse(x0 < 2, 2, x0)
select <- 8   #�����ϐ���


#�p�����[�^�̐ݒ�
alpha0 <- rep(0.2, k)   #�����̃f�B���N�����O���z�̃p�����[�^
alpha1 <- rep(0.15, v)   #�P��̃f�B���N�����O���z�̃p�����[�^
alpha2 <- rep(0.1, a)   #�⏕�f�[�^�̃f�B�N�������O���z�̃p�����[�^

#�f�B���N�������̔���
theta <- extraDistr::rdirichlet(d, alpha0)   #�����̃g�s�b�N���z���f�B���N���������甭��
phi <- extraDistr::rdirichlet(k, alpha1)   #�P��̃g�s�b�N���z���f�B���N���������甭��
omega <- extraDistr::rdirichlet(k, alpha2)   #�⏕�f�[�^�̃g�s�b�N���z���f�B�N�����������甭��


#�������z�̗�������f�[�^�𔭐�
WX <- matrix(0, nrow=d, ncol=v)
AX <- matrix(0, nrow=d, ncol=a)
Z1 <- list()
Z2 <- list()

for(i in 1:d){
  print(i)
  
  #�����̃g�s�b�N���z�𔭐�
  z1 <- t(rmultinom(w[i], 1, theta[i, ]))   #�����̃g�s�b�N���z�𔭐�
  
  #�����̃g�s�b�N���z����P��𔭐�
  zn <- z1 %*% c(1:k)   #0,1�𐔒l�ɒu��������
  zdn <- cbind(zn, z1)   #apply�֐��Ŏg����悤�ɍs��ɂ��Ă���
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi[x[1], ])))   #�����̃g�s�b�N����P��𐶐�
  wdn <- colSums(wn)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  WX[i, ] <- wdn  
  
  #�����̃g�s�b�N���z����⏕�ϐ��𔭐�
  z2 <- t(rmultinom(x[i], 1, theta[i, ]))
  zx <- z2 %*% 1:k
  zax <- cbind(zx, z2)
  an <- t(apply(zax, 1, function(x) rmultinom(1, 1, omega[x[1], ])))
  adn <- colSums(an)
  AX[i, ] <- adn
  
  #�����g�s�b�N����ѕ⏕���g�s�b�N���i�[
  Z1[[i]] <- z1
  Z2[[i]] <- z2
}

####�����ϐ��̔���####
#�����ϐ��̊i�[�p�z��
y <- matrix(0, nrow=d, ncol=select)
Pr <- matrix(0, nrow=d, ncol=select)
Pr0 <- matrix(0, nrow=d, ncol=select)

##�Ó��ȉ����ϐ�����������܂Ŕ���������
for(j in 1:5000){
  ##�p�����[�^�̐ݒ�
  #�g�s�b�N���f���̃p�����[�^
  sparse1 <- matrix(rbinom((select-1)*k, 1, 0.4), nrow=k, ncol=select-1)   #�p�����[�^�̃X�p�[�X�s��
  b00 <- runif(select-1, -0.5, 0.5)
  b01 <- (matrix(runif((select-1)*k, -4.0, 4.0), nrow=k, ncol=select-1)) * sparse1
  b02 <- (b01 + mvrnorm(k, rep(0, select-1), diag(0.2, select-1))) * sparse1
  b0 <- rbind(b00, b01, b02)
  rownames(b0) <- NULL
  
  #�P��̕ϗʌ��ʂ̃p�����[�^
  sparse2 <- matrix(rbinom((select-1)*v, 1, 0.3), nrow=v, ncol=select-1)   #�p�����[�^�̃X�p�[�X�s��
  cov0 <- diag(runif(select-1, 0.025, 0.25))
  a0 <- mvrnorm(v, rep(0, select-1), cov0) * sparse2

  ##�������ƂɊm���Ɖ����ϐ��𔭐�
  for(i in 1:d){
    logit <- c(c(1, log(colSums(Z1[[i]])+1), log(colSums(Z2[[i]])+1)) %*% b0, 0)
    Pr[i, ] <- exp(logit) / sum(exp(logit))
    y[i, ] <- t(rmultinom(1, 1, Pr[i, ]))
  }
  
  t1 <- sum(apply(Pr, 1, which.max)==y %*% 1:select)/d
  print(round(c(t1, min(colSums(y))), 3))
  
  if(t1 > 0.85 &  min(colSums(y)) > 250) break
}


####EM�A���S���Y���Ńg�s�b�N���f���𐄒�####
####�g�s�b�N���f���̂��߂̃f�[�^�Ɗ֐��̏���####
##���ꂼ��̕������̒P��̏o������ѕ⏕���̏o�����x�N�g���ɕ��ׂ�
##�f�[�^����pID���쐬
ID1_list <- list()
wd_list <- list()
ID2_list <- list()
ad_list <- list()

#���l���Ƃɋ��lID����ђP��ID���쐬
for(i in 1:nrow(WX)){
  print(i)
  
  #�P���ID�x�N�g�����쐬
  ID1_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
  num2 <- subset(num1, num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
  
  #�⏕����ID�x�N�g�����쐬
  ID2_list[[i]] <- rep(i, x[i])
  num1 <- (AX[i, ] > 0) * (1:a)
  num2 <- subset(num1, num1 > 0)
  A1 <- AX[i, (AX[i, ] > 0)]
  number <- rep(num2, A1)
  ad_list[[i]] <- number
}

#���X�g���x�N�g���ɕϊ�
ID1_d <- unlist(ID1_list)
ID2_d <- unlist(ID2_list)
wd <- unlist(wd_list)
ad <- unlist(ad_list)

##�C���f�b�N�X���쐬
doc1_list <- list()
word_list <- list()
doc2_list <- list()
aux_list <- list()
for(i in 1:length(unique(ID1_d))) {doc1_list[[i]] <- subset(1:length(ID1_d), ID1_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- subset(1:length(wd), wd==i)}
for(i in 1:length(unique(ID2_d))) {doc2_list[[i]] <- subset(1:length(ID2_d), ID2_d==i)}
for(i in 1:length(unique(ad))) {aux_list[[i]] <- subset(1:length(ad), ad==i)}
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

#omega�̏����l
freq_v <- matrix(colSums(AX)/sum(AX), nrow=k, ncol=a, byrow=T)
rand_v <- matrix(trunc(rnorm(k*a, 0, (colSums(AX)/2))), nrow=k, ncol=a, byrow=T)   #�����_����
omega_r <- abs(freq_v + rand_v) / rowSums(abs(freq_v + rand_v))   #�g�s�b�N���Ƃ̏o�����������_���ɏ�����

###�p�����[�^�̍X�V
##�P�ꃌ�x���̕��S���ƃp�����[�^��������
#�P�ꃌ�x���̕��S���̍X�V
word_fr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, w=w, k=k)
Bw <- word_fr$Br   #���S��
r1 <- word_fr$r   #������

#theta�̍X�V
wsum <- (data.frame(id=ID1_d, Br=Bw) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_all(funs(sum)))[, 2:(k+1)]
theta_r <- wsum / matrix(w, nrow=d, ncol=k)   #�p�����[�^���v�Z

##phi�̍X�V
vf <- (data.frame(id=wd, Br=Bw) %>%
         dplyr::group_by(id) %>%
         dplyr::summarize_all(funs(sum)))[, 2:(k+1)]
phi_r <- t(vf) / matrix(colSums(vf), nrow=k, ncol=v)


##�⏕��񃌃x���̕��S���ƃp�����[�^��������
#�⏕��񃌃x���̕��S���̍X�V
aux_fr <- burden_fr(theta=theta_r, phi=omega_r, wd=ad, w=x, k=k)
Ba <- aux_fr$Br   #���S��
r2 <- aux_fr$r   #������

##omega�̍X�V
af <- (data.frame(id=ad, Br=Ba) %>%
         dplyr::group_by(id) %>%
         dplyr::summarize_all(funs(sum)))[, 2:(k+1)]
omega_r <- t(af) / matrix(colSums(af), nrow=k, ncol=a)


#�ΐ��ޓx�̌v�Z
LLw <- sum(log(rowSums(word_fr$Bur)))   #�P�ꃌ�x���̑ΐ��ޓx
LLa <- sum(log(rowSums(aux_fr$Bur)))   #�⏕��񃌃x���̑ΐ��ޓx
LL <- LLw + LLa


####EM�A���S���Y���Ńp�����[�^���X�V####
#�X�V�X�e�[�^�X
iter <- 1
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 0.1
LL1 <- LL   #�ΐ��ޓx�̏����l
LLs <- c()

##EM�A���S���Y������������܂Ŕ���������
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  
  ##�P�ꃌ�x���̃p�����[�^���Ŗސ���
  #�P�ꃌ�x���̕��S���̍X�V
  word_fr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, w=w, k=k)
  Bw <- word_fr$Br   #���S��
  r1 <- word_fr$r   #������
  
  #theta�̍X�V
  wsum <- (data.frame(id=ID1_d, Br=Bw) %>%
             dplyr::group_by(id) %>%
             dplyr::summarize_all(funs(sum)))[, 2:(k+1)]
  theta_r <- wsum / matrix(w, nrow=d, ncol=k)   #�p�����[�^���v�Z
  
  ##phi�̍X�V
  vf <- (data.frame(id=wd, Br=Bw) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_all(funs(sum)))[, 2:(k+1)]
  phi_r <- t(vf) / matrix(colSums(vf), nrow=k, ncol=v)
  
  
  ##�⏕��񃌃x���̃p�����[�^���Ŗސ���
  #�⏕��񃌃x���̕��S���̍X�V
  aux_fr <- burden_fr(theta=theta_r, phi=omega_r, wd=ad, w=x, k=k)
  Ba <- aux_fr$Br   #���S��
  r2 <- aux_fr$r   #������
  
  ##omega�̍X�V
  af <- (data.frame(id=ad, Br=Ba) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_all(funs(sum)))[, 2:(k+1)]
  omega_r <- t(af) / matrix(colSums(af), nrow=k, ncol=a)
  
  
  ##�ϑ��f�[�^�̑ΐ��ޓx�̌v�Z
  LLw <- sum(log(rowSums(word_fr$Bur)))   #�P�ꃌ�x���̑ΐ��ޓx
  LLa <- sum(log(rowSums(aux_fr$Bur)))   #�⏕��񃌃x���̑ΐ��ޓx
  LL <- LLw + LLa
  
  ##�A���S���Y���̍X�V
  iter <- iter+1
  dl <- LL1 - LL
  LL1 <- LL
  LLs <- c(LLs, LL)
  print(LL)
}


####���茋�ʂƓ��v��####
plot(1:length(LLs), LLs, type="l", xlab="iter", ylab="LL", main="�ΐ��ޓx�̕ω�", lwd=2)

(PHI <- data.frame(round(t(phi_r), 3), t=round(t(phi), 3)))   #phi�̐^�̒l�Ɛ��茋�ʂ̔�r
(OMEGA <- data.frame(round(t(omega_r), 3), t=round(t(omega), 3)))   #omega�̐^�̒l�Ɛ��茋�ʂ̔�r
(THETA <- data.frame(w, round(theta_r, 3), t=round(theta, 3)))   #theta�̐^�̒l�Ɛ��茋�ʂ̔�r
r   #�������̐��茋��

round(colSums(THETA[, 2:(k+1)]) / sum(THETA[, 2:(k+1)]), 3)   #���肳�ꂽ�������̊e�g�s�b�N�̔䗦
round(colSums(THETA[, (k+1):(2*k)]) / sum(THETA[, (k+1):(2*k)]), 3)   #�^�̕������̊e�g�s�b�N�̔䗦

#AIC��BIC
tp <- dim(theta_r)[1]*dim(theta_r)[2] 
pp <- dim(phi_r)[1]*dim(phi_r)[2] + dim(omega_r)[1]*dim(omega_r)[2]

(AIC <- -2*LL + 2*(tp+pp)) 
(BIC <- -2*LL + log(nrow(WX))*(tp+pp))

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

#omega�̃v���b�g
barplot(omega_r[1, ], col=c(1:ncol(omega_r)), density=50)
abline(h=mean(omega_r[1, ]), col=10, lty=5)
barplot(omega_r[2, ], col=c(1:ncol(omega_r)), density=50)
abline(h=mean(omega_r[2, ]), col=10, lty=5)
barplot(omega_r[3, ], col=c(1:ncol(omega_r)), density=50)
abline(h=mean(omega_r[3, ]), col=10, lty=5)
barplot(omega_r[4, ], col=c(1:ncol(omega_r)), density=50)
abline(h=mean(omega_r[4, ]), col=10, lty=5)
barplot(omega_r[5, ], col=c(1:ncol(omega_r)), density=50)
abline(h=mean(omega_r[5, ]), col=10, lty=5)


####�K�w�x�C�Y�������W�b�g���f���ŕ��ރ��f�����쐬####
####�f�[�^�̐ݒ�####
##���肳�ꂽ�g�s�b�N���牞���ϐ����쐬
#�����g�s�b�N�̏o�����W�v
W0 <- as.matrix((data.frame(id=ID1_d, Br=Bw) %>%
                   dplyr::group_by(id) %>%
                   dplyr::summarize_all(funs(sum))))[, 2:(k+1)]

#�^�C�g���g�s�b�N�̏o�����W�v
A0 <- as.matrix((data.frame(id=ID2_d, Ba=Ba) %>%
                   dplyr::group_by(id) %>%
                   dplyr::summarize_all(funs(sum))))[, 2:(k+1)]

#�f�[�^���������Đ����ϐ��Ƃ���
Data1 <- cbind(1, log(W0 + 1), log(A0 + 1))
Data <- matrix(0, nrow=d, ncol=3*k+1)
round(cbind(theta_r, theta), 3)

for(i in 1:d){
  Data[i, ] <- c(1, colSums(Z1[[i]]), log(colSums(Z1[[i]])+1), log(colSums(Z2[[i]])+1))
}
round(cbind(W0, Data[, 2:21]), 0)

####�}���R�t�A�������e�J�����@�̐ݒ�####
##�������W�b�g���f���̑ΐ��ޓx
fr <- function(beta, y, x, hh, select){
  
  #���W�b�g�Ɗm���̌v�Z
  logit <- t(x %*% t(beta))
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #�ΐ��ޓx��ݒ�
  LLi <- rowSums(y*log(Pr)) 
  return(LLi)
}

##�w�K�f�[�^�ƃe�X�g�f�[�^�ɕ���
index_test <- sample(1:nrow(Data), 1000)
n1 <- d-length(index_test)
n2 <- length(index_test)
Data_train <- Data[-index_test, ]
Data_test <- Data[index_test, ]
y_train <- y[-index_test, ]
y_test <- y[index_test, ]
y_vec <- y_test %*% 1:select

##�A���S���Y���̐ݒ�
R <- 20000
keep <- 4
sbeta <- 1.5
iter <- 0

##�����l�̐ݒ�
par <- 2*k+1
beta0 <- scale(colSums(y_train))
oldbeta <- mvrnorm(n1, beta0[-select], diag(0.2, select-1))
oldtheta <- matrix(0, nrow=par, ncol=select-1)

oldcov <- diag(0.1, select-1)
inv_cov <- solve(oldcov)
mu <- Data_train %*% oldtheta

##�T���v�����O���ʂ̕ۑ��p�z��
THETA <- array(0, dim=c(par, select-1, R/keep))
COV <- array(0, dim=c(select-1, select-1, R/keep))

##�A���S���Y������p�z��
lognew <- rep(0, n1)
logold <- rep(0, n1)
logpnew <- rep(0, n1)
logpold <- rep(0, n1)
lambda <- c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001)
x_diag <- diag(select)[, -select]
er_new <- matrix(0, nrow=n1, select-1)
er_old <- matrix(0, nrow=n1, select-1)


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�������W�b�g���f���̃p�����[�^���T���v�����O
  #�V�����p�����[�^���T���v�����O
  betad <- oldbeta
  betan <- betad + mvrnorm(n1, rep(0, select-1), diag(0.015, select-1))

  #�덷��ݒ�
  er_new <- betan - mu
  er_old <- betad - mu

  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew <- fr(betan, y_train, x_diag, n1, select)
  logold <- fr(betad, y_train, x_diag, n1, select)
  logpnew <- -0.5 * rowSums(er_new %*% inv_cov * er_new)
  logpold <- -0.5 * rowSums(er_old %*% inv_cov * er_old)  
  
  #���g���|���X�w�C�X�e�B���O�@�Ńp�����[�^�̍̑�������
  rand <- runif(n1)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew + logpnew - logold - logpold)   #�̑𗦂��v�Z
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=n1, ncol=select-1)
  oldbeta <- flag*betan + (1-flag)*betad   #alpha��rand�������Ă�����̑�

  ##lasso�ŊK�w���f���̉�A�p�����[�^���T���v�����O
  for(j in 1:(select-1)){
    if(lambda[j] > 0.5){
      lambda[j] <- 0.001
    }
    res <- blasso(X=Data_train[, -1], y=oldbeta[, j], beta=oldtheta[-1, j], lambda2=lambda[j], s2=diag(oldcov)[j], 
                  normalize=TRUE, T=2)
    oldtheta[, j] <- c(res$mu[2], res$beta[2, ])
    lambda[j] <- res$lambda2[2]
  }
 
  ##�K�w���f���̕��U�����U�s��𐄒�
  mu <- Data_train %*% oldtheta
  oldcov <- var(oldbeta - mu)
  inv_cov <- solve(oldcov)
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- oldtheta
    COV[, , mkeep] <- oldcov
    print(rp)
    print(round(lambda, 3))
    print(round(mean(alpha), 3))
    print(sum(lognew))
    print(round(cbind(oldtheta[1:15, ], b0[1:15, ]), 3))

    ##�\�����z�𐄒�
    logit <- Data_test %*% oldtheta
    Pr <- exp(logit) / rowSums(exp(logit))
    print(mean(y_vec==apply(Pr, 1, which.max)))
  }
}

colSums(y)
round(cbind(oldtheta, b0), 3)



      