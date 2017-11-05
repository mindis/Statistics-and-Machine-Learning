#####�p�`���R�z�����f��#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
detach("package:gtools", unload=TRUE)
library(bayesm)
library(ExtDist)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(76432)

####�f�[�^�̔���####
#set.seed(423943)
#�����f�[�^�̐ݒ�
d <- 2500   #������
v <- 250   #��b��
w <- rpois(d, rgamma(d, 160, 1.0))   #1����������̒P�ꐔ
hist(w, col="grey")

#�g�s�b�N��ݒ�
k1 <- 3   #��ʃg�s�b�N��
k2 <- 11   #���ʃg�s�b�N����
k21 <- 5   #���ʃg�s�b�N��1
k22 <- 4   #���ʃg�s�b�N��2
k23 <- 4   #���ʃg�s�b�N��3
k2v <- c(k21, k22, k23)

#�l�X�g�\����ݒ�
nest <- rbind(c(rep(1, k21), rep(0, k2-k21)),
              c(rep(0, 4), rep(1, k22), rep(0, 3)),
              c(rep(0, k2-k23), rep(1, k23)))


#ID�̐ݒ�
word_id <- rep(1:d, w)

##�p�����[�^�̐ݒ�
alpha01 <- runif(k1, 0.4, 0.7)   #�����̏�ʃg�s�b�N�̃f�B���N�����O���z�̃p�����[�^
alpha021 <- rep(0.25, k21)   #�����̉��ʃg�s�b�N1�̃f�B���N�����O���z�̃p�����[�^
alpha022 <- rep(0.3, k22)   #�����̉��ʃg�s�b�N2�̃f�B���N�����O���z�̃p�����[�^
alpha023 <- rep(0.3, k23)   #�����̉��ʃg�s�b�N2�̃f�B���N�����O���z�̃p�����[�^
alpha1 <- rep(0.15, v)   #�P��̃f�B���N�����O���z�̃p�����[�^

#�f�B���N�������̔���
theta01 <- extraDistr::rdirichlet(d, alpha01)
theta021 <- extraDistr::rdirichlet(d, alpha021)
theta022 <- extraDistr::rdirichlet(d, alpha022)
theta023 <- extraDistr::rdirichlet(d, alpha023)
theta02 <- list(theta021, theta022, theta023)
phi0 <- extraDistr::rdirichlet(k2, alpha1)   #�P��̃g�s�b�N�����f�B���N�����z���甭��

##�������z����g�s�b�N����ђP��f�[�^�𔭐�
WX <- matrix(0, nrow=d, ncol=v)
Z1 <- list()
Z2 <- list()

#�������ƂɃg�s�b�N�ƒP��𒀎�����
for(i in 1:d){
  print(i)
  
  #�����̏�ʃg�s�b�N���z�𔭐�
  z1 <- t(rmultinom(w[i], 1, theta01[i, ]))   #�����̏�ʃg�s�b�N���z�𔭐�
  z1_vec <- as.numeric(z1 %*% 1:k1)
  
  #�����̉��ʃg�s�b�N���z�𔭐�
  z2 <- matrix(0, nrow=length(z1_vec), ncol=k2)
  for(j in 1:length(z1_vec)){
    z2[j, nest[z1_vec[j], ]==1] <- as.numeric(t(rmultinom(1, 1, theta02[[z1_vec[j]]][i, ])))
  }
  z2
  #�����̉��ʃg�s�b�N���z����P��𔭐�
  zn <- z2 %*% c(1:k2)   #0,1�𐔒l�ɒu��������
  zdn <- cbind(zn, z2)   #apply�֐��Ŏg����悤�ɍs��ɂ��Ă���
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi0[x[1], ])))   #�����̃g�s�b�N����P��𐶐�
  wdn <- colSums(wn)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  WX[i, ] <- wdn  
  
  #�����������g�s�b�N���i�[
  Z1[[i]] <- z1_vec
  Z2[[i]] <- zdn[, 1]
}

#�����������f�[�^�̊m�F
colSums(WX); round(colMeans(WX), 3)
table(unlist(Z1)); table(unlist(Z2))
storage.mode(WX) <- "integer"   #�f�[�^�s��𐮐��^�s��ɕύX


####�g�s�b�N���f���̂��߂̃f�[�^�Ɗ֐��̏���####
##���ꂼ��̕������̒P��̏o������ѕ⏕���̏o�����x�N�g���ɕ��ׂ�
##�f�[�^����pID���쐬
ID_list <- list()
wd_list <- list()

#���l���Ƃɋ��lID����ђP��ID���쐬
for(i in 1:nrow(WX)){
  print(i)
  
  #�P���ID�x�N�g�����쐬
  ID_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
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
for(i in 1:length(unique(ID_d))) {doc_list[[i]] <- which(ID_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- which(wd==i)}
gc(); gc()


####�l�X�g�����ʂ��邽�߂̏����l��ݒ肷��####
##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, k){
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

##EM�A���S���Y���̏����l��ݒ肷��
##�����l�������_���ɐݒ�
#phi�̏����l
freq_v <- matrix(colSums(WX), nrow=k2, ncol=v, byrow=T)   #�P��̏o����
rand_v <- matrix(trunc(rnorm(k2*v, 0, (colSums(WX)/2))), nrow=k2, ncol=v, byrow=T)   #�����_����
phi_r <- abs(freq_v + rand_v) / rowSums(abs(freq_v + rand_v))   #�g�s�b�N���Ƃ̏o�����������_���ɏ�����

#theta�̏����l
theta_r <- rdirichlet(d, runif(k2, 0.2, 4))   #�f�B���N�����z���珉���l��ݒ�


###�p�����[�^�̍X�V
##���S���̌v�Z
bfr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, k=k2)
Br <- bfr$Br   #���S��
r <- bfr$r   #������

##theta�̍X�V
tsum <- (data.frame(id=ID_d, Br=Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_all(funs(sum)))[, 2:(k2+1)]
theta_r <- tsum / matrix(w, nrow=d, ncol=k2, byrow=T)   #�p�����[�^���v�Z

##phi�̍X�V
vf <- (data.frame(id=wd, Br=Br) %>%
         dplyr::group_by(id) %>%
         dplyr::summarize_all(funs(sum)))[, 2:(k2+1)]
phi_r <- t(vf) / matrix(colSums(vf), nrow=k2, ncol=v)

#�ΐ��ޓx�̌v�Z
(LLS <- sum(log(rowSums(bfr$Bur))))


####EM�A���S���Y���Ńp�����[�^���X�V####
#�X�V�X�e�[�^�X
iter <- 1
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 0.5
LLo <- LLS   #�ΐ��ޓx�̏����l
LLw <- LLS

###�p�����[�^�̍X�V
##���S���̌v�Z
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  bfr <- burden_fr(theta=theta_r, phi=phi_r, wd=wd, k=k2)
  Br <- bfr$Br   #���S��
  r <- bfr$r   #������
  
  ##theta�̍X�V
  tsum <- (data.frame(id=ID_d, Br=Br) %>%
             dplyr::group_by(id) %>%
             dplyr::summarize_all(funs(sum)))[, 2:(k2+1)]
  theta_r <- tsum / matrix(w, nrow=d, ncol=k2, byrow=T)   #�p�����[�^���v�Z
  
  ##phi�̍X�V
  vf <- (data.frame(id=wd, Br=Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_all(funs(sum)))[, 2:(k2+1)]
  phi_r <- t(vf) / matrix(colSums(vf), nrow=k2, ncol=v)
  
  #�ΐ��ޓx�̌v�Z
  LLS <- sum(log(rowSums(bfr$Bur)))
  
  iter <- iter+1
  dl <- LLS-LLo
  LLo <- LLS
  LLw <- c(LLw, LLo)
  print(LLo)
}


####�}���R�t�A�������e�J�����@�Ńp�`���R�z�����f���𐄒�####
##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta1, theta2, phi, wd, w, k1, k2, nest){
  Bur <- matrix(0, nrow=length(wd), ncol=k2)   #���S�W���̊i�[�p�z��
  for(j in 1:k1){
    for(k in 1:k2){
      if(nest[j, k]==0) next
      theta_nest <- theta2[[j]]
      Bur[, k] <- Bur[, k] + rep(theta1[, j], w) * rep(theta_nest[, sum(nest[j, 1:k]==1)], w) * phi[k, wd]
    }
  }
  Br <- Bur / rowSums(Bur)   #���S���̌v�Z
  r <- colSums(Br) / sum(Br)   #�������̌v�Z
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}


##�A���S���Y���̐ݒ�
R <- 10000   #�T���v�����O��
keep <- 2   #2���1��̊����ŃT���v�����O���ʂ��i�[
iter <- 0

##�l�X�g�̐ݒ�
sort_vec <- rep(0, k2)
for(i in 1:k2){
  er_sq <- rep(0, k2)
  for(j in 1:k2){
    er_sq[j] <- sum((phi_r[j, ] - phi0[i, ])^2)
  }
  er_sq
  sort_vec[i] <- which.min(er_sq)
}


##���O���z�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z
alpha01 <- rep(1.0, k1)
alpha02 <- rep(1.0, k2)
alpha01m <- matrix(alpha01, nrow=d, ncol=k1, byrow=T)
alpha02m <- matrix(alpha02, nrow=d, ncol=k2, byrow=T)
beta0 <- rep(0.5, v)
beta0m <- matrix(beta0, nrow=v, ncol=k2)


##�p�����[�^�̏����l
#��ʃg�s�b�N����щ��ʃg�s�b�N�̏����p�����[�^��ݒ�
Zi0 <- Br[, sort_vec]
theta_ho <- matrix(0, nrow=length(wd), ncol=k1)
wsum1 <- matrix(0, nrow=d, ncol=k1)
theta2 <- list()
theta1 <- matrix(1/k1, nrow=d, ncol=k1)   #�����̏�ʃg�s�b�N�̃p�����[�^�̏����l

for(j in 1:k1){
  theta2[[j]] <- (theta_r[, sort_vec])[, nest[j, ]==1]
  theta_ho[, j] <- rowSums((matrix(theta1[, j], nrow=d, ncol=sum(nest[j, ])) * theta2[[j]])[ID_d, ] * Zi0[, nest[j, ]==1])
}

theta_rate <- theta_ho / rowSums(theta_ho)   #�����m��

#�������z�����ʃg�s�b�N���T���v�����O
Zi1 <- rmnom(nrow(theta_rate), 1, theta_rate)

#��ʃg�s�b�N���z�̃p�����[�^���X�V
for(i in 1:d){wsum1[i, ] <- colSums(Zi1[doc_list[[i]], ])}
theta1 <- extraDistr::rdirichlet(d, wsum1 + alpha01m)

#�P��g�s�b�N�̃p�����[�^�̏����l
phi <- phi_r[sort_vec, ]   


##�p�����[�^�̊i�[�p�z��
THETA1 <- array(0, dim=c(d, k1, R/keep))
THETA21 <- array(0, dim=c(d, k21, R/keep))
THETA22 <- array(0, dim=c(d, k22, R/keep))
THETA23 <- array(0, dim=c(d, k23, R/keep))
PHI <- array(0, dim=c(k2, v, R/keep))
Z0_SEG <- matrix(0, nrow=length(wd), ncol=k2)
Z1_SEG <- matrix(0, nrow=length(wd), ncol=k1)
storage.mode(Z0_SEG) <- "integer"
storage.mode(Z1_SEG) <- "integer"
gc(); gc()


##MCMC����p�z��
wsum0 <- matrix(0, nrow=d, ncol=k2)
wsum1 <- matrix(0, nrow=d, ncol=k1)
wsum2 <- list()
for(j in 1:k1){
  wsum2[[j]] <- matrix(0, nrow=d, ncol=sum(nest[j, ]))
}
theta_ho <- matrix(0, nrow=length(ID_d), ncol=k1)
vf0 <- matrix(0, nrow=v, ncol=k2)


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�P�ꂲ�ƂɃg�s�b�N���T���v�����O
  #�P�ꂲ�ƂɃg�s�b�N�̏o�������v�Z
  fr <- burden_fr(theta1, theta2, phi, wd, w, k1, k2, nest)
  word_rate <- fr$Br
  
  #�������z����P��g�s�b�N���T���v�����O
  Zi0 <- rmnom(nrow(word_rate), 1, word_rate)
  
  #�������ƂɃg�s�b�N�o���m�����v�Z
  for(i in 1:d){wsum0[i, ] <- colSums(Zi0[doc_list[[i]], ])}
  
  
  ##��ʃg�s�b�N�p�����[�^���T���v�����O
  #���ʃg�s�b�N�̏o���������Ӊ��������S��
  for(j in 1:k1){
    theta_ho[, j] <- rowSums((matrix(theta1[, j], nrow=d, ncol=sum(nest[j, ])) * theta2[[j]])[ID_d, ] * Zi0[, nest[j, ]==1])
  }
  theta_rate <- theta_ho / rowSums(theta_ho)
  
  #�������z�����ʃg�s�b�N���T���v�����O
  Zi1 <- rmnom(nrow(theta_rate), 1, theta_rate)
  
  #��ʃg�s�b�N���z�̃p�����[�^���X�V
  for(i in 1:d){wsum1[i, ] <- colSums(Zi1[doc_list[[i]], ])}
  theta1 <- extraDistr::rdirichlet(d, wsum1 + alpha01m)
  
  ##���ʃg�s�b�N�p�����[�^���T���v�����O
  for(j in 1:k1){
    for(i in 1:d){
      #�f�B�N�������z�̃p�����[�^
      wsum2[[j]][i, ] <- colSums(Zi0[doc_list[[i]], nest[j, ]==1] * 
                                   matrix(Zi1[doc_list[[i]], j], nrow=length(doc_list[[i]]), ncol=sum(nest[j, ]))) + 1
    }
    
    #���ʃg�s�b�N���z�̃p�����[�^���X�V
    theta2[[j]] <- extraDistr::rdirichlet(d, wsum2[[j]])
  }
  
  ##�g�s�b�N���ƂɒP�ꕪ�z�̃p�����[�^���T���v�����O
  for(i in 1:v){vf0[i, ] <- colSums(Zi0[word_list[[i]], ])}
  vf <- vf0 + beta0m
  phi <- extraDistr::rdirichlet(k2, t(vf))
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA1[, , mkeep] <- theta1
    THETA21[, , mkeep] <- theta2[[1]]
    THETA22[, , mkeep] <- theta2[[2]]
    THETA23[, , mkeep] <- theta2[[3]]
    PHI[, , mkeep] <- phi
    
    #�g�s�b�N�����̓T���v�����O���Ԃ̔����𒴂�����i�[����
    if(rp >= R/2){
      Z0_SEG <- Z0_SEG + Zi0
      Z1_SEG <- Z1_SEG + Zi1
    }
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(round(cbind(theta1[1:12, ], theta01[1:12, ], wsum2[[1]][1:12, ], theta2[[1]][1:12, ], theta02[[1]][1:12, ]), 3))
    print(round(cbind(phi[, 1:10], phi0[, 1:10]), 3))
  }
}


####�T���v�����O���ʂ̉����Ɨv��####
burnin <- 2000   #�o�[���C������

##�T���v�����O���ʂ̉���
#�����̏�ʃg�s�b�N���z�̃T���v�����O����
matplot(t(THETA1[1, , ]), type="l", ylab="�p�����[�^", main="����1�̏�ʃg�s�b�N���z�̃T���v�����O����")
matplot(t(THETA1[2, , ]), type="l", ylab="�p�����[�^", main="����2�̏�ʃg�s�b�N���z�̃T���v�����O����")
matplot(t(THETA1[3, , ]), type="l", ylab="�p�����[�^", main="����3�̏�ʃg�s�b�N���z�̃T���v�����O����")
matplot(t(THETA1[4, , ]), type="l", ylab="�p�����[�^", main="����4�̏�ʃg�s�b�N���z�̃T���v�����O����")

#�����̉��ʃg�s�b�N���z�̃T���v�����O����
matplot(t(THETA21[1, , ]), type="l", ylab="�p�����[�^", main="����1�̉��ʃg�s�b�N���z�̃T���v�����O����")
matplot(t(THETA21[2, , ]), type="l", ylab="�p�����[�^", main="����2�̉��ʃg�s�b�N���z�̃T���v�����O����")
matplot(t(THETA21[3, , ]), type="l", ylab="�p�����[�^", main="����3�̉��ʃg�s�b�N���z�̃T���v�����O����")
matplot(t(THETA21[4, , ]), type="l", ylab="�p�����[�^", main="����4�̉��ʃg�s�b�N���z�̃T���v�����O����")

#�P��̏o���m���̃T���v�����O����
matplot(t(PHI[1, 1:10, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N1�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[2, 11:20, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N2�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[3, 21:30, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N3�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[4, 31:40, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N4�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[5, 41:50, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N5�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[6, 51:60, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N6�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[7, 61:70, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N7�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[8, 71:80, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N8�̒P��̏o�����̃T���v�����O����")


##�T���v�����O���ʂ̗v�񐄒��
#��ʃg�s�b�N���z�̎��㐄���
topic_mu1 <- apply(THETA1[, , burnin:(R/keep)], c(1, 2), mean)   #�g�s�b�N���z�̎��㕽��
round(cbind(topic_mu1, theta01), 3)
round(topic_sd1 <- apply(THETA1[, , burnin:(R/keep)], c(1, 2), sd), 3)   #�g�s�b�N���z�̎���W���΍�

#���ʃg�s�b�N���z�̎��㐄���
topic_mu21 <- apply(THETA21[, , burnin:(R/keep)], c(1, 2), mean)   #�g�s�b�N���z�̎��㕽��
topic_mu22 <- apply(THETA22[, , burnin:(R/keep)], c(1, 2), mean)   #�g�s�b�N���z�̎��㕽��
topic_mu23 <- apply(THETA23[, , burnin:(R/keep)], c(1, 2), mean)   #�g�s�b�N���z�̎��㕽��
round(cbind(topic_mu21, theta021), 3)
round(cbind(topic_mu22, theta022), 3)
round(cbind(topic_mu23, theta023), 3)
round(topic_sd21 <- apply(THETA21[, , burnin:(R/keep)], c(1, 2), sd), 3)   #�g�s�b�N���z�̎���W���΍�
round(topic_sd22 <- apply(THETA22[, , burnin:(R/keep)], c(1, 2), sd), 3)   #�g�s�b�N���z�̎���W���΍�
round(topic_sd23 <- apply(THETA23[, , burnin:(R/keep)], c(1, 2), sd), 3)   #�g�s�b�N���z�̎���W���΍�

#�P��o���m���̎��㐄���
word_mu <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #�P��̏o�����̎��㕽��
round(rbind(word_mu, phi0)[, 1:50], 3)



  