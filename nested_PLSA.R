#####����q�^�g�s�b�N���f��#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
detach("package:bayesm", unload=TRUE)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####�f�[�^�̔���####
#set.seed(423943)
#�f�[�^�����ݒ�
hh0 <- 500   #���[�U�[��
item0 <- 150   #�A�C�e����

##ID�ƃ��r���[�����𔭐�
#ID�����ݒ�
u.id0 <- rep(1:hh0, rep(item0, hh0))
i.id0 <- rep(1:item0, hh0)

#���r���[�����𔭐�
hist <- rep(0, hh0*item0)
for(i in 1:item0){
  p <- runif(1, 0.25, 0.5)
  hist[i.id0==i] <- rbinom(hh0, 1, p)
}

#���r���[��������ID���Đݒ�
index <- subset(1:length(hist), hist==1)
u.id <- u.id0[index]
i.id <- i.id0[index]
ID <- data.frame(no=1:length(u.id), id=u.id, item=i.id)

#�f�[�^�̍Đݒ�
k <- 8   #�g�s�b�N��
hh <- length(unique(u.id))   #���[�U�[��
item <- length(unique(i.id))   #�A�C�e����
d <- length(u.id)   #������
v <- 250   #��b��

#1����������̒P�ꐔ
freq <- as.numeric(table(u.id))
w <- c()
for(i in 1:hh){
  par <- rgamma(1, 200, 1.10)
  w <- c(w, rpois(freq[i], par))   #�w����
}

####bag of word�`���̕����s��𔭐�####
#�p�����[�^�̐ݒ�
alpha0 <- round(runif(k, 0.15, 0.25), 3)   #���[�U�[�̃f�B���N�����O���z�̃p�����[�^
alpha1 <- round(runif(k, 0.15, 0.25), 3)   #�A�C�e���̃f�B���N�����O���z�̃p�����[�^
alpha2 <- rep(0.2, v)   #�P��̃f�B���N�����O���z�̃p�����[�^

#�f�B���N�������̔���
beta <- rdirichlet(hh, alpha0)   #���[�U�[�̃g�s�b�N���z���f�B���N���������甭��
gamma <- rdirichlet(item, alpha1)   #�A�C�e���̃g�s�b�N���z���f�B���N���������甭�� 
phi <- rdirichlet(k, alpha2)   #�P��̃g�s�b�N���z���f�B���N���������甭��
pi <- rbeta(sum(w), 0.15, 0.15)


#�X�C�b�`���O���x���I���̂��߂̃C���f�b�N�X���쐬
index_list <- list()
index_pi <- list()
for(i in 1:d){index_list[[i]] <- rep(i, w[i])}
index <- unlist(index_list)
for(i in 1:d){
  print(i)
  index_pi[[i]] <- which(index==i)
}


#�������z�̗�������f�[�^�𔭐�
WX <- matrix(0, nrow=d, ncol=v)
Pi <- matrix(0, nrow=d, ncol=v)
theta1 <- matrix(0, nrow=d, ncol=k)
theta2 <- matrix(0, nrow=d, ncol=k)
Z <- list()
y <- list()

for(i in 1:hh){
  
  r1 <- i.id[u.id==i]
  u1 <- u.id[u.id==i]
  index1 <- subset(1:length(u.id), u.id==i)
  
  for(j in 1:sum(u.id==i)){
    
    #�C���f�b�N�X��ݒ�
    r2 <- r1[j]   #�A�C�e���C���f�b�N�X
    u2 <- u1[j]   #���[�U�[�C���f�b�N�X
    index2 <- index1[j]   #�����C���f�b�N�X
    
    #�������z��蕶���g�s�b�N�𐶐�
    #�����̃g�s�b�N���z�𔭐�
    z1 <- t(rmultinom(w[index2], 1, beta[u2, ]))   #���[�U�[�̃g�s�b�N�𔭐�   
    z2 <- t(rmultinom(w[index2], 1, gamma[r2, ]))   #�A�C�e���̃g�s�b�N�̔���
    
    #�ǂ���̃g�s�b�N���甭��������������
    y[[index2]] <- rbinom(w[index2], 1, pi[index_pi[[index2]]])
    Y <- matrix(y[[index2]], nrow=w[index2], ncol=k)
    z <- Y*z1 + (1-Y)*z2
    
    #�g�s�b�N���x�N�g���`���ɕϊ�
    zn <- z %*% c(1:k)   #0,1�𐔒l�ɒu��������
    zdn <- cbind(zn, z)   #apply�֐��Ŏg����悤�ɍs��ɂ��Ă���
    
    #�����g�s�b�N���P��𐶐�
    wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi[x[1], ])))   #�����̃g�s�b�N����P��𐶐�
    wdn <- colSums(wn)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
    WX[index2, ] <- wdn  
    Z[[index2]] <- zdn[, 1]
    
    #�g�s�b�N�̔����m�����i�[���Ă���
    freq_w0 <- colSums(wn)
    freq_w <- ifelse(freq_w0==0, 1, freq_w0)
    Pi[index2, ] <- colSums(matrix(pi[index_pi[[index2]]], nrow=length(pi[index_pi[[index2]]]), ncol=v) * wn) / freq_w
  }
  print(i)
}

storage.mode(WX) <- "integer"
gc(); gc()

####EM�A���S���Y���Ńg�s�b�N���f���𐄒�####
####�g�s�b�N���f���̂��߂̃f�[�^�Ɗ֐��̏���####
##���ꂼ��̕������̒P��̏o�����x�N�g���ɕ��ׂ�
##�f�[�^����pID���쐬
ID_list <- list()
pi_list <- list()
y_list <- list()
user_list <- list()
item_list <- list()
wd_list <- list()

#�������ƂɃ��[�U�[ID����ђP��ID���쐬
for(i in 1:nrow(WX)){
  print(i)
  
  #ID���i�[
  ID_list[[i]] <- rep(i, w[i])
  user_list[[i]] <- rep(ID$id[i], w[i])
  item_list[[i]] <- rep(ID$item[i], w[i])
  
  #�o���P����x�N�g�������Ċi�[
  num11 <- (WX[i, ] > 0) * c(1:v) 
  num12 <- subset(num11, num11 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number1 <- rep(num12, W1)
  wd_list[[i]] <- number1
  
  #�g�s�b�N�����m�����x�N�g�������Ċi�[
  num21 <- Pi[i, ] * rep(1, v) 
  num22 <- subset(num21, num21 > 0)
  P1 <- WX[i, (WX[i, ] > 0)]
  number2 <- rep(num22, P1)
  pi_list[[i]] <- number2
}

#���X�g���x�N�g���ɕϊ�
ID_d <- unlist(ID_list)
ID_u <- unlist(user_list)
ID_i <- unlist(item_list)
wd <- unlist(wd_list)
pd <- unlist(pi_list)

##�C���f�b�N�X���쐬
doc_list <- list()
user_list <- list()
item_list <- list()
word_list <- list()

for(i in 1:length(unique(ID_d))) {doc_list[[i]] <- which(ID_d==i)}
for(i in 1:length(unique(ID_u))) {user_list[[i]] <- which(ID_u==i)}
for(i in 1:length(unique(ID_i))) {item_list[[i]] <- which(ID_i==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- which(wd==i)}
gc(); gc()

##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(beta, gamma, phi, r, ID_u, ID_i, wd, w, k){
  #�ޓx�̊i�[�p�z��
  Bur1 <- matrix(0, nrow=length(wd), ncol=k)
  Bur2 <- matrix(0, nrow=length(wd), ncol=k)
  
  #���[�U�[�A�A�C�e�����ƂɃg�s�b�N���z�̖ޓx���v�Z
  for(kk in 1:k){
    #���S�W�����v�Z
    phi_vec <- phi[kk, wd]
    Bi1 <- beta[ID_u, kk] * phi_vec   #���[�U�[�̖ޓx
    Bi2 <- gamma[ID_i, kk] * phi_vec   #�A�C�e���̖ޓx
    
    #���[�U�A�A�C�e�����Ƃɕ��S�W�����i�[
    Bur1[, kk] <- Bi1   
    Bur2[, kk] <- Bi2
  }
  
  #���ꂼ��̖ޓx�����l�ϐ��̕��S��
  y1 <- r[1]*rowSums(Bur1)
  y2 <- r[2]*rowSums(Bur2)
  y_rate <- y1/(y1 + y2)   #���S��
  r1 <- c(mean(y_rate), 1-mean(y_rate))
  
  #�g�s�b�N���z�̏d�ݕt���ޓx
  Bur <- y_rate*Bur1 + (1-y_rate)*Bur2   #�ϑ��f�[�^�̏d�ݕt���ޓx
  Br <- Bur / rowSums(Bur)   #���S���̌v�Z
  r2 <- colSums(Br) / sum(Br)   #�������̌v�Z
  bval <- list(Br=Br, Bur=Bur, y_rate=y_rate, r1=r1, r2=r2)
  return(bval)
}

####EM�A���S���Y���̏����l��ݒ肷��####
##�����l�������_���ɐݒ�
#phi�̏����l
freq_v <- matrix(colSums(WX), nrow=k, ncol=v, byrow=T)   #�P��̏o����
rand_v <- matrix(trunc(rnorm(k*v, 0, (colSums(WX)/2))), nrow=k, ncol=v, byrow=T)   #�����_����
phi_r <- abs(freq_v + rand_v) / rowSums(abs(freq_v + rand_v))   #�g�s�b�N���Ƃ̏o�����������_���ɏ�����

#beta�̏����l
beta_r <- rdirichlet(hh, runif(k, 0.2, 3))   #�f�B���N�����z���珉���l��ݒ�

#gamma�̏����l
gamma_r <- rdirichlet(item, runif(k, 0.2, 3))   #�f�B�N�������z���珉���l��ݒ�

#�g�s�b�N�������̍������̏����l
r <- c(0.5, 0.5)

###�p�����[�^�̍X�V
##���S���̌v�Z
bfr <- burden_fr(beta_r, gamma_r, phi_r, r, ID_u, ID_i, wd, w, k)
Br <- bfr$Br   #���S��
r <- bfr$r1   #������
y_rate <- bfr$y_rate

#beta�̍X�V
usum <- (data.frame(id=ID_u, Br=y_rate*Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
beta_r <- usum / matrix(rowSums(usum), nrow=hh, ncol=k)   #�p�����[�^���v�Z

#gamma�̍X�V
isum <- (data.frame(id=ID_i, Br=(1-y_rate)*Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
gamma_r <- isum / matrix(rowSums(isum), nrow=item, ncol=k)   #�p�����[�^���v�Z

#phi�̍X�V
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


###�p�����[�^�̍X�V
##���S���̌v�Z
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  ##���S���̌v�Z
  bfr <- burden_fr(beta_r, gamma_r, phi_r, r, ID_u, ID_i, wd, w, k)
  Br <- bfr$Br   #���S��
  r <- bfr$r1   #������
  y_rate <- bfr$y_rate
  
  #beta�̍X�V
  usum <- (data.frame(id=ID_u, Br=y_rate*Br) %>%
             dplyr::group_by(id) %>%
             dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
  beta_r <- usum / matrix(rowSums(usum), nrow=hh, ncol=k)   #�p�����[�^���v�Z
  
  #gamma�̍X�V
  isum <- (data.frame(id=ID_i, Br=(1-y_rate)*Br) %>%
             dplyr::group_by(id) %>%
             dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
  gamma_r <- isum / matrix(rowSums(isum), nrow=item, ncol=k)   #�p�����[�^���v�Z
  
  #phi�̍X�V
  vf <- (data.frame(id=wd, Br=Br) %>%
           dplyr::group_by(id) %>%
           dplyr::summarize_each(funs(sum)))[, 2:(k+1)]
  phi_r <- t(vf) / matrix(colSums(vf), nrow=k, ncol=v)
  
  #�ΐ��ޓx�̌v�Z
  (LLS <- sum(log(rowSums(bfr$Bur))))
  
  iter <- iter+1
  dl <- LLS-LLo
  LLo <- LLS
  LLw <- c(LLw, LLo)
  print(LLo)
}

(b <- round(cbind(beta_r, beta), 3))
(g <- round(cbind(gamma_r, gamma), 3))
(p <- round(cbind(t(phi_r), t(phi)), 3))