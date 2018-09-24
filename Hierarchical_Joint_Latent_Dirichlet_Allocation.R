#####Hierarchical Joint Latent Dirichlet Allocation#####
options(warn=0)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(2506787)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
k1 <- 20   #���[�U�[�P�ʂ̃g�s�b�N��
k2 <- 30   #�Z�b�V�����P�ʂ̃g�s�b�N��
k3 <- 30   #�A�C�e���P�ʂ̃g�s�b�N��
hh <- 3000   #���[�U�[��
item <- 1500   #�A�C�e����
pages <- 25   #�y�[�W��
pt <- rtpois(hh, rgamma(hh, 20, 0.25), a=1, b=Inf)   #�Z�b�V������
hhpt <- sum(pt)   #���Z�b�V������
w <- extraDistr::rtpois(hhpt, rgamma(hhpt, 2.75, 0.5), a=0, b=Inf)   #�Z�b�V����������̉{���A�C�e����
f <- sum(w)   #���{���A�C�e����

##ID�ƃC���f�b�N�X�̐ݒ�
#�Z�b�V�����P�ʂ�ID��ݒ�
user_id <- rep(1:hh, pt)

#�A�C�e���P�ʂ�ID�̐ݒ�
user_id_vec <- rep(rep(1:hh, pt), w)
session_id <- rep(1:hhpt, w)
user_no <- as.numeric(unlist(tapply(1:f, user_id_vec, rank)))
session_no <- as.numeric(unlist(tapply(1:f, session_id, rank)))

#�C���f�b�N�X�̐ݒ�
user_index <- session_index <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id_vec==i)
}
vec <- c(0, cumsum(w))
for(i in 1:hhpt){
  session_index[[i]] <- (1:f)[(vec[i]+1):vec[i+1]]
}

##�p�����[�^�𐶐�
#�f�B���N�����z�̃p�����[�^��ݒ�
alpha01 <- rep(0.1, k1)
alpha02 <- rep(0.1, k3)
beta01 <- rep(0.1, k2)
beta02 <- rep(0.2, pages)
beta03 <- matrix(0.015, nrow=k3, ncol=item, byrow=T)
for(j in 1:k2){
  beta03[j, matrix(1:item, nrow=k3, ncol=item/k3, byrow=T)[j, ]] <- 0.1
}


##���ׂẴA�C�e�����o������܂Ńf�[�^�̐����𑱂���
rp <- 0
repeat {
  rp <- rp + 1
  
  #�f�B���N�����z����p�����[�^�𐶐�
  theta1 <- thetat1 <- extraDistr::rdirichlet(hh, alpha01)
  theta2 <- thetat2 <- extraDistr::rdirichlet(k2, alpha02)
  gamma <- gammat <- extraDistr::rdirichlet(k1, beta01)
  omega <- omegat <- extraDistr::rdirichlet(k1, beta02)
  phi <- phit <- extraDistr::rdirichlet(k3, beta03)

  #�A�C�e���o���m�����Ⴂ�g�s�b�N�����ւ���
  index <- which(colMaxs(phi) < (k2*5)/f)
  for(j in 1:length(index)){
    phi[as.numeric(rmnom(1, 1, extraDistr::rdirichlet(1, beta01)) %*% 1:k2), index[j]] <- (k2*5)/f
  }
  
  ##�����ϐ��𐶐�
  y_list <- v_list <- p_list <- z_list <- w_list <- list()
  w_sums <- rep(0, item)
  
  for(i in 1:hh){
    #�C���f�b�N�X�𒊏o
    u_index <- user_index[[i]]
    
    #�Z�b�V�����̃g�s�b�N�𐶐�
    y <- rmnom(pt[i], 1, theta1[i, ])
    y_vec <- as.numeric(y %*% 1:k1)
    
    #�Z�b�V�����g�s�b�N����Z�b�V�����𐶐�
    v <- rmnom(pt[i], 1, gamma[y_vec, ])
    v_vec <- as.numeric(v %*% 1:k2)
    
    #�Z�b�V��������A�C�e���g�s�b�N�ƃy�[�W�𐶐�
    s_index <- session_id[u_index] - (min(session_id[u_index])-1); n <- length(s_index)
    p <- rmnom(n, 1, omega[y_vec[s_index], ])
    p_vec <- as.numeric(p %*% 1:pages)
    z <- rmnom(n, 1, theta2[v_vec[s_index], ])
    z_vec <- as.numeric(z %*% 1:k2)
    
    #�A�C�e���g�s�b�N����A�C�e���𐶐�
    w <- rmnom(n, 1, phi[z_vec, ])
    w_vec <- as.numeric(w %*% 1:item)
  
    #���������f�[�^���i�[
    y_list[[i]] <- y
    v_list[[i]] <- v
    p_list[[i]] <- p_vec
    z_list[[i]] <- z_vec
    w_list[[i]] <- w_vec
    w_sums <- w_sums + colSums(w)
  }
  #break����
  print(c(rp, sum(w_sums==0)))
  if(sum(w_sums==0)==0){
    break
  }
}

#���X�g��ϊ�
y <- do.call(rbind, y_list); y_vec <- (y %*% 1:k1)
v <- do.call(rbind, v_list); v_vec <- (v %*% 1:k2)[session_id]
z_vec <- unlist(z_list)
pd <- unlist(p_list)
wd <- unlist(w_list)
page_data <- sparseMatrix(1:f, pd, x=rep(1, f), dims=c(f, pages))
page_data_T <- t(page_data)
item_data <- sparseMatrix(1:f, wd, x=rep(1, f), dims=c(f, item))
item_data_T <- t(item_data)
hist(w_sums, col="grey", breaks=50, xlab="�A�C�e���o���p�x", main="�A�C�e���o�����z")


####�M�u�X�T���v�����O��Hierarchical Joint Latent Dirichlet Allocation�𐄒�####
##�A�C�e�����Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, w, k){
  #���S�W�����v�Z
  Bur <- theta[w, ] * t(phi)[wd, ]   #�ޓx
  Br <- Bur / rowSums(Bur)   #���S��
  r <- colSums(Br) / sum(Br)   #������
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##�A���S���Y���̐ݒ�
R <- 2000
keep <- 2  
iter <- 0
burnin <- 300
disp <- 10

##���O���z�̐ݒ�
#�f�B���N�����z�̎��O���z
alpha01 <- 0.1
beta01 <- 0.1
er <- 0.0001

##�p�����[�^�̐^�l
theta1 <- thetat1
theta2 <- thetat2
gamma <- gammat
omega <- omegat
phi <- phit

##�p�����[�^�̏����l��ݒ�
#�g�s�b�N���z�̏����l
theta1 <- extraDistr::rdirichlet(hh, rep(2.0, k1))
theta2 <- extraDistr::rdirichlet(k2, rep(2.0, k3))

#�A�C�e�����z�̏����l
gamma <- extraDistr::rdirichlet(k1, rep(2.0, k2))
omega <- extraDistr::rdirichlet(k1, rep(2.0, pages))
phi <- extraDistr::rdirichlet(k3, rep(2.0, item))

#�g�s�b�N�̏����l
y_vec <- as.numeric(rmnom(hhpt, 1, rep(1, k1)) %*% 1:k1)


##�p�����[�^�̊i�[�p�z��
THETA1 <- array(0, dim=c(hh, k1, R/keep))
THETA2 <- array(0, dim=c(k2, k3, R/keep))
GAMMA <- array(0, dim=c(k1, k2, R/keep))
OMEGA <- array(0, dim=c(k1, pages, R/keep))
PHI <- array(0, dim=c(k3, item, R/keep))
SEG1 <- matrix(0, nrow=hhpt, ncol=k1)
SEG2 <- matrix(0, nrow=hhpt, ncol=k2)
SEG3 <- matrix(0, nrow=f, ncol=k3)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"
storage.mode(SEG3) <- "integer"

##�f�[�^�ƃC���f�b�N�X�̐ݒ�
#�C���f�b�N�X�̐ݒ�
user_dt <- sparseMatrix(user_id, 1:hhpt, x=rep(1, hhpt), dims=c(hh, hhpt))
session_dt <- sparseMatrix(session_id, 1:f, x=rep(1, f), dims=c(hhpt, f))


##�ΐ��ޓx�̊�l
LLst <- sum(item_data %*% log(colMeans(item_data)))
LLbest <- sum(log(as.numeric((thetat2[as.numeric(v %*% 1:k2)[session_id], ] * t(phit)[wd, ]) %*% rep(1, k3))))


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�Z�b�V�����̐��ݕϐ��𐶐�
  #�f�[�^�̐ݒ�
  theta1_vec <- theta1[user_id, ]   #�A�C�e���g�s�b�N
  phi_vec <- t(phi)[wd, ]   #�P�ꕪ�z�@
  
  #�Z�b�V�����̊����m���𐶐�
  Lho_topic <- exp(session_dt %*% log(phi_vec %*% t(theta2)))
  Lho_session <- gamma[y_vec, ] * Lho_topic
  session_rate <- Lho_session / as.numeric(Lho_session %*% rep(1, k2))
  session_rate <- (session_rate + er) / (session_rate + er) %*% rep(1, k2)
  
  #�������z����Z�b�V�����𐶐�
  V <- extraDistr::rmnom(hhpt, 1, session_rate)
  v_vec <- as.numeric(V %*% 1:k2)
  
  ##�A�C�e���g�s�b�N�𐶐�
  #�A�C�e���g�s�b�N�̊����m��
  v_session <- v_vec[session_id]
  Lho_w <- theta2[v_session, ] * phi_vec   #�A�C�e���g�s�b�N�̊��Җޓx
  topic_rate <- Lho_w / as.numeric(Lho_w %*% rep(1, k2))
  
  #�������z����A�C�e���g�s�b�N�𐶐�
  Z <- extraDistr::rmnom(f, 1, topic_rate)
  z_vec <- as.numeric(Z %*% 1:k2)
  
  ##�Z�b�V�����g�s�b�N�𐶐�
  #�Z�b�V�����g�s�b�N�̊����m��
  Lho_s <- theta1_vec * exp(session_dt %*% t(log(omega))[pd, ]) * t(gamma)[v_vec, ]   #�Z�b�V�����g�s�b�N�̊��Җޓx
  topic_rate <- Lho_s / as.numeric(Lho_s %*% rep(1, k1))
  topic_rate <- (topic_rate + er) / (topic_rate + er) %*% rep(1, k1)
  
  #�������z����Z�b�V�����g�s�b�N�𐶐�
  Y <- extraDistr::rmnom(hhpt, 1, topic_rate)
  y_vec <- as.numeric(Y %*% 1:k1)
  
  ##�g�s�b�N���z�̃p�����[�^���T���v�����O
  #�f�B���N�����z�̃p�����[�^��ݒ�
  v_dt <- sparseMatrix(v_session, 1:f, x=rep(1, f), dims=c(k2, f)) 
  y_sums <- user_dt %*% Y + alpha01
  z_sums <- v_dt %*% Z + alpha01
  
  #�f�B���N�����z����g�s�b�N���z�𐶐�
  theta1 <- extraDistr::rdirichlet(hh, y_sums)
  theta2 <- extraDistr::rdirichlet(k2, z_sums)
  
  
  ##�o���m�����z�̃p�����[�^���T���v�����O
  #�f�B���N�����z�̃p�����[�^��ݒ�
  y_dt <- sparseMatrix(1:f, y_vec[session_id], x=rep(1, f), dims=c(f, k1))
  p_sums <- t(page_data_T %*% y_dt) + beta01
  s_sums <- t(v_dt %*% y_dt) + beta01
  w_sums <- t(item_data_T %*% Z) + beta01
    
  #�f�B���N�����z����o���m�����z���T���v�����O
  omega <- extraDistr::rdirichlet(k1, p_sums)
  gamma <- extraDistr::rdirichlet(k1, s_sums)
  phi <- extraDistr::rdirichlet(k3, w_sums)

  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA1[, , mkeep] <- theta1
    THETA2[, , mkeep] <- theta2
    PHI[, , mkeep] <- phi
    GAMMA[, , mkeep] <- gamma
    OMEGA[, , mkeep] <- omega
  }  
  
  #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
  if(rp%%keep==0 & rp >= burnin){
    SEG1 <- SEG1 + Y
    SEG2 <- SEG2 + V
    SEG3 <- SEG3 + Z
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    LL <- sum(log(as.numeric((theta2[v_vec[session_id], ] * t(phi)[wd, ]) %*% rep(1, k3))))
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(c(LL, LLbest, LLst))
  }
}
