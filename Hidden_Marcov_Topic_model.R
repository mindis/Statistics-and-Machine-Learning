#####Hidden Marcov Topic Model#####
options(warn=2)
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

#set.seed(5723)

####�f�[�^�̔���####
k1 <- 15   #HMM�̍�����
k2 <- 15   #���ʂ̃g�s�b�N��
d <- 2000   #������
v <- 500   #��b��
s <- rpois(d, 15)   #���͐�
s[s < 5] <- ceiling(runif(sum(s < 5), 5, 10))
a <- sum(s)   #�����͐�
w <- rpois(a, 12)   #���͂�����̒P�ꐔ
w[w < 5] <- ceiling(runif(sum(w < 5), 5, 10))
f <- sum(w)   #���P�ꐔ

#����ID�̐ݒ�
u_id <- rep(1:d, s)
t_id <- c()
for(i in 1:d){t_id <- c(t_id, 1:s[i])}
words <- as.numeric(tapply(w, u_id, sum))

#���͋�؂�̃x�N�g�����쐬
ID_d <- rep(1:d, words)
td_d <- c()
for(i in 1:d){
  td_d <- c(td_d, rep(1:s[i], w[u_id==i]))
}
nd_d <- rep(1:a, w)
x_vec <- rep(0, f)
x_vec[c(1, cumsum(w[-a])+1)] <- 1

#�C���f�b�N�X��ݒ�
s_list <- list()
vec_list <- list()
for(i in 1:a){
  if(i%%1000==0){
    print(i)
  }
  s_list[[i]] <- which(nd_d==i)
  vec_list[[i]] <- rep(1, length(s_list[[i]]))
}

##�p�����[�^�̐ݒ�
#�f�B���N�����z�̃p�����[�^
alpha01 <- rep(1, k1)
alpha02 <- matrix(0.2, nrow=k1, ncol=k1)
diag(alpha02) <- 2.25
alpha03 <- rep(0.15, k2)
alpha11 <- rep(0.1, v)


for(l in 1:100){
  print(l)
  #�p�����[�^�𐶐�
  theta1 <- thetat1 <- extraDistr::rdirichlet(1, alpha01)
  theta2 <- thetat2 <- extraDistr::rdirichlet(k1, alpha02)
  theta3 <- thetat3 <- extraDistr::rdirichlet(k1, alpha03)
  phi0 <- t(extraDistr::rdirichlet(v, rep(0.01, k2))) * 
                          (matrix(extraDistr::rdirichlet(1, rep(2.0, v)), nrow=k2, ncol=v, byrow=T))
  phi <- phit <- phi0 / rowSums(phi0)

  ##���f���ɂ��ƂÂ��P��𐶐�����
  WX <- matrix(0, nrow=a, ncol=v)
  Z1_list <- list()
  Z2_list <- list()
  wd_list <- list()
  
  for(i in 1:a){
    #���͂��Ƃ̃Z�O�����g�𐶐�
    if(t_id[i]==1){
      z1 <- rmnom(1, 1, theta1)
      Z1_list[[i]] <- as.numeric(z1 %*% 1:k1)
    } else {
      z1 <- rmnom(1, 1, theta2[Z1_list[[i-1]], ])
      Z1_list[[i]] <- as.numeric(z1 %*% 1:k1)
    }
    
    #�P�ꂲ�ƂɃg�s�b�N�ƒP��𐶐�
    z2 <- rmnom(w[i], 1, theta3[Z1_list[[i]], ])
    Z2_list[[i]] <- as.numeric(z2 %*% 1:k2)
    
    #�g�s�b�N�Ɋ�Â��P��𐶐�
    word <- rmnom(w[i], 1, phi[Z2_list[[i]], ])
    WX[i, ] <- colSums(word)
    wd_list[[i]] <- as.numeric(word %*% 1:v)
  }
  if(min(colSums(WX)) > 0) break
}
colSums(WX)


#���X�g��ϊ�
wd <- unlist(wd_list)
z1 <- unlist(Z1_list)
z2 <- unlist(Z2_list)
Data <- matrix(as.numeric(table(1:f, wd)), nrow=f, ncol=v) 
sparse_data <- as(Data, "CsparseMatrix")
rm(Data)


##�C���f�b�N�X���쐬
doc_list <- list()
td_list <- s_list
word_list <- list()
for(i in 1:d){doc_list[[i]] <- which(ID_d==i)}
for(i in 1:v){word_list[[i]] <- which(wd==i)}


####�}���R�t�A�������e�J�����@��MHMM�g�s�b�N���f���𐄒�####
##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #���S�W���̊i�[�p
  for(j in 1:k){
    #���S�W�����v�Z
    Bi <- rep(theta[, j], w) * phi[j, wd]   #�ޓx
    Bur[, j] <- Bi   
  }
  Br <- Bur / rowSums(Bur)   #���S���̌v�Z
  bval <- list(Br=Br, Bur=Bur)
  return(bval)
}

##�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z���v�Z���邽�߂̊֐�
LLobz <- function(Data, phi, r, const, hh, k){
  
  #�������z�̑ΐ��ޓx
  log_phi <- log(t(phi))
  LLi <- const + Data %*% log_phi
  
  #logsumexp�̖ޓx
  LLi_max <- matrix(apply(LLi, 1, max), nrow=hh, ncol=k)
  r_matrix <- matrix(r, nrow=hh, ncol=k, byrow=T)
  
  #�����m���̃p�����[�^��ݒ�
  expl <- exp(LLi - LLi_max)
  z <- expl / rowSums(expl)   #�Z�O�����g�����m��
  
  #�ϑ��f�[�^�̑ΐ��ޓx
  r_log <- matrix(log(r), nrow=hh, ncol=k, byrow=T)
  LLosum <- sum(log(rowSums(exp(r_log + LLi))))   #�ϑ��f�[�^�̑ΐ��ޓx
  rval <- list(LLob=LLosum, z=z, LL=LLi)
  return(rval)
}


####HMM�g�s�b�N���f����MCMC�A���S���Y���̐ݒ�####
##�A���S���Y���̐ݒ�
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##�p�����[�^�̐^�l
theta1 <- thetat1
theta2 <- thetat2
theta3 <- thetat3
phi <- phit
z1_vec <- z1


##MHMT���f���̏����l��ݒ�
##�����������z�ŃZ�O�����g������������
const <- lfactorial(w) - rowSums(lfactorial(WX))   #�������z�̖��x�֐��̑ΐ��ޓx�̒萔

#�p�����[�^�̏����l
#phi�̏����l
alpha0 <- colSums(WX) / sum(WX) + 0.001
phi <- extraDistr::rdirichlet(k1, alpha0*v)

#�������̏����l
r <- rep(1/k1, k1)

#�ϑ��f�[�^�̑ΐ��ޓx�̏�����
L <- LLobz(WX, phi, r, const, a, k1)
LL1 <- L$LLob
z <- L$z

#�X�V�X�e�[�^�X
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 1
iter <- 0 

##EM�A���S���Y���őΐ��ޓx���ő剻
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  #E�X�e�b�v�̌v�Z
  z <- L$z   #���ݕϐ�z�̏o��
  
  #M�X�e�b�v�̌v�Z�ƍœK��
  #phi�̐���
  df0 <- matrix(0, nrow=k1, ncol=v)
  for(j in 1:k1){
    #���S�f�[�^�̑ΐ��ޓx����phi�̐���ʂ��v�Z
    phi[j, ] <- colSums(matrix(z[, j], nrow=a, ncol=v) * WX) / sum(z[, j] * w)   #�d�ݕt���������z�̍Ŗސ���
  }
  
  #�������𐄒�
  r <- apply(z, 2, sum) / a
  
  #�ϑ��f�[�^�̑ΐ��ޓx���v�Z
  phi[phi==0] <- min(phi[phi > 0])
  L <- LLobz(WX, phi, r, const, a, k1)
  LL <- L$LLob   #�ϑ��f�[�^�̑ΐ��ޓx
  iter <- iter+1   
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

#�����l��ݒ�
theta1 <- extraDistr::rdirichlet(1, rep(1, k1))
alpha <- matrix(0.3, nrow=k1, ncol=k1)
diag(alpha) <- 1.5
theta2 <- extraDistr::rdirichlet(k1, alpha)
theta3 <- extraDistr::rdirichlet(k1, rep(0.3, k2))
z1_vec <- as.numeric(rmnom(a, 1, z) %*% 1:k1)


##���O���z�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z
alpha01 <- 0.1
alpha02 <- 0.1
beta01 <- 1
beta02 <- 1


##�p�����[�^�̊i�[�p�z��
THETA1 <- matrix(0, nrow=R/keep, ncol=k1)
THETA2 <- array(0, dim=c(k1, k1, R/keep))
THETA3 <- array(0, dim=c(k1, k2, R/keep))
PHI <- array(0, dim=c(k2, v, R/keep))
SEG1 <- matrix(0, nrow=a, ncol=k1)
SEG2 <- matrix(0, nrow=f, ncol=k2)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"


##MCMC����p�z��
max_time <- max(t_id)
max_word <- max(words)
index_t11 <- which(t_id==1)
index_t21 <- list()
index_t22 <- list()
for(j in 2:max_word){
  index_t21[[j]] <- which(t_id==j)-1
  index_t22[[j]] <- which(t_id==j)
}

#�ΐ��ޓx�̊�l
LLst <- sum(WX %*% log(colSums(WX)/f))


####�M�u�X�T���v�����O��HTM���f���̃p�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�P�ꂲ�Ƃɕ��͒P�ʂŃg�s�b�N�𐶐�
  #�g�s�b�N�ޓx���v�Z
  z1_indicate <- rep(z1_vec, w)
  word_par <- matrix(0, nrow=f, ncol=k2)
  for(j in 1:k2){
    word_par[, j] <- theta3[z1_indicate, j] * phi[j, wd]   #�g�s�b�N�ޓx
  }
  

  #�g�s�b�N�̊����m������g�s�b�N�𐶐�
  topic_rate <- word_par / rowSums(word_par)
  Zi2 <- rmnom(f, 1, topic_rate)
  z2_vec <- as.numeric(Zi2 %*% 1:k2)


  ##HMM�ŕ��͒P�ʂ̃Z�O�����g�𐶐�
  #���͒P�ʂł̃g�s�b�N�p�x�s����쐬
  HMM_data <- matrix(0, nrow=a, ncol=k2)
  for(i in 1:a){
    HMM_data[i, ] <- vec_list[[i]] %*% Zi2[s_list[[i]], , drop=FALSE]
  }
  
  #���ݕϐ����Ƃɖޓx�𐄒�
  theta_log <- log(t(theta3))
  LLi0 <- HMM_data %*% theta_log   #�ΐ��ޓx
  LLi_max <- rowMaxs(LLi0)
  LLi <- exp(LLi0 - LLi_max)   #�ޓx�ɕϊ�
  
  #�Z�O�����g�����m���̐���ƃZ�O�����g�̐���
  z_rate1 <- matrix(0, nrow=a, ncol=k1)
  Zi1 <- matrix(0, nrow=a, ncol=k1)
  z1_vec <- rep(0, a)
  rf02 <- matrix(0, nrow=k1, ncol=k1) 

  for(j in 1:max_time){
    if(j==1){
      #�Z�O�����g�̊����m��
      LLs <- matrix(theta1, nrow=length(index_t11), ncol=k1, byrow=T) * LLi[index_t11, ]   #�d�ݕt���ޓx
      z_rate1[index_t11, ] <- LLs / rowSums(LLs)   #�����m��
      
      #�������z���Z�O�����g�𐶐�
      Zi1[index_t11, ] <- rmnom(length(index_t11), 1, z_rate1[index_t11, ])
      z1_vec[index_t11] <- as.numeric(Zi1[index_t11, ] %*% 1:k1)
      
      #�������̃p�����[�^���X�V
      rf01 <- colSums(Zi1[index_t11, ])
      
    } else {
      
      #�Z�O�����g�̊����m��
      index <- index_t22[[j]]
      LLs <- theta2[z1_vec[index_t21[[j]]], , drop=FALSE] * LLi[index, , drop=FALSE]   #�d�ݕt���ޓx
      z_rate1[index, ] <- LLs / rowSums(LLs)   #�����m��
      
      #�������z���Z�O�����g�𐶐�
      Zi1[index, ] <- rmnom(length(index), 1, z_rate1[index, ])
      z1_vec[index] <- as.numeric(Zi1[index, ] %*% 1:k1)
      
      #�������̃p�����[�^���X�V
      rf02 <- rf02 + t(Zi1[index_t21[[j]], , drop=FALSE]) %*% Zi1[index, , drop=FALSE]   #�}���R�t����
    }
  }
  
  ##�p�����[�^���T���v�����O
  #�f�B�N�������z����HMM�̍��������T���v�����O
  rf11 <- colSums(Zi1[index_t11, ]) + beta01
  rf12 <- rf02 + alpha01
  theta1 <- extraDistr::rdirichlet(1, rf11)
  theta2 <- extraDistr::rdirichlet(k1, rf12)

  #�g�s�b�N���z�̃p�����[�^���T���v�����O
  wf0 <- matrix(0, nrow=k1, ncol=k2)
  for(j in 1:k1){
    wf0[j, ] <- colSums(HMM_data * Zi1[, j])
  }
  wf <- wf0 + beta02
  theta3 <- extraDistr::rdirichlet(k1, wf)
  
  #�P�ꕪ�zphi���T���v�����O
  vf0 <- matrix(0, nrow=k2, ncol=v)
  for(j in 1:v){
    vf0[, j] <- colSums(Zi2[word_list[[j]], , drop=FALSE])
  }
  vf <- vf0 + alpha02
  phi <- extraDistr::rdirichlet(k2, vf)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA1[mkeep, ] <- theta1
    THETA2[, , mkeep] <- theta2
    THETA3[, , mkeep] <- theta3
    PHI[, , mkeep] <- phi
    
    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(mkeep >= burnin & rp%%keep==0){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
    }
    
    #�T���v�����O���ʂ��m�F
    if(rp%%disp==0){
      print(rp)
      print(c(sum(log(rowSums(word_par))), LLst))
      print(round(rbind(theta1, thetat1), 3))
      print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
    }
  }
}


####�T���v�����O���ʂ̉����Ɨv��####
burnin <- 2000/keep   #�o�[���C������
RS <- R/keep

##�T���v�����O���ʂ̉���
#HMM�̏������z�̃T���v�����O����
matplot(THETA1, type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

#HMM�̃p�����[�^�̃T���v�����O����
matplot(t(THETA2[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA2[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA2[15, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

#�����̃g�s�b�N���z�̃T���v�����O����
matplot(t(THETA3[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA3[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA3[15, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

#�P��̏o���m���̃T���v�����O����
matplot(t(PHI[, 1, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N1�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[, 100, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N2�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[, 200, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N2�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[, 300, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N2�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[, 400, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N3�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[, 500, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N4�̒P��̏o�����̃T���v�����O����")


##�T���v�����O���ʂ̗v�񐄒��
#�g�s�b�N���z�̎��㐄���
topic_mu <- apply(THETA[, , burnin:(R/keep)], c(1, 2), mean)   #�g�s�b�N���z�̎��㕽��
round(cbind(topic_mu, thetat), 3)
round(topic_sd <- apply(THETA[, , burnin:(R/keep)], c(1, 2), sd), 3)   #�g�s�b�N���z�̎���W���΍�

#�P��o���m���̎��㐄���
word_mu <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #�P��̏o�����̎��㕽��
word <- round(t(rbind(word_mu, phit)), 3)
colnames(word) <- 1:ncol(word)
word

##�g�s�b�N�̎��㕪�z�̗v��
round(cbind(z1, seg1_mu <- SEG1 / length(burnin:RS)), 3)
round(cbind(z2, seg2_mu <- SEG2 / rowSums(SEG2)), 3)



