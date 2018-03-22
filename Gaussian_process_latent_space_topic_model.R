#####�K�E�X�ߒ��A����ԃg�s�b�N���f��#####
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
gc()
#set.seed(5723)

####�f�[�^�̔���####
k <- 10   #�g�s�b�N��
d <- 2000   #���[�U�[��
v <- 500   #��b��
w <- rpois(d, rgamma(d, 45, 0.3))   #1�l������̃y�[�W�{����
f <- sum(w)   #����b��

#ID��ݒ�
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}


##�p�����[�^�̐ݒ�
G0 <- GT0 <- extraDistr::rdirichlet(1, rep(2.0, v))   #�P��o����
u <- ut <- mvrnorm(d, rep(0, k), diag(0.75, k))
phi <- phit <- mvrnorm(v, rep(0, k), diag(0.75, k))


##�f�[�^�̐���
word_list <- list()
word_vec_list <- list()
WX <- matrix(0, nrow=d, ncol=v)

for(i in 1:d){

  #�f�B�N����-�������z����P��𐶐�
  alpha <- G0 * exp(u[i, ] %*% t(phi))
  words <- extraDistr::rdirmnom(w[i], 1, alpha)
  words_vec <- as.numeric(words %*% 1:v)
  
  #���������P����i�[
  WX[i, ] <- colSums(words)
  word_list[[i]] <- words
  word_vec_list[[i]] <- words_vec
}

#���X�g��ϊ�
WX_T <- t(WX)
word_vec <- unlist(word_vec_list)
word_data <- do.call(rbind, word_list)
sparse_data <- as(word_data, "CsparseMatrix")
storage.mode(WX) <- "integer"
storage.mode(WX) <- "integer"
rm(word_data); rm(word_list)


##�C���f�b�N�X���쐬
dw_list <- list()
for(j in 1:v){
  index <- which(word_vec==j)
  dw_list[[j]] <- d_id[index]
}


####�}���R�t�A�������e�J�����@�ŘA����ԃg�s�b�N���f���𐄒�####
##�A���S���Y���̐ݒ�
R <- 10000
keep <- 5 
iter <- 0
burnin <- 2500/keep
disp <- 10

#�f�[�^�̐ݒ�
rej2 <- rep(0, v)
WX_vec <- as.numeric(WX)
WX_T_vec  <- as.numeric(t(WX))

#�C���f�b�N�X��ݒ�
index_nzeros1 <- which(as.numeric(WX) > 0)
index_nzeros2 <- which(as.numeric(WX_T) > 0)
index_word <- list()
for(j in 1:v){
  index_word[[j]] <- which(WX[, j] > 0)
}

##���O���z�̐ݒ�
cov <- diag(k)
inv_cov2 <- inv_cov1 <- solve(cov)
mu <- rep(0, k)
sigma <- diag(k)


##�p�����[�^�̐^�l
G0 <- GT0
u <- ut
phi <- phit

##�p�����[�^�̏����l
G0 <- colSums(WX) / sum(WX)
K0 <- rowSums(WX) / sum(WX)
u <- mvrnorm(d, rep(0, k), diag(0.01, k))
phi <- mvrnorm(v, rep(0, k), diag(0.01, k))


##�p�����[�^�̕ۑ��p�z��
logl <- rep(0, R/keep)
U <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(v, k, R/keep))

##�ΐ��ޓx�̊�l
#���j�O�������f���̑ΐ��ޓx
par <- colSums(WX)/sum(WX)
LLst <- sum(WX %*% log(par))

#�x�X�g�ȑΐ��ޓx
alpha <- matrix(GT0, nrow=d, ncol=v, byrow=T) * exp(ut %*% t(phit))
LLbest <- sum(lgamma(rowSums(alpha)) - lgamma(rowSums(alpha + WX))) + sum(lgamma(alpha + WX) - lgamma(alpha))



####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##���[�U�[�g�s�b�N�s��U���T���v�����O
  #�V�����p�����[�^���T���v�����O
  u_old <- u
  u_new <- u_old + mvrnorm(d, rep(0, k), diag(0.01, k))
  alphan <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u_new %*% t(phi))
  alphad <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u_old %*% t(phi))
  
  #Polya���z�̃p�����[�^���v�Z
  dirn_vec <- dird_vec <- rep(0, d*v)
  alphan_vec <- as.numeric(alphan)
  alphad_vec <- as.numeric(alphad)
  dirn_vec[index_nzeros1] <- lgamma(alphan_vec[index_nzeros1] + WX_vec[index_nzeros1]) - lgamma(alphan_vec[index_nzeros1])
  dird_vec[index_nzeros1] <- lgamma(alphad_vec[index_nzeros1] + WX_vec[index_nzeros1]) - lgamma(alphad_vec[index_nzeros1])
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  dir_new <- lgamma(rowSums(alphan)) - lgamma(rowSums(alphan + WX))
  dir_old <- lgamma(rowSums(alphad)) - lgamma(rowSums(alphad + WX))
  lognew1 <- dir_new + rowSums(matrix(dirn_vec, nrow=d, ncol=v))
  logold1 <- dir_old + rowSums(matrix(dird_vec, nrow=d, ncol=v))
  logpnew1 <- -0.5 * rowSums(u_new %*% inv_cov1 * u_new)
  logpold1 <- -0.5 * rowSums(u_old %*% inv_cov1 * u_old)
  
  ##MH�T���v�����O
  rand <- runif(d)   #��l���z���痐���𔭐�
  LL_diff <- lognew1 + logpnew1 - logold1 - logpold1
  LL_diff[LL_diff > 100] <- 100
  LLind_diff <- exp(LL_diff)   #�̑𗦂��v�Z
  tau <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #tau�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((tau >= rand)*1 + (tau < rand)*0), nrow=d, ncol=k)
  rej1 <- mean(flag[, 1])
  u <- flag*u_new + (1-flag)*u_old   #alpha��rand�������Ă�����̑�


  ##�P�ꕪ�z�̃p�����[�^phi���T���v�����O
  index_vec <- index_word[[j]]
  
  #�P�ꂲ�Ƃ�MH�T���v�����O�����s
  #�V�����p�����[�^���T���v�����O
  phid <- phi
  phin <- phid + mvrnorm(v, rep(0, k), diag(0.05, k))
  alphan <- matrix(K0, nrow=v, ncol=d) * exp(phin %*% t(u))
  alphad <- matrix(K0, nrow=v, ncol=d) * exp(phid %*% t(u))
  
  #Polya���z�̃p�����[�^���v�Z
  dirn_vec <- dird_vec <- rep(0, d*v)
  alphan_vec <- as.numeric(alphan)
  alphad_vec <- as.numeric(alphad)
  dirn_vec[index_nzeros2] <- lgamma(alphan_vec[index_nzeros2] + WX_T_vec[index_nzeros2]) - lgamma(alphan_vec[index_nzeros2])
  dird_vec[index_nzeros2] <- lgamma(alphad_vec[index_nzeros2] + WX_T_vec[index_nzeros2]) - lgamma(alphad_vec[index_nzeros2])
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  dir_new <- lgamma(rowSums(alphan)) - lgamma(rowSums(alphan + WX_T))
  dir_old <- lgamma(rowSums(alphad)) - lgamma(rowSums(alphad + WX_T))
  lognew2 <- dir_new + rowSums(matrix(dirn_vec, nrow=v, ncol=d))
  logold2 <- dir_old + rowSums(matrix(dird_vec, nrow=v, ncol=d))
  logpnew2 <- -0.5 * rowSums(phin %*% inv_cov2 * phin)
  logpold2 <- -0.5 * rowSums(phid %*% inv_cov2 * phid)
  
  
  ##MH�T���v�����O
  rand <- runif(v)   #��l���z���痐���𔭐�
  LL_diff <- lognew2 + logpnew2 - logold2 - logpold2
  LL_diff[LL_diff > 100] <- 100
  LLind_diff <- exp(LL_diff)   #�̑𗦂��v�Z
  tau <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #tau�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((tau >= rand)*1 + (tau < rand)*0), nrow=v, ncol=k)
  rej2 <- mean(flag[, 1])
  phi <- flag*phin + (1-flag)*phid   #alpha��rand�������Ă�����̑�

  #�p�����[�^�𐳋K��
  u <- scale(u)
  phi <- scale(phi)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    logl[mkeep] <- sum(lognew1)
    U[, , mkeep] <- u
    PHI[, , mkeep] <- phi
    
    #�T���v�����O���ʂ��m�F
    if(rp%%disp==0){
      alpha <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u %*% t(phi))
      LL <- sum(lgamma(rowSums(alpha)) - lgamma(rowSums(alpha + WX))) + sum(lgamma(alpha + WX) - lgamma(alpha))
   
      print(rp)
      print(c(LL, LLbest, LLst))
      print(round(c(rej1, mean(rej2)), 3))
      print(round(cbind(u[1:5, ], ut[1:5, ]), 2))
      print(round(cbind(phi[1:5, ], phit[1:5, ]), 2))
    }
  }
}

####�T���v�����O���ʂ̉����Ɨv��####
##�T���v�����O���ʂ��v���b�g
matplot(t(U[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�g�s�b�N���z�̃T���v�����O����")
matplot(t(U[2, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�g�s�b�N���z�̃T���v�����O����")
matplot(t(U[3, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�g�s�b�N���z�̃T���v�����O����")
matplot(t(U[4, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�g�s�b�N���z�̃T���v�����O����")
matplot(t(U[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�g�s�b�N���z�̃T���v�����O����")

matplot(t(PHI[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�g�s�b�N���z�̃T���v�����O����")
matplot(t(PHI[2, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�g�s�b�N���z�̃T���v�����O����")
matplot(t(PHI[3, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�g�s�b�N���z�̃T���v�����O����")
matplot(t(PHI[4, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�g�s�b�N���z�̃T���v�����O����")
matplot(t(PHI[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�g�s�b�N���z�̃T���v�����O����")
