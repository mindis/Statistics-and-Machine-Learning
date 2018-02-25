#####�A����ԃg�s�b�N���f��#####
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
u <- ut <- mvrnorm(d, rep(0, k), diag(k))
phi <- phit <- mvrnorm(v, rep(0, k), diag(k))

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
word_vec <- unlist(word_vec_list)
word_data <- do.call(rbind, word_list)
sparse_data <- as(word_data, "CsparseMatrix")
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
R <- 2000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 4

#�f�[�^�̐ݒ�
rej2 <- rep(0, v)
WX_vec <- as.numeric(WX)

#�C���f�b�N�X��ݒ�
index_nzeros <- which(as.numeric(WX) > 0)
index_word <- list()
for(j in 1:v){
  index_word[[j]] <- which(WX[, j] > 0)
}

##���O���z�̐ݒ�
cov <- diag(k)
inv_cov <- solve(cov)
mu <- rep(0, k)
sigma <- diag(k)


##�p�����[�^�̐^�l
G0 <- GT0
u <- ut
phi <- phit


##�p�����[�^�̏����l
G0 <- colSums(WX) / sum(WX)
u <- mvrnorm(d, rep(0, k), diag(k))
phi <- mvrnorm(v, rep(0, k), diag(k))


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
  
  ##�g�s�b�N���z�̃p�����[�^u���T���v�����O
  #�V�����p�����[�^���T���v�����O
  u_old <- u
  u_new <- u_old + mvrnorm(d, rep(0, k), diag(0.01, k))
  alphan <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u_new %*% t(phi))
  alphad <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u_old %*% t(phi))
  
  #Polya���z�̃p�����[�^���v�Z
  dirn_vec <- dird_vec <- rep(0, d*v)
  alphan_vec <- as.numeric(alphan)
  alphad_vec <- as.numeric(alphad)
  dirn_vec[index_nzeros] <- lgamma(alphan_vec[index_nzeros] + WX_vec[index_nzeros]) - lgamma(alphan_vec[index_nzeros])
  dird_vec[index_nzeros] <- lgamma(alphad_vec[index_nzeros] + WX_vec[index_nzeros]) - lgamma(alphad_vec[index_nzeros])
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew1 <- lgamma(rowSums(alphan)) - lgamma(rowSums(alphan + WX)) + rowSums(matrix(dirn_vec, nrow=d, ncol=v))
  logold1 <- lgamma(rowSums(alphad)) - lgamma(rowSums(alphad + WX)) + rowSums(matrix(dird_vec, nrow=d, ncol=v))
  logpnew1 <- -0.5 * rowSums(u_new %*% inv_cov * u_new)
  logpold1 <- -0.5 * rowSums(u_old %*% inv_cov * u_old)
  
  
  #MH�T���v�����O
  rand <- runif(d)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew1 + logpnew1 - logold1 - logpold1)   #�̑𗦂��v�Z
  tau <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #tau�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((tau >= rand)*1 + (tau < rand)*0), nrow=d, ncol=k)
  rej1 <- mean(flag[, 1])
  u <- flag*u_new + (1-flag)*u_old   #alpha��rand�������Ă�����̑�
  u <- scale(u)   #���K��
  
  
  ##�P�ꕪ�z�̃p�����[�^phi���T���v�����O
  for(j in 1:v){
    index_vec <- index_word[[j]]
    
    if(j==1){
      #�V�����p�����[�^���T���v�����O
      phid <- phi 
      alphan <- alphad <- matrix(G0, nrow=d, ncol=v, byrow=T) * exp(u %*% t(phid))
      phin <- phid[j, ] + mvrnorm(1, rep(0, k), diag(0.01, k))
      alphan[, j] <- G0[j] * exp(u %*% phin)
      
      #�J��Ԃ�1��ڂ̃p�����[�^���X�V
      dir_vec <- rep(0, d*v)
      alpha_vec <- as.numeric(alphad)
      dir_vec[index_nzeros] <- lgamma(alpha_vec[index_nzeros] + WX_vec[index_nzeros]) - lgamma(alpha_vec[index_nzeros])
      
      #Polya���z�̃p�����[�^���X�V
      dir_mnd1 <- lgamma(rowSums(alphad)) - lgamma(rowSums(alphad + WX))
      dir_mnn2 <- dir_mnd2 <- matrix(dir_vec, nrow=d, ncol=v)
      dir_mnn1 <- lgamma(rowSums(alphan)) - lgamma(rowSums(alphan + WX))
      dir_mnn2[index_vec, j] <- lgamma(alphan[index_vec, j] + WX[index_vec, j]) - lgamma(alphan[index_vec, j])
    
    } else {
        
      #�V�����p�����[�^���T���v�����O
      alphan <- alphad
      phin <- phi[j, ] + mvrnorm(1, rep(0, k), diag(0.01, k))
      alphan[, j] <- G0[j] * exp(u %*% phin)
      
      #�J��Ԃ�2��ڈȍ~�̃p�����[�^���X�V
      dir_mnn2 <- dir_mnd2
      dir_mnn1 <- lgamma(rowSums(alphan)) - lgamma(rowSums(alphan + WX))
      dir_mnn2[index_vec, j] <- lgamma(alphan[index_vec, j] + WX[index_vec, j]) - lgamma(alphan[index_vec, j])
    }
    
    #�ΐ��ޓx�Ƒΐ����O���z���v�Z
    lognew2 <- sum(dir_mnn1) + sum(dir_mnn2)
    logold2 <- sum(dir_mnd1) + sum(dir_mnd2)
    logpnew2 <- lndMvn(phin, mu, cov)
    logpold2 <- lndMvn(phid[j, ], mu, cov)
    
    #MH�T���v�����O
    rand <- runif(1)   #��l���z���痐���𔭐�
    tau <- min(c(1, exp(lognew2 + logpnew2 - logold2 - logpold2)))   #�̑𗦂��v�Z
    
    #tau�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
    rej2[j] <- flag <- as.numeric(tau >= rand)
    phi[j, ] <- flag*phin + (1-flag)*phid[j, ]
    alphad[, j] <- flag*alphan[, j] + (1-flag)*alphad[, j]
    dir_mnd1 <- flag*dir_mnn1 + (1-flag)*dir_mnd1
    dir_mnd2 <- flag*dir_mnn2 + (1-flag)*dir_mnd2
  }
  phi <- scale(phi)   #���K��
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    logl[mkeep] <- lognew2
    U[, , mkeep] <- u
    PHI[, , mkeep] <- phi
    
    #�T���v�����O���ʂ��m�F
    if(rp%%disp==0){
      print(rp)
      print(c(lognew2, LLbest, LLst))
      print(round(c(rej1, mean(rej2)), 3))
      print(round(cbind(u[1:5, ], ut[1:5, ]), 2))
      print(round(cbind(phi[1:5, ], phit[1:5, ]), 2))
    }
  }
}

