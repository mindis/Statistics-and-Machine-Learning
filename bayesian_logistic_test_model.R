#####�x�C�W�A����l���W�X�e�B�b�N�e�X�g���f��####
library(irtoys)
library(bayesm)
library(extraDistr)
library(matrixStats)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(43587)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 5000   #�팱�Ґ�
k <- 100   #���ڐ�

##�p�����[�^�̐ݒ�
theta0 <- rnorm(hh, 0, 1)   #�팱�ҕꐔ
beta0 <- rnorm(k, 0.75, 1.25)   #����x�ꐔ
alpha0 <- runif(k, 0.3, 2.0)   #���ʗ͕ꐔ
c0 <- runif(k, 0.1, 0.3)   #���Đ��ʕꐔ 

##�����ϐ��̔���
Pr0 <- matrix(0, nrow=hh, ncol=k)
Data <- matrix(0, nrow=hh, ncol=k)

for(j in 1:k){
  Pr0[, j] <- c0[j] + (1-c0[j]) / (1+exp(-alpha0[j]*(theta0-beta0[j])))   #�������̐ݒ�
  Data[, j] <- rbinom(hh, 1, Pr0[, j])   #�����L���̐���
}
colMeans(Data)   #���ڂ��Ƃ̕��ϐ�����
mean(Data)   #�S���ڂł̐�����


####�}���R�t�A�������e�J�����@�œ�l���W�X�e�B�b�N�e�X�g���f���𐄒�####
##��l���W�X�e�B�b�N�e�X�g���f���̑ΐ��ޓx���`
loglike <- function(Data, theta, beta, alpha, c, hh, k){
  
  #�p�����[�^�̐ݒ�
  theta_m <- matrix(theta, nrow=hh, ncol=k)
  gamma_m <- matrix(c, nrow=hh, ncol=k, byrow=T)
  alpha_m <- matrix(alpha, nrow=hh, ncol=k, byrow=T)
  beta_m <- matrix(beta, nrow=hh, ncol=k, byrow=T)
  
  #3�p�����[�^���W�X�e�B�b�N�e�X�g���f���̔����m�����`
  pr <- gamma_m + (1-gamma_m) / (1+exp(-alpha_m*(theta_m-beta_m)))   
  
  #�ΐ��ޓx���`
  LLi <- Data*log(pr) + (1-Data)*log(1-pr)
  LL <- sum(LLi)
  val <- list(LLi=LLi, LL=LL)
  return(val)
}

##�A���S���Y���̐ݒ�
R <- 10000
keep <- 4
iter <- 0

##���O���z�̐ݒ�
mu0 <- 0
sigma0 <- 1
mu0_vec <- rep(0, 3)
cov0 <- diag(0.01, 3)

##�����l�̐ݒ�
#�팱�ҕꐔ�̏����l
r <- as.integer(rank(rowMeans(Data)))
rand <- sort(rnorm(hh, 0, 1), decreasing=TRUE)
oldtheta <- rand[r]

#���ڕꐔ�̏����l
oldalpha <- runif(k, 0.3, 0.8)
r <- as.integer(rank(colMeans(Data)))
rand <- sort(rnorm(k, 0.5, 1), decreasing=TRUE)
oldbeta <- rand[r]
oldgamma <- rep(0.2, k)

##�T���v�����O���ʂ̕ۑ��p�z��
THETA <- matrix(0, nrow=R/keep, ncol=hh)
ALPHA <- matrix(0, nrow=R/keep, ncol=k)
BETA <- matrix(0, nrow=R/keep, ncol=k)
GAMMA <- matrix(0, nrow=R/keep, ncol=k)

####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�팱�ҕꐔtheta��MH�@�ŃT���v�����O
  #�V�����p�����[�^���T���v�����O
  thetad <- oldtheta
  thetan <- thetad + rnorm(hh, 0, 0.25)
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew1 <- rowSums(loglike(Data, thetan, oldbeta, oldalpha, oldgamma, hh, k)$LLi)
  logold1 <- rowSums(loglike(Data, thetad, oldbeta, oldalpha, oldgamma, hh, k)$LLi)
  logpnew1 <- dnorm(thetan, mu0, sigma0, log=TRUE)
  logpold1 <- dnorm(thetad, mu0, sigma0, log=TRUE)
  
  #MH�T���v�����O
  rand <- runif(hh)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew1 + logpnew1 - logold1 - logpold1)   #�̑𗦂��v�Z
  alpha1 <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- (alpha1 >= rand)*1 + (alpha1 < rand)*0
  oldtheta <- flag*thetan + (1-flag)*thetad #alpha��rand�������Ă�����̑�
  
  
  ##���ڕꐔ��MH�@�ŃT���v�����O
  #�V�����p�����[�^���T���v�����O
  alphad <- oldalpha
  betad <- oldbeta
  gammad <- oldgamma
  oldpar <- cbind(alphad, betad, gammad)
  alphan <- alphad + rnorm(k, 0, 0.1)
  betan <- betad + rnorm(k, 0, 0.1)
  gamman0 <- gammad + rnorm(k, 0, 0.01)
  gamman <- ifelse(gamman0 < 0, 0, gamman0)
  newpar <- cbind(alphan, betan, gamman)
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew2 <- colSums(loglike(Data, oldtheta, betan, alphan, gamman, hh, k)$LLi)
  logold2 <- colSums(loglike(Data, oldtheta, betad, alphad, gammad, hh, k)$LLi)
  logpnew2 <- apply(newpar, 1, function(x) lndMvn(x, mu0_vec, cov0))
  logpold2 <- apply(oldpar, 1, function(x) lndMvn(x, mu0_vec, cov0))
  
  #MH�T���v�����O
  rand <- runif(k)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew2 + logpnew2 - logold2 - logpold2)   #�̑𗦂��v�Z
  alpha2 <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- (alpha2 >= rand)*1 + (alpha2 < rand)*0
  par <- flag*newpar + (1-flag)*oldpar   #alpha��rand�������Ă�����̑�
  oldalpha <- par[, 1]; oldbeta <- par[, 2]; oldgamma <- par[, 3]
  
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[mkeep, ] <- oldtheta
    ALPHA[mkeep, ] <- oldalpha
    BETA[mkeep, ] <- oldbeta
    GAMMA[mkeep, ] <- oldgamma
    print(rp)
    print(sum(lognew2))
    print(c(mean(alpha1), mean(alpha2)))
    print(round(rbind(oldtheta, theta0)[, 1:15], 3))
    print(round(cbind(newpar[1:10, ], cbind(alpha0, beta0, c0)[1:10, ]), 3))
  }
}

####�T���v�����O���ʂ̉����Ɨv��####
burnin <- 500
RS <- R/keep

##�T���v�����O���ʂ̉���
matplot(ALPHA[, 1:5], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(ALPHA[, 6:10], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(BETA[, 1:5], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(BETA[, 6:10], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(GAMMA[, 1:5], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(GAMMA[, 6:10], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(THETA[, 1:5], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(THETA[, 6:10], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(THETA[, 11:15], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

##�T���v�����O���ʂ̗v��
round(cbind(colMeans(ALPHA[burnin:RS, ]), alpha0), 3)
round(apply(ALPHA[burnin:RS, ], 1, sd), 3)
round(cbind(colMeans(BETA[burnin:RS, ]), beta0), 3)
round(apply(BETA[burnin:RS, ], 1, sd), 3)
round(cbind(colMeans(GAMMA[burnin:RS, ]), c0), 3)
round(apply(GAMMA[burnin:RS, ], 1, sd), 3)
round(cbind(colMeans(THETA[burnin:RS, ]), theta0), 3)
round(apply(THETA[burnin:RS, ], 1, sd), 3)
