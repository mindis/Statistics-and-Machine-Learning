#####�K�w�x�C�Y�v���r�b�g���f��#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)
source("bdiag_m.R")

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 2000   #���[�U�[��
pt <- 90   #�ϑ����Ԑ�
hhpt <- hh*pt   #�����R�[�h��
w <- 7   #��������
k1 <- 7   #�p�����[�^��
k2 <- 2   #���ݕϐ��̃p�����[�^��

##ID�̐ݒ�
u_id <- rep(1:hh, rep(pt, hh))
t_id <- rep(1:pt, hh)


##�����ϐ��̔���
#�T�����̐���
x <- matrix(diag(w), nrow=pt, ncol=w, byrow=T)
week <- matrix(as.numeric(t(x)), nrow=hhpt, ncol=w, byrow=T)[, -1]

#�ϑ��ϐ��̐���
Z1 <- log(rgamma(pt, 5, 1)) * rbinom(pt, 1, 0.3)   #�~���ʂ̑ΐ�
Z2 <- runif(pt, 0, 1)   #�`���V���i�̕��ϒl����
Z3 <- log(rpois(pt, rgamma(pt, 25, 0.5)))   #�`���V�f�ڏ��i���̑ΐ�
Z3 <- Z3 - min(Z3)
Z4 <- rbinom(pt, 1, 0.5)   #���X�`���V��
Zi1 <- matrix(as.numeric(t(cbind(Z1, Z2, Z3, Z4))), nrow=hhpt, ncol=4, byrow=T)


##���ݕϐ��̔���
#�ƒ���݌ɏ󋵂̏����l
Z5 <- Z6 <- rep(0, hhpt)
Z5[t_id==1] <- round(rgamma(hh, 1, 1), 2)
Z6[t_id==1] <- round(rgamma(hh, 1, 1), 2)
Zi2 <- cbind(Z5, Z6)
C1 <- rgamma(hh, 20, 45)
C2 <- rgamma(hh, 15, 40)
C_INV1 <- Z5[t_id==1]
C_INV2 <- Z6[t_id==1]


##�p�����[�^�̐ݒ�
#��A�W���̃p�����[�^
beta0 <- rnorm(hh, 0.2, 0.4)
betaw <- mvrnorm(hh, rep(0, w-1), diag(0.1, w-1))
colnames(betaw) <- c("betaw1", "betaw2", "betaw3", "betaw4", "betaw5", "betaw6")
beta1 <- rnorm(hh, -0.5, 0.3)
beta2 <- rnorm(hh, 0.6, 0.25)
beta3 <- rnorm(hh, 0.4, 0.25)
beta4 <- rnorm(hh, -0.8, 0.2)
beta5 <- rnorm(hh, -0.6, 0.15)
beta6 <- rnorm(hh, -0.5, 0.15)
beta <- betat <- cbind(beta0, betaw, beta1, beta2, beta3, beta4, beta5, beta6)
colnames(betat) <- NULL

#���ݕϐ��̃p�����[�^
delta1 <- deltat1 <- exp(rnorm(hh, -0.4, 0.25))
delta2 <- deltat2 <- exp(rnorm(hh, -0.7, 0.25))
lambda1 <- lambdat1 <- rnorm(hh, 0.3, 0.2)
lambda2 <- lambdat2 <- rnorm(hh, 0.6, 0.2)
theta <- thetat <- cbind(delta1, delta2, lambda1, lambda2)

##�����ϐ��̔���
y <- rep(0, hhpt); mu_vec <- rep(0, hhpt)
value <- matrix(0, nrow=hhpt, ncol=2)
INV1 <- rep(0, hhpt); INV2 <- rep(0, hhpt)

for(j in 1:pt){

  #���݌��p�̕��ς𐶐�
  index <- which(t_id==j)
  Zi <- cbind(1, week[index, ], Zi1[index, ], Zi2[index, ])
  mu <- beta0 + rowSums(Zi * beta)
  mu_vec[index] <- mu
  
  #�����ϐ��𐶐�
  U <- rnorm(hh, mu, 1)   #���p�֐�
  y[index] <- ifelse(U > 0, 1, 0)   #���X�L���𐶐�
  
  if(j < pt){
    #�w������΁A�ƒ���݌ɗʂɒǉ�
    index1 <- which(y==0 & t_id==j) + 1; index2 <- which(y==1 & t_id==j) + 1
    value[index2-1, 1] <- rgamma(length(index2), 10, 30); value[index2-1, 2] <- rgamma(length(index2), 10, 35)   #�w�����z
    INV1[index1] <- Zi2[index1-1, 1]
    INV2[index1] <- Zi2[index1-1, 2]
    INV1[index2] <- Zi2[index2-1, 1] + value[index2-1, 1]
    INV2[index2] <- Zi2[index2-1, 2] + value[index2-1, 2]
    
    #�ƒ���݌ɗʂ��X�V
    index <- which(t_id==(j+1))
    Zi2[index, 1] <- INV1[index] - INV1[index]*(delta1*C1 / (delta1*C1 + INV1[index]^lambda1))
    Zi2[index, 2] <- INV2[index] - INV2[index]*(delta2*C1 / (delta2*C1 + INV2[index]^lambda2))
  }
}
Zit2 <- Zi2
mean(y)   #���X�m��
z <- round(cbind(u_id, y, value, Zi2), 3)   #�ƒ���݌ɗʂƍw���ʂ��m�F



####�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y�v���r�b�g���f���𐄒�####
##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##�A���S���Y���̐ݒ�
R <- 10000
keep <- 4
burnin <- 2000/keep
disp <- 8
k1 <- ncol(cbind(1, week, Zi1, Zi2))
k2 <- 4

##���O���z�̐ݒ�
sigma <- 1
beta0 <- 0
tau0 <- 1/100

#�ϗʌ��ʂ̎��O���z

Bbar1 <- rep(0, k1)
Bbar2 <- rep(0, k2)
A1 <- 0.01 * diag(1, k1)
nu1 <- k1*2
nu2 <- k2*2
V1 <- nu1 * diag(k1)
V2 <- nu2 * diag(k2)


#��A�W���̏����l
beta_fixed <- c(as.numeric(glm(y ~ cbind(week, Z1, Z2, Z3, Z4), family=binomial(link="probit"))$coef), -0.1, -0.1)
beta <- mvrnorm(hh, beta_fixed, diag(0.1, k1))
delta1 <- exp(rnorm(hh, -0.5, 0.1))
delta2 <- exp(rnorm(hh, -0.5, 0.1))
lambda1 <- rnorm(hh, 0.5, 0.1)
lambda2 <- rnorm(hh, 0.5, 0.1)
theta <- cbind(delta1, delta2, lambda1, lambda2)

#�K�w���f���̏����l
alpha1 <- rep(0, k1)
alpha2 <- colMeans(cbind(deltat1, deltat2, lambdat1, lambdat2))
alpha2 <- rep(1, k2)
Cov1 <- diag(0.1, k1)
Cov2 <- diag(0.1, k2)
Cov_inv1 <- solve(Cov1)
Cov_inv2 <- solve(Cov2)

#���ݕϐ��̏����l
INV1 <- rep(0, hhpt); INV2 <- rep(0, hhpt)
Zi2 <- matrix(0, nrow=hhpt, ncol=2)
Zi2[t_id==1, ] <- Zit2[t_id==1, ]

for(j in 1:pt){
  if(j < pt){
    #�w������΁A�ƒ���݌ɗʂɒǉ�
    index1 <- which(y==0 & t_id==j) + 1; index2 <- which(y==1 & t_id==j) + 1
    INV1[index1] <- Zi2[index1-1, 1]
    INV2[index1] <- Zi2[index1-1, 2]
    INV1[index2] <- Zi2[index2-1, 1] + value[index2-1, 1]
    INV2[index2] <- Zi2[index2-1, 2] + value[index2-1, 2]
    
    #�ƒ���݌ɗʂ��X�V
    index <- which(t_id==j+1)
    Zi2[index, 1] <- INV1[index] - INV1[index]*(delta1*C1 / (delta1*C1 + INV1[index]^lambda1))
    Zi2[index, 2] <- INV2[index] - INV2[index]*(delta2*C2 / (delta2*C2 + INV2[index]^lambda2))
  }
}

##�p�����[�^�̊i�[�p�z��
BETA <- array(0, dim=c(hh, k1, R/keep))
THETA <- array(0, dim=c(hh, k2, R/keep))
ALPHA1 <- matrix(0, nrow=R/keep, ncol=k1)
ALPHA2 <- matrix(0, nrow=R/keep, ncol=k2)
COV1 <- array(0, dim=c(k1, k1, R/keep))
COV2 <- array(0, dim=c(k2, k2, R/keep))
Zi2_mu <- matrix(0, nrow=hhpt, ncol=2)

##�C���f�b�N�X��ݒ�
index_t <- list(); index_t1 <- list(); index_t2 <- list()
index_u <- list()
for(j in 1:pt){
  index_t[[j]] <- which(t_id==j)
  index_t1[[j]] <- which(t_id==j & y==0)
  index_t2[[j]] <- which(t_id==j & y==1)
}
for(i in 1:hh){
  index_u[[i]] <- which(u_id==i)
}

##�ؒf�̈���`
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�ؒf���K���z������݌��p�𐶐�
  Zi <- cbind(1, week, Zi1, Zi2)
  mu <- rowSums(Zi * beta[u_id, ])   #���p�֐��̕���
  U <- rtnorm(mu, sigma, a, b)   #���݌��p�𐶐�

  ##�l�ʂ�beta���T���v�����O
  for(i in 1:hh){
    #�f�[�^�𒊏o
    index <- index_u[[i]]
    x <- Zi[index, ]
    
    #��A�W���̎��㕪�z�̃p�����[�^
    Xy <- t(x) %*% U[index]
    XXV <- solve(t(x) %*% x + Cov_inv1)
    beta_mu <- XXV %*% (Xy + Cov_inv1 %*% alpha1)
  
    #���ϗʐ��K���z����beta���T���v�����O
    beta[i, ] <- mvrnorm(1, beta_mu, sigma*XXV)
  }
  
  ##MH�@�Ő��ݕϐ��̃p�����[�^���T���v�����O
  #�V�����p�����[�^���T���v�����O
  thetad <- theta
  thetan <- thetad + mvrnorm(hh, rep(0, k2), diag(0.01, k2))
  thetan[, 1:2] <- abs(thetan[, 1:2])
  Zi_n2 <- Zi_d2 <- Zi2
  
  #�ƒ���݌ɗʂ��X�V
  Zi_n2[-index_t[[1]], ] <- 0
  INV1 <- rep(0, hhpt); INV2 <- rep(0, hhpt)

  for(j in 1:pt){
    if(j < pt){
      #�w������΁A�ƒ���݌ɗʂɒǉ�
      index1 <- index_t1[[j]]; index2 <- index_t2[[j]]
      index <- index_t[[j+1]] 
      INV1[index1+1] <- Zi_n2[index1, 1]
      INV2[index1+1] <- Zi_n2[index1, 2]
      INV1[index2+1] <- Zi_n2[index2, 1] + value[index2, 1]
      INV2[index2+1] <- Zi_n2[index2, 2] + value[index2, 2]
      
      #�ƒ���݌ɗʂ��X�V
      Zi_n2[index, 1] <- INV1[index] - INV1[index]*(thetan[, 1]*C1 / (thetan[, 1]*C1 + INV1[index]^thetan[, 3]))
      Zi_n2[index, 2] <- INV2[index] - INV2[index]*(thetan[, 2]*C2 / (thetan[, 2]*C2 + INV2[index]^thetan[, 4]))
    }
  }
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew <- logold <- rep(0, hh)
  Zi_n <- cbind(1, week, Zi1, Zi_n2); Zi_d <- Zi
  mu1 <- rowSums(Zi_n * beta[u_id, ])
  mu2 <- rowSums(Zi_d * beta[u_id, ])

  for(i in 1:hh){
    er1 <- (thetan[i, ] - alpha2)
    er2 <- (thetad[i, ] - alpha2)
    lognew[i] <- -1/2 * (sum((U[index_u[[i]]] - mu1[index_u[[i]]])^2) + er1 %*% Cov_inv2 %*% er1)
    logold[i] <- -1/2 * (sum((U[index_u[[i]]] - mu2[index_u[[i]]])^2) + er2 %*% Cov_inv2 %*% er2)
  }

  #MH�T���v�����O
  rand <- runif(hh)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew - logold)   #�̑𗦂��v�Z
  omega <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((omega >= rand)*1 + (omega < rand)*0), nrow=hh, ncol=k2)
  theta <- flag*thetan + (1-flag)*thetad   #alpha��rand�������Ă�����̑�
  Zi2 <- flag[u_id, 1:(k2/2)]*Zi_n2 + (1-flag[u_id, 1:(k2/2)])*Zi_d2
  
  ##��A�W���̊K�w���f���̃p�����[�^���T���v�����O
  #�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
  V_par <- V1 + t(beta) %*% beta
  Sn <- nu1 + hh
  Cov1 <- bayesm::rwishart(Sn, solve(V_par))$IW
  Cov_inv1 <- solve(Cov1)
  
  #���ϗʐ��K���z���畽�σx�N�g�����T���v�����O
  beta_mu <- hh/(hh + tau0) * colMeans(beta)
  alpha1 <- mvrnorm(1, beta_mu, Cov1/(hh + tau0))
  
  
  ##���ݕϐ��̊K�w���f���̃p�����[�^���T���v�����O
  #�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
  V_par <- V2 + t(theta) %*% theta
  Sn <- nu2 + hh
  Cov2 <- bayesm::rwishart(Sn, solve(V_par))$IW
  Cov_inv2 <- solve(Cov2)

  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    BETA[, , mkeep] <- beta
    THETA[, , mkeep] <- theta
    ALPHA1[mkeep, ] <- alpha1
    ALPHA2[mkeep, ] <- alpha2
    COV1[, , mkeep] <- Cov1
    COV2[, , mkeep] <- Cov2
    
    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(rp%%keep==0 & rp >= burnin){
      Zi2_mu <- Zi2_mu + Zi2
    }
    
    ##�T���v�����O���ʂ��m
    if(rp%%disp==0){
      print(rp)
      print(sum(dbinom(y, 1, pnorm(mu, 0, 1), log=TRUE)))
      print(round(diag(Cov1), 3))
      print(round(rbind(alpha1, alphat1=colMeans(betat)), 3))
    }
  }
}


####�T���v�����O���ʂ̗v��Ɖ���####
RS <- R/keep
burnin <- 2000/keep

##�T���v�����O���ʂ̉���
matplot(t(BETA[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(BETA[10, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(ALPHA1, type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

##���㕽��
round(apply(BETA[, , burnin:RS], c(1, 2), mean), 3)   #�l�ʂ�beta�̎��㕽��
round(colMeans(ALPHA1[burnin:RS, ]), 3)   #�K�w���f���̉�A�W���̎��㕽��
round(apply(COV1[, , burnin:RS], c(1, 2), mean), 3)   #�K�w���f���̕��U�����U�s��̎��㕽��
