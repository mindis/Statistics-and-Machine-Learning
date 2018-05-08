#####mixed effect poisson block model#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(flexmix)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

set.seed(506832)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
d <- 150   #�A�C�e����
k <- 7   #���ݕϐ���
N <- d*(d-1)/2   #���T���v����
vec <- rep(1, k)

##ID��ݒ�
id1 <- id2 <- c()
for(i in 1:(d-1)){
  id1 <- c(id1, rep(i, length((i+1):d)))
  id2 <- c(id2, (i+1):d)
}

##���ݕϐ��̐���
#�f�B���N�����z����p�����[�^�𐶐�
alpha0 <- rep(0.1, k)
theta <- thetat <- extraDistr::rdirichlet(d, alpha0)
Z1 <- rmnom(N, 1, theta[id1, ])
Z2 <- rmnom(N, 1, theta[id2, ])
z1_vec <- as.numeric(Z1 %*% 1:k); z2_vec <- as.numeric(Z2 %*% 1:k)


##�����ϐ��̐���
#�p�����[�^�𐶐�
cov <- covt <- 0.75
beta <- betat <- 0.8
alphat <- alpha <- rnorm(d, beta, cov)   #�ϗʌ��ʂ̃p�����[�^
phi0 <- matrix(rnorm(k*k, 0, 0.85), nrow=k, ncol=k)   #���ݕϐ��̃p�����[�^
phi0[upper.tri(phi0)] <- 0
phi <- phi0 + t(phi0)
diag(phi) <- diag(phi0)
phit <- phi

#�|�A�\�����z�̕��ύ\��
mu <- alpha[id1] + alpha[id2] + (phi[z1_vec, ] * Z2) %*% vec
lambda <- exp(mu)


#�|�A�\�����z���牞���ϐ��𐶐�
y <- rpois(N, lambda)
sum(y); mean(y)
hist(y, xlab="�p�x", main="�A�C�e���Ԃ̏o���p�x", col="grey", breaks=25)


####�}���R�t�A�������e�J�����@��mixed effect poisson block model�𐄒�####
##�ϗʌ��ʃ|�A�\����A���f���̑ΐ��ޓx
loglike1 <- function(alpha, theta, y, y_factorial, z1, Z2, vec, id1, id2){
  
  #�ޓx���`����
  lambda <- exp(alpha[id1] + alpha[id2] + (phi[z1, ] * Z2) %*% vec)   #���ύ\��
  LLi <- as.numeric(y*log(lambda)-lambda - y_factorial)   #�ΐ��ޓx
  LL <- sum(LLi)   #�ΐ��ޓx�̘a
  
  #���ʂ�Ԃ�
  LL_value <- list(LLi=LLi, LL=LL)
  return(LL_value)
}

loglike2 <- function(alpha, phi, y, y_factorial, index, id1, id2){
  #�ޓx���`����
  lambda <- exp(alpha[id1[index]] + alpha[id2[index]] + phi)
  LL <- sum(y[index]*log(lambda)-lambda - y_factorial[index])
  return(LL)
}


##�A���S���Y���̐ݒ�
R <- 10000
keep <- 4  
iter <- 0
burnin <- 1000/keep
disp <- 100

##���O���z�̐ݒ�
#��A�W���̎��O���z
alpha01 <- 0
beta01 <- 0
tau01 <- 0.01

#�ϗʌ��ʂ̎��O���z
alpha02 <- 0
s02 <- 0.01
v02 <- 0.01

#�f�B���N�����z�̎��O���z
alpha03 <- 0.1


##�p�����[�^�̐^�l
theta <- thetat
phi <- phit
alpha <- alphat
beta <- betat
cov <- covt
Zi1 <- Z1; Zi2 <- Z2
z_vec1 <- as.numeric(Zi1 %*% 1:k)
z_vec2 <- as.numeric(Zi2 %*% 1:k)


##�����l�̐ݒ�
#�ϗʌ��ʂ̏����l
cov <- 0.5   #�K�w���f���̕W���΍��̏����l
beta <- 0.8   #�K�w���f���̕��ς̏����l
mu <- rep(0, d)
for(i in 1:d){
  mu[i] <- mean(y[which(id1==i | id2==i)])
}
alpha <- rnorm(d, 0, 0.5)
rank_mu <- ceiling(rank(mu))
alpha <- sort(x, decreasing=TRUE)[rank_mu]   #�ϗʌ��ʂ̏����l

#���ݕϐ��̃p�����[�^
phi0 <- matrix(rnorm(k*k, 0, 0.5), nrow=k, ncol=k)   #���ݕϐ��̃p�����[�^
phi0[upper.tri(phi0)] <- 0
phi <- phi0 + t(phi0)

#���ݕϐ��̏����l
theta <- extraDistr::rdirichlet(d, rep(0.2, k))
Zi1 <- rmnom(N, 1, theta[id1, ])
Zi2 <- rmnom(N, 1, theta[id2, ])
z1_vec <- as.numeric(Zi1 %*% 1:k); z2_vec <- as.numeric(Zi2 %*% 1:k)


##�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, k, R/keep))
ALPHA <- matrix(0, nrow=R/keep, ncol=d)
BETA <- rep(0, R/keep)
COV <- rep(0, R/keep)
SEG2 <- SEG1 <- matrix(0, nrow=N, ncol=k)

##�萔���v�Z
y_factorial <- lfactorial(y)
Y_factorial <- matrix(y_factorial, nrow=N, ncol=k*k)
upper_tri <- matrix(as.logical(upper.tri(phi) + diag(1, k)), nrow=k, ncol=k)
vec <- rep(1, k)

##�C���f�b�N�X��ݒ�
#���̓f�[�^�̃C���f�b�N�X
item_list <- list()
seg_list1 <- seg_list2 <- list()
index_list1 <- index_list2 <- list()
for(i in 1:d){
  item_list[[i]] <- which(id1==i | id2==i)
  seg_list1[[i]] <- as.numeric(id1[item_list[[i]]]!=i)
  seg_list2[[i]] <- as.numeric(id2[item_list[[i]]]!=i)
  index_list1[[i]] <- which(item_list[[i]] * seg_list1[[i]] != 0)
  index_list2[[i]] <- which(item_list[[i]] * seg_list2[[i]] != 0)
}
item_vec <- rep(1, d-1)

#�p�����[�^�̃C���f�b�N�X
par <- c()
for(i in 1:k){
  for(j in i:k){
    par <- c(par, i, j)
  }
}
index_par <- matrix(par, nrow=k*(k-1)/2+k, ncol=2, byrow=T)


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
    
  ##���g���|���X�w�C�X�e�B���O�@�ŕϗʌ��ʂ��T���v�����O
  #�V�����p�����[�^���T���v�����O
  alphad <- alpha
  alphan <- alphad + rnorm(d, 0, 0.1)
  
  #���O���z�̌덷
  er_new <- alphan - beta 
  er_old <- alphad - beta
  
  #�ΐ��ޓx�Ƒΐ����O���z��ݒ�
  lognew0 <- loglike1(alphan, phi, y, y_factorial, z1_vec, Zi2, vec, id1, id2)$LLi
  logold0 <- loglike1(alphad, phi, y, y_factorial, z1_vec, Zi2, vec, id1, id2)$LLi
  logpnew1 <- -0.5 * (er_new^2 / cov)
  logpold1 <- -0.5 * (er_old^2 / cov)
  
  #�A�C�e�����Ƃɑΐ��ޓx�̘a�����
  lognew1 <- logold1 <- rep(0, d)
  for(i in 1:d){
    lognew1[i] <- sum(lognew0[item_list[[i]]])
    logold1[i] <- sum(logold0[item_list[[i]]])
  }
  
  #MH�T���v�����O
  rand <- runif(d)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew1 + logpnew1 - logold1 - logpold1)   #�̑𗦂��v�Z
  LLind_diff <- ifelse(LLind_diff==Inf, 1, ifelse(LLind_diff==-Inf, 0, LLind_diff))
  alpha2 <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- (alpha2 >= rand)*1 + (alpha2 < rand)*0
  alpha <- flag*alphan + (1-flag)*alphad   #alpha��rand�������Ă�����̑�
  
  
  ##�K�w���f���̃p�����[�^���T���v�����O
  #���K���z���畽�σp�����[�^���T���v�����O
  #beta_mu <- d/(d + tau01) * mean(alpha)
  #beta <- rnorm(1, beta_mu, cov/(d + tau01))
  
  #�t�K���}���z����W���΍����T���v�����O
  s <- s02 + sum((alpha - mean(alpha))^2)
  v <- v02 + d
  cov <- sqrt(1/rgamma(1, v/2, s/2))
  
  
  ##�M�u�X�T���v�����O�Ő��ݕϐ����T���v�����O
  #�������z������ݕϐ����T���v�����O
  for(i in 1:d){
    
    #�C���f�b�N�X�𒊏o
    index <- item_list[[i]]
    z1_allocation <- seg_list1[[i]]
    z2_allocation <- seg_list2[[i]]
    
    #�|�A�\�����z�̕��ύ\��
    Zi_pairs <- Zi1[index, ] * z1_allocation + Zi2[index, ] * z2_allocation   #�΂ƂȂ���ݕϐ�
    Zi_mu <- phi[as.numeric(Zi_pairs %*% 1:k), ]
    lambda <- exp(matrix(alpha[id1[index]] + alpha[id2[index]], nrow=length(index), ncol=k) + Zi_mu)
    
    #�ΐ��ޓx�Ɛ��ݕϐ��̊����m��
    LLi <- y[index]*log(lambda)-lambda - y_factorial[index]   #�ΐ��ޓx
    z_par <- exp(LLi - rowMaxs(LLi)) * matrix(theta[i, ], nrow=length(index), ncol=k, byrow=T)
    z_rate <- z_par / rowSums(z_par)
    
    #�������z�����ݕϐ����T���v�����O
    index1 <- index_list1[[i]]; index2 <- index_list2[[i]]
    if(length(index1) > 0){
      Zi1[index1, ] <- rmnom(length(index1), 1, z_rate[index1, ])
    } 
    if(length(index2) > 0){
      Zi2[index2, ] <- rmnom(length(index2), 1, z_rate[index2, ])
    }
  }
  #���ݕϐ��s���ϊ�
  z1_vec <- as.numeric(Zi1 %*% 1:k)
  z2_vec <- as.numeric(Zi2 %*% 1:k)
  Zi1_T <- t(Zi1); Zi2_T <- t(Zi2)
  
  
  ##MH�@�Ő��ݕϐ��̃p�����[�^���T���v�����O
  for(i in 1:nrow(index_par)){
    
    #�C���f�b�N�X��ݒ�
    index_phi <- index_par[i, ]
    index_z <- which(Zi1[, index_phi[1]]*Zi2[, index_phi[2]]==1)
    
    #�V�����p�����[�^���T���v�����O
    phid <- phi[index_phi[1], index_phi[2]]
    phin <- phid + rnorm(1, 0, 0.1)
    
    #�ΐ��ޓx�Ƒΐ����O���z��ݒ�
    lognew2 <- loglike2(alpha, phin, y, y_factorial, index_z, id1, id2)
    logold2 <- loglike2(alpha, phid, y, y_factorial, index_z, id1, id2)
    logpnew2 <- -0.5 * (phin^2 / (1/tau01))
    logpold2 <- -0.5 * (phid^2 / (1/tau01))
    
    #MH�T���v�����O
    rand <- runif(1)   #��l���z���痐���𔭐�
    LLind_diff <- exp(lognew2 + logpnew2 - logold2 - logpold2)   #�̑𗦂��v�Z
    if(LLind_diff==Inf) LLind_diff <- 1; if(LLind_diff==-Inf) LLind_diff <- 0   #Inf�̏ꍇ�̏���
    alpha2 <- (LLind_diff >= 1)*1 + (LLind_diff < 1)*LLind_diff
    
    #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
    flag <- (alpha2 >= rand)*1 + (alpha2 < rand)*0
    phi[index_phi[1], index_phi[2]] <- flag*phin + (1-flag)*phid   #alpha��rand�������Ă�����̑�
  }
  phi[lower.tri(phi)] <- phi[upper.tri(phi)]   #�p�����[�^��Ώۍs��ɕύX
  phi <- phi - mean(phi)   #�p�����[�^�𒆉���
  
  
  ##�f�B���N�����z����g�s�b�N���z���T���v�����O
  wsum0 <- matrix(0, nrow=d, ncol=k)
  for(i in 1:d){
    index1 <- item_list[[i]] * seg_list1[[i]]
    index2 <- item_list[[i]] * seg_list2[[i]]
    wsum0[i, ] <- cbind(Zi1_T[, index2], Zi2_T[, index1]) %*% item_vec
  }
  wsum <- wsum0 + alpha03   #�f�B���N�����z�̃p�����[�^
  theta <- extraDistr::rdirichlet(d, wsum)   #�f�B���N�����z����g�s�b�N�𐶐�
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    ALPHA[mkeep, ] <- alpha
    BETA[mkeep] <- beta
    COV[mkeep] <- cov
  }
  
  #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
  if(rp%%keep==0 & rp >= burnin){
    SEG1 <- SEG1 + Zi1
    SEG2 <- SEG2 + Zi2
  }
  
  if(rp%%disp==0){
    #�T���v�����O���ʂ̕\��
    print(rp)
    print(sum(lognew1))
    print(round(cbind(phi, phit), 3))
  }
}
matplot(ALPHA[, 1:10], type="l")
matplot(t(PHI[1, , ]), type="l")
plot(1:(R/keep), COV, type="l")