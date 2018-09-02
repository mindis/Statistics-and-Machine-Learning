#####HMC�K�w�x�C�Y���W�X�e�B�b�N��A���f��#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(FAdist)
library(bayesm)
library(extraDistr)
library(condMVNorm)
library(actuar)
library(gtools)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(34027)
##�f�[�^�̐ݒ�
hh <- 5000   #����Ґ�
pt <- rtpois(hh, rgamma(hh, 15.0, 0.25), a=1, b=Inf)   #�w���ڐG��
hhpt <- sum(pt)   #�S�T���v����


#ID��ݒ�
u_id <- rep(1:hh, pt)
t_id <- as.numeric(unlist(tapply(1:hhpt, u_id, rank)))
ID <- data.frame(no=1:hhpt, id=u_id, t=t_id)   #�f�[�^�̌���

##�f�[�^�̔���
##�̓����f���̐����ϐ��̔���
#�A���ϐ�
cont <- 4 
x1 <- matrix(runif(hhpt*cont, 0, 1), nrow=hhpt, ncol=cont) 

#��l�ϐ�
bin <- 3
x2 <- matrix(0, nrow=hhpt, ncol=bin)
for(j in 1:bin){
  prob <- runif(1, 0.2, 0.8)
  x2[, j] <- rbinom(hhpt, 1, prob)  
}

#���l�ϐ�
multi <- 4
x3 <- rmnom(hhpt, 1, extraDistr::rdirichlet(1, rep(2.5, multi))); x3 <- x3[, -which.min(colSums(x3))]

#�f�[�^�̌���
x <- cbind(1, x1, x2, x3)
k1 <- ncol(x)

##�̊ԃ��f���̐����ϐ��̔���
#�A���ϐ�
cont <- 3
u1 <- matrix(runif(hh*cont, 0, 1), nrow=hh, ncol=cont) 

#��l�ϐ�
bin <- 4
u2 <- matrix(0, nrow=hh, ncol=bin)
for(j in 1:bin){
  prob <- runif(1, 0.2, 0.8)
  u2[, j] <- rbinom(hh, 1, prob)  
}

#���l�ϐ�
multi <- 5
u3 <- rmnom(hh, 1, extraDistr::rdirichlet(1, rep(2.5, multi))); u3 <- u3[, -which.min(colSums(u3))]

#�f�[�^�̌���
u <- cbind(1, u1, u2, u3)
k2 <- ncol(u)


##��A�W���̐ݒ�
##�̊ԉ�A�W���̐ݒ�
#�Ó��Ȕ����ϐ����o����܂ŉ�A�W����ݒ肵����
repeat {

  #�K�w���f���̃p�����[�^�𐶐�
  Cov <- Covt <- diag(runif(k1, 0.01, 0.2))
  theta <- thetat <-  matrix(rnorm(k1*k2, runif(k1*k2, -0.4, 0.3), 0.75), nrow=k2, ncol=k1)
  
  #�̓���A�W���𐶐�
  beta <- betat <- u %*% theta + mvrnorm(hh, rep(0, k1), Cov)
  
  #�m���̔���
  logit <- as.numeric((x * beta[u_id, ]) %*% rep(1, k1))
  prob <- exp(logit) / (1 + exp(logit))
  
  #�����ϐ��𐶐�
  y <- rbinom(hhpt, 1, prob)
  if(mean(y) > 0.20 & mean(y) < 0.5) break 
}

#�����ϐ��̗v��
summary(prob)   #�����ϐ��̗v�񓝌v��
mean(y)   #�����m��
hist(prob, breaks=25, col="grey", xlab="�����m��", main="�����m���̕��z")


##�e�X�g�f�[�^���쐬
##�̓����f���̐����ϐ��̔���
#�A���ϐ�
cont <- 4 
x1 <- matrix(runif(hhpt*cont, 0, 1), nrow=hhpt, ncol=cont) 

#��l�ϐ�
bin <- 3
x2 <- matrix(0, nrow=hhpt, ncol=bin)
for(j in 1:bin){
  prob <- runif(1, 0.2, 0.8)
  x2[, j] <- rbinom(hhpt, 1, prob)  
}

#���l�ϐ�
multi <- 4
x3 <- rmnom(hhpt, 1, extraDistr::rdirichlet(1, rep(2.5, multi))); x3 <- x3[, -which.min(colSums(x3))]

#�f�[�^�̌���
x_test <- cbind(1, x1, x2, x3)
k1 <- ncol(x)

##�̊ԃ��f���̐����ϐ��̔���
#�A���ϐ�
cont <- 3
u1 <- matrix(runif(hh*cont, 0, 1), nrow=hh, ncol=cont) 

#��l�ϐ�
bin <- 4
u2 <- matrix(0, nrow=hh, ncol=bin)
for(j in 1:bin){
  prob <- runif(1, 0.2, 0.8)
  u2[, j] <- rbinom(hh, 1, prob)  
}

#���l�ϐ�
multi <- 5
u3 <- rmnom(hh, 1, extraDistr::rdirichlet(1, rep(2.5, multi))); u3 <- u3[, -which.min(colSums(u3))]

#�f�[�^�̌���
u_test <- cbind(1, u1, u2, u3)


##�����ϐ��𐶐�
#�̓���A�W���𐶐�
beta_test <- u_test %*% theta + mvrnorm(hh, rep(0, k1), Cov)

#�m���̔���
logit <- as.numeric((x_test * beta_test[u_id, ]) %*% rep(1, k1))
prob <- exp(logit) / (1 + exp(logit))
  
#�����ϐ��𐶐�
y_test <- rbinom(hhpt, 1, prob)


####�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y���W�X�e�B�b�N��A���f���𐄒�####
##�ΐ����㕪�z���v�Z����֐�
loglike <- function(beta, y, x, u_mu, inv_Cov, u_id, u_index, hh, k1){
  
  #���W�b�g���f���̑ΐ��ޓx
  mu <- exp(as.numeric((x * beta[u_id, ]) %*% rep(1, k1))) 
  prob <- mu / (1 + mu)   #�m���̌v�Z
  LLi_logit <- y*log(prob) + (1-y)*log(1-prob)   #���W�b�g���f���̑ΐ��ޓx
  
  #���ϗʐ��K���z�̑ΐ��ޓx
  er <- beta - u_mu  #�덷
  LLi_mvn <- -1/2 * as.numeric((er %*% inv_Cov * er) %*% rep(1, k1))
  
  #���[�U�[���Ƃ̑ΐ����㕪�z
  LLi <- rep(0, hh)
  for(i in 1:hh){
    LLi[i] <- sum(LLi_logit[u_index[[i]]]) + LLi_mvn[i] 
  }
  return(LLi)
}

##���W�X�e�B�b�N-���K���z�̑ΐ����㕪�z�̔����֐�
dloglike <- function(beta, y, x, u_mu, inv_Cov, hh, u_id, u_index, k1){
  #�����m���̐ݒ�
  mu <- exp(as.numeric((x * beta[u_id, ]) %*% rep(1, k1)))   #���W�b�g�̎w���֐�
  prob <- mu / (1 + mu)   #�m���̌v�Z
  
  #�����֐��̐ݒ�
  er <- beta - u_mu
  dlogit <- y*x - x*prob   #���W�X�e�B�b�N��A�̑ΐ��ޓx�̔����֐�
  dmvn <- -t(inv_Cov %*% t(er))   #���ϗʐ��K���z�̑ΐ����O���z�̔����֐�

  #�ΐ����㕪�z�̔����֐��̘a
  LLd <- matrix(0, nrow=hh, ncol=k1)
  for(i in 1:hh){
    LLd[i, ] <- -(colSums(dlogit[u_index[[i]], ]) + dmvn[i, ])
  }
  return(LLd)
}

##���[�v�t���b�O�@�������֐�
leapfrog <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, y, x, u_mu, inv_Cov, hh, u_id, u_index, k1) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, y, x, u_mu, inv_Cov, hh, u_id, u_index, k1) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##�A���S���Y���̐ݒ�
R <- 10000
keep <- 2
disp <- 10
burnin <- 2000/keep
iter <- 0
e <- 0.1
L <- 3

##�C���f�b�N�X�ƃf�[�^��ݒ�
#�C���f�b�N�X�̐ݒ�
u_index <- list()
for(i in 1:hh){
  u_index[[i]] <- which(u_id==i)
}


##���O���z�̐ݒ�
Deltabar <- matrix(0, nrow=k2, ncol=k1)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta <- 0.01 * diag(k2)   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu <- k1   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(k1) #�t�E�B�V���[�g���z�̃p�����[�^

##�^�l�̐ݒ�
beta <- betat
theta <- thetat
Cov <- Covt

##�����l�̐ݒ�
beta <- mvrnorm(hh, rep(0, k1), 0.01 * diag(k1))
theta <- matrix(rnorm(k1*k2, 0, 0.1), nrow=k2, ncol=k1); u_mu <- u %*% theta
Cov <- 0.1 * diag(k1); inv_Cov <- solve(Cov)


##�T���v�����O���ʂ̕ۑ��p�z��
BETA <- array(0, dim=c(hh, k1, R/keep))
THETA <- array(0, dim=c(k2, k1,  R/keep))
COV <- matrix(0, nrow=R/keep, ncol=k1)

##�ΐ��ޓx�̊�l
#�^�l�ł̑ΐ��ޓx
mu <- exp(as.numeric((x * betat[u_id, ]) %*% rep(1, k1)))
prob <- mu / (1 + mu)   #�����m��
LLbest <- sum(y*log(prob) + (1-y)*log(1-prob))

#1�p�����[�^�ł̑ΐ��ޓx
LLst <- sum(y*log(mean(y)) + (1-y)*log(1-mean(y)))

##�e�X�g�f�[�^�̑ΐ��ޓx
#�^�l�ł̑ΐ��ޓx
mu <- exp(as.numeric((x_test * (u_test %*% thetat)[u_id, ]) %*% rep(1, k1)))
prob <- mu / (1 + mu)   #�����m��
LLbest_new <- sum(y_test*log(prob) + (1-y_test)*log(1-prob))

#1�p�����[�^�ł̑ΐ��ޓx
LLst_new <- sum(y_test*log(mean(y)) + (1-y_test)*log(1-mean(y)))


####�n�~���g�j�A�������e�J�����@�ɂ��p�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##HMC�ɂ�胆�[�U�[�ʂ̃p�����[�^���T���v�����O
  #HMC�̐V�����p�����[�^�𐶐�
  rold <- mvrnorm(hh, rep(0, k1), diag(k1))   #�W�����ϗʐ��K���z����p�����[�^�𐶐�
  betad <- beta
  
  #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
  res <- leapfrog(rold, betad, dloglike, e, L)   
  rnew <- res$r
  betan <- res$z
  
  #�ړ��O�ƈړ���̃n�~���g�j�A��
  Hnew <- -loglike(betan, y, x, u_mu, inv_Cov, u_id, u_index, hh, k1) + as.numeric(rnew^2 %*% rep(1, k1))/2
  Hold <- -loglike(betad, y, x, u_mu, inv_Cov, u_id, u_index, hh, k1) + as.numeric(rold^2 %*% rep(1, k1))/2
  
  #HMC�@�ɂ��p�����[�^�̍̑�������
  rand <- runif(hh) #��l���z���痐���𔭐�
  alpha <- rowMins(cbind(1, exp(Hold - Hnew)))   #�̑𗦂�����
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- as.numeric(alpha > rand)
  beta <- flag*betan + (1-flag)*betad
  

  ##�K�w���f���̃p�����[�^���T���v�����O
  #���ϗʉ�A���f������K�w���f���̉�A�p�����[�^�ƕ��U���T���v�����O
  out <- rmultireg(beta, u, Deltabar, ADelta, nu, V)
  theta <- out$B
  Cov <- diag(diag(out$Sigma))
  u_mu <- u %*% theta
  inv_Cov <- solve(Cov)
  
  ##�T���v�����O���ʂ̕ۑ��ƕ\��
  #�T���v�����O���ʂ̕ۑ�
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[, , mkeep] <- beta
    THETA[, , mkeep] <- theta
    COV[mkeep, ] <- diag(Cov)
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���Z�o
    mu <- exp(as.numeric((x * beta[u_id, ]) %*% rep(1, k1)))
    prob <- mu / (1 + mu)   #�����m��
    LL <- sum(y*log(prob) + (1-y)*log(1-prob))
    
    #�e�X�g�f�[�^�̑ΐ��ޓx
    mu <- exp(as.numeric((x_test * (u_test %*% theta)[u_id, ]) %*% rep(1, k1)))
    prob <- mu / (1 + mu)   #�����m��
    LL_new <- sum(y_test*log(prob) + (1-y_test)*log(1-prob))
    
    #�T���v�����O���ʂ�\��
    print(rp)
    print(mean(alpha))
    print(c(LL, LLbest, LLst))
    print(c(LL_new, LLbest_new, LLst_new))
    print(round(rbind(diag(Cov), diag(Covt)), 3))
  }
}

####�T���v�����O���ʂ̊m�F�ƓK���x�̊m�F####
#�T���v�����O���ꂽ�p�����[�^���v���b�g
burnin <- 2000/keep
RS <- R/keep

#�T���v�����O���ꂽ�p�����[�^���v���b�g
matplot(t(BETA[1, , 1:RS]), type="l", ylab="parameter")
matplot(t(BETA[100, , 1:RS]), type="l", ylab="parameter")
matplot(t(BETA[1000, , 1:RS]), type="l", ylab="parameter")
matplot(t(BETA[2500, , 1:RS]), type="l", ylab="parameter")
matplot(t(BETA[5000, , 1:RS]), type="l", ylab="parameter")
matplot(t(THETA[, 1, 1:RS]), type="l", ylab="parameter")
matplot(t(THETA[, 5, 1:RS]), type="l", ylab="parameter")
matplot(t(THETA[, 10, 1:RS]), type="l", ylab="parameter")
matplot(THETA[1:RS, 6:9], type="l", ylab="parameter")


##�K�w���f���̉�A�W���̃p�����[�^
round(rbind(apply(THETA[, , burnin:RS], c(1, 2), mean), thetat), 3)
round(matrix(apply(THETA[, , burnin:(R/keep)], c(1, 2), function(x) quantile(x, 0.05)), nrow=k2, ncol=k1), 3)
round(matrix(apply(THETA[, , burnin:(R/keep)], c(1, 2), function(x) quantile(x, 0.95)), nrow=k2, ncol=k1), 3)

##�l�ʂ̃p�����[�^
i <- 20; sum(ID$id==i)   #�lid�𒊏o
round(cbind(beta_mu <- apply(BETA[, , burnin:RS], c(1, 2), mean), betat), 2)   #�l�ʂ̃p�����[�^����l�̎��㕽��
round(apply(BETA[, , burnin:RS], c(1, 2), summary), 3)   #�l�ʂ̃p�����[�^����l�̗v�񓝌v
round(apply(BETA[, , burnin:RS], c(1, 2), function(x) quantile(x, c(0.05, 0.95))), 3)   #����M�p���

#���ʂ��v���b�g
hist(BETA[i, 1, burnin:RS], col="grey", xlab="beta", main="beta�̌l���̎��㕪�z", breaks=20)
hist(BETA[, 3, RS], col="grey", xlab="beta", main="beta�̌l�ʂ̎��㕪�z", breaks=20)

##����\�����z�ōw���m����\��
logit.pred <- as.numeric((x * beta_mu[u_id, ]) %*% rep(1, k1))   #���W�b�g�̌v�Z
prob.pred <- as.numeric(exp(logit.pred) / (1 + exp(logit.pred)))   #�m���̌v�Z
summary(prob.pred)   #����\�����z�̗v��
hist(prob.pred, col="grey", xlab="�\���m��", main="�l�ʂ̎���\�����z", breaks=25)