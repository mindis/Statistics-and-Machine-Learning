#####�K���}��A���f��#####
options(warn=0)
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(2506787)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 100000   #�T���v����
k <- 11   #�����ϐ���


##�f���x�N�g���𐶐�
k1 <- 4; k2 <- 5; k3 <- 5
x1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
x2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hh, 1, pr)
}
x3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #�f�[�^������
k <- ncol(x)   #�����ϐ���

##�����ϐ��̐���
repeat {
  #�p�����[�^�̐���
  beta <- betat <- c(0.75, rnorm(k-1, 0, 0.5))
  alpha <- alphat <- 0.5
  
  #�ؒf�|�A�\�����z���牞���ϐ��𐶐�
  lambda <- as.numeric(exp(x %*% beta))   #���Ғl
  y <- rgamma(hh, lambda, alpha)
  
  if(max(y) > 15 & max(y) < 50){
    break
  }
}
hist(y, breaks=25, col="grey", main="�������Ԃ̕��z", xlab="��������")


####�Ŗޖ@�ŃK���}��A���f���𐄒�####
##�K���}��A���f���̐���̂��߂̊֐�
#�K���}��A���f���̑ΐ��ޓx
loglike <- function(theta, y, y_log, x){
  #�p�����[�^�̐ݒ�
  alpha <- theta[1]
  beta <- theta[-1]
  lambda <- as.numeric(exp(x %*% beta))   #���Ғl

  #�ΐ��ޓx�̘a
  LL <- sum(alpha * as.numeric(-y/lambda - x %*% beta) + alpha*log(alpha) - lgamma(alpha) + (alpha-1)*y_log)
  return(LL)
}

#�K���}��A���f���̑ΐ��ޓx�̔����֐�
dloglike <- function(theta, y, y_log, x){ 
  #�p�����[�^�̐ݒ�
  alpha <- theta[1]
  beta <- theta[-1]
  mu <- as.numeric(x %*% beta)
  lambda <- exp(mu)   #���Ғl
  
  #���z�x�N�g���̌v�Z
  sc1 <- hh*(log(alpha) - digamma(alpha)) + sum(1 - y/lambda + log(y/lambda))   #�`��p�����[�^�̌��z�x�N�g��
  sc2 <- colSums((y-lambda) / (lambda^2/alpha) * lambda * x)   #�ړx�p�����[�^�̌��z�x�N�g��
  sc <- c(sc1, sc2)
  return(sc)
}


##�K���}��A���f�������j���[�g���@�ōŖސ���
#�f�[�^�̐ݒ�
y_log <- log(y)   #y�̑ΐ�

#�p�����[�^�𐄒�
theta <- c(1.0, rep(0, k))   #�����l
res <- optim(theta, loglike, gr=dloglike, y, y_log, x, method="BFGS", hessian=TRUE,   #���j���[�g���@
             control=list(fnscale=-1, trace=TRUE))

#���茋��
LLmle <- res$value
theta <- res$par
alpha_mle <- theta[1]; beta_mle <- theta[-1]
c(alpha_mle, alphat)
rbind(beta_mle, betat)   #�ړx�p�����[�^
(tval <- theta/sqrt(-diag(solve(res$hessian))))   #t�l
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(hh)*length(beta)) #BIC

#�ϑ����ʂƊ��Ғl�̔�r
lambda <- as.numeric(exp(x %*% beta_mle))
round(data.frame(y, lambda), 3)   #�ϑ����ʂƂ̔�r

##�|�A�\����A�Ƃ̔�r
out <- glm(y ~ x[, -1], family=Gamma(link="log"))
rbind(mle1=beta_mle, mle2=as.numeric(out$coefficients))   #��A�W��
res$value; as.numeric(logLik(out))   #�ΐ��ޓx
sum((y - as.numeric(out$fitted.values))^2)   #���덷


####�n�~���g�j�A�������e�J�����@�ŃK���}��A���f���𐄒�####
##�ΐ����㕪�z���v�Z����֐�
loglike <- function(beta, alpha, y, y_log, x){
  
  #�K���}��A���f���̑ΐ��ޓx
  lambda <- as.numeric(exp(x %*% beta))   #���Ғl
  Lho <- sum(alpha * as.numeric(-y/lambda - x %*% beta) + alpha*log(alpha) - lgamma(alpha) + (alpha-1)*y_log)   #�ΐ��ޓx�֐�
  
  #���ϗʐ��K���z�̑ΐ����O���z
  log_mvn <- -1/2 * as.numeric(beta %*% inv_tau %*% beta)
  
  #�ΐ����㕪�z
  LL <- Lho + log_mvn
  return(list(LL=LL, Lho=Lho))
}

##HMC�Ŏړx�p�����[�^���T���v�����O���邽�߂̊֐�
#�K���}��A�̑ΐ����㕪�z�̔����֐�
dloglike_beta <- function(beta, alpha, y, y_log, x){ 
  
  #���Ғl�̐ݒ�
  mu <- as.numeric(x %*% beta)
  lambda <- exp(mu)   #���Ғl
  
  #�����֐��̐ݒ�
  dlgamma <- colSums((y-lambda) / (lambda^2/alpha) * lambda * x)   #�ړx�p�����[�^�̌��z�x�N�g��
  dmvn <- as.numeric(-inv_tau %*% beta)
  
  #�ΐ����㕪�z�̔����֐��̘a
  LLd <- -(dlgamma + dmvn)
  return(LLd)
}

#���[�v�t���b�O�@�������֐�
leapfrog_beta <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, alpha, y, y_log, x) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, alpha, y, y_log, x) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##HMC�Ō`��p�����[�^���T���v�����O���邽�߂̊֐�
#�K���}��A���f���̑ΐ����㕪�z�̔����֐�
dloglike_alpha <- function(alpha, beta, y, y_log, x){ 
  #���Ғl�̐ݒ�
  mu <- as.numeric(x %*% beta)
  lambda <- exp(mu)   #���Ғl
  
  #���z�x�N�g���̌v�Z
  dlgamma <- hh*(log(alpha) - digamma(alpha)) + sum(1 - y/lambda + log(y/lambda))   #�`��p�����[�^�̌��z�x�N�g��
  return(dlgamma)
}

#���[�v�t���b�O�@�������֐�
leapfrog_alpha <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, beta, y, y_log, x) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, beta, y, y_log, x) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}


##�A���S���Y���̐ݒ�
R <- 5000
keep <- 2
disp <- 20
burnin <- 1000/keep
iter <- 0
e <- 0.001
L <- 3

#���O���z�̐ݒ�
gamma <- rep(0, k)
inv_tau <- solve(100 * diag(k))

#�p�����[�^�̐^�l
alpha <- alphat
beta <- betat

#�����l�̐ݒ�
alpha <- 1.0
beta <- rep(0, k)


#�p�����[�^�̊i�[�p�z��
ALPHA <- rep(0, R/keep)
BETA <- matrix(0, nrow=R/keep, ncol=k)

#�ΐ��ޓx�֐��̊�l
LLbest <- loglike(betat, alphat, y, y_log, x)$Lho


####HMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�ړx�p�����[�^���T���v�����O
  #HMC�̐V�����p�����[�^�𐶐�
  rold <- as.numeric(mvrnorm(1, rep(0, k), diag(k)))   #�W�����ϗʐ��K���z����p�����[�^�𐶐�
  betad <- beta
  
  #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
  res <- leapfrog_beta(rold, betad, dloglike_beta, e, L)
  rnew <- res$r
  betan <- res$z
  
  #�ړ��O�ƈړ���̃n�~���g�j�A��
  Hnew <- -loglike(betan, alpha, y, y_log, x)$LL + as.numeric(rnew^2 %*% rep(1, k))/2
  Hold <- -loglike(betad, alpha, y, y_log, x)$LL + as.numeric(rold^2 %*% rep(1, k))/2
  
  #�p�����[�^�̍̑�������
  rand <- runif(1)   #��l���z���痐���𔭐�
  gamma <- min(c(1, exp(Hold - Hnew)))   #�̑𗦂�����
  gamma_beta <- gamma
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- as.numeric(gamma > rand)
  beta <- flag*betan + (1-flag)*betad
  
  
  ##�`��p�����[�^���T���v�����O
  #MH�@�̐V�����p�����[�^�𐶐�
  d <- sign(dloglike_alpha(alpha, beta, y, y_log, x))   #���z�̕���
  alphad <- alpha
  alphan <- alphad + d*abs(rnorm(1, 0, 0.01))  

  #�Ɨ�MH�@�̑ΐ��ޓx
  lognew <- loglike(beta, alphan, y, y_log, x)$Lho
  logold <- loglike(beta, alphad, y, y_log, x)$Lho 
  
  #�p�����[�^�̍̑�������
  rand <- runif(1)   #��l���z���痐���𔭐�
  gamma <- min(c(1, exp(lognew - logold)))   #�̑𗦂�����
  gamma_alpha <- gamma
  
  #gamma�̒l�Ɋ�Â��V����alpha���̑����邩�ǂ���������
  flag <- as.numeric(gamma > rand)
  alpha <- flag*alphan + (1-flag)*alphad
  
  ##�T���v�����O���ʂ̕ۑ��ƕ\��
  #�T���v�����O���ʂ̕ۑ�
  if(rp%%keep==0){
    mkeep <- rp/keep
    ALPHA[mkeep] <- alpha
    BETA[mkeep, ] <- beta
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���Z�o
    LL <- loglike(beta, alpha, y, y_log, x)$Lho
    
    #�T���v�����O���ʂ�\��
    print(rp)
    print(c(gamma_alpha, gamma_beta))
    print(c(LL, LLmle, LLbest))
    print(round(c(alpha, alpha_mle, alphat), 3))
    print(round(rbind(beta=beta, betaml=beta_mle, betat=betat), 3))
  }
}

####���茋�ʂ̊m�F�Ɨv��####
#�T���v�����O���ʂ̃v���b�g
plot(1:(R/keep), ALPHA, type="l", main="alpha�̃T���v�����O���ʂ̃v���b�g", ylab="alpha�̐���l", xlab="�T���v�����O��")
matplot(BETA, type="l", main="beta�̃T���v�����O���ʂ̃v���b�g", ylab="beta�̐���l", xlab="�T���v�����O��")