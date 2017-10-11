#####�p�����ԂƗ��U���Ԃ̓������͂̂��߂̊K�w�x�C�Y���f��#####
library(MASS)
library(survival)
library(frailtypack)
library(nlme)
library(glmm)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

##�f�[�^�̐ݒ�
hh <- 1500
pt <- rpois(hh, 15.0)
pt <- ifelse(pt==0, 1, pt)
hhpt <- sum(pt)
dt <- 100

##ID�̐ݒ�
id <- rep(1:hh, pt)
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}
ID <- data.frame(no=1:hhpt, id=id, t=t)

####�����ϐ��̐ݒ�####
##�Œ���ʂ̐����ϐ��̐ݒ�
#�l���ł̋��ʕϐ��̔���
k1 <- 4
X1 <- matrix(0, nrow=hhpt, ncol=k1)
for(i in 1:hh){
  X1[ID$id==i, 1:2] <- matrix(rnorm(2, 0, 1), nrow=sum(ID$id==i), ncol=2, byrow=T)
  X1[ID$id==i, 3:4] <- matrix(rbinom(2, 1, runif(1, 0.4, 0.6)), nrow=sum(ID$id==i), ncol=2, byrow=T)
}

#�l�A���_�ŕω�����A���ϐ��̔���
k2 <- 3
X2 <- matrix(runif(hhpt*(k2), 0, 1), hhpt, (k2))

#�l�A���_�ŕω������l�ϐ�
k3 <- 3
X3 <- matrix(0, hhpt, k3)
for(i in 1:k3){
  bin <- rbinom(hhpt, 1, runif(1, 0.3, 0.7))
  X3[, i] <- bin
}

#�f�[�^�̌���
X <- cbind(1, X1, X2, X3)

##�ϗʌ��ʂ̐����ϐ��̐ݒ�
k <- 1   #�ϗʌ��ʂ̕ϐ���
Z <- matrix(0, nrow=hhpt, ncol=hh*k)
for(i in 1:hh){
  r <- ((i-1)*k+1):((i-1)*k+k)
  
  Z[ID$id==i, r] <- 1
}

####�p�����Ԃ̉����ϐ��̐ݒ�####
##�p�����[�^�̐ݒ�
for(i in 1:10000){
  print(i)
  alpha <- runif(1, 0.8, 1.2)   #�`��p�����[�^
  b1 <- c(runif(1, 0, 1.2), runif(k1/2, 0, 0.7), runif(k1/2, -0.4, 0.9), runif(k2+k3, -0.5, 1.0))   #�Œ���ʂ̃p�����[�^
  v.par <- runif(1, 0.6, 1.0)
  f <- rnorm(hh, 0, v.par)   #�t���C���e�B�[�̃p�����[�^
 
  #���C�u�����z����C�x���g���Ԃ𔭐�
  lambda <- exp(X %*% b1 + Z %*% f)   #���`����
  y <- rweibull(nrow(lambda), alpha, lambda)
  
  if(min(y) > 0.001) break
}

##�ł��؂�̐ݒ�
#�ϐ��̊i�[�p���X�g
ID.list <- list()
y.list <- list()
X.list <- list()
Z.list <- list()
z.list <- list()

#�l���Ƃɑł��؂�ϐ���ݒ�
for(i in 1:hh){
  print(i)
  y_ind <- y[id==i]
  z <- rep(0, length(y_ind))
  c_sum <- cumsum(y_ind)
  
  #�ݐώ��Ԃ�100�ȏ�̃C�x���g�͑ł��؂�
  index1 <- subset(1:length(c_sum), c_sum <= 100)
  
  if(max(c_sum) <= 100){
    index2 <- index1
  } else {
    index2 <- c(index1, length(index1)+1)
  }
  
  #�����ϐ��̑ł��؂��ݒ�
  if(max(c_sum) > 100 & length(index1) > 0){
    print(1)
    y_vec <- c(y_ind[index1], dt-c_sum[length(index1)])
    z[length(y_vec)] <- 1
  } else if(max(c_sum) > 100 & length(index1)==0) {
    print(2)
    y_vec <- 100
    z <- 1
  } else {
    print(3)
    y_vec <- y_ind[index2]
  }
  
  #�ł��؂�ꂽ�ϐ����i�[
  y.list[[i]] <- y_vec[index2]
  ID.list[[i]] <- ID[id==i, ][index2, ]
  X.list[[i]] <- X[id==i, ][index2, ]
  Z.list[[i]] <- Z[id==i, ][index2, ]
  z.list[[i]] <- z[index2]
}

#���X�g���s�񂠂邢�̓x�N�g����
y1 <- unlist(y.list)
ID <- do.call(rbind, ID.list)
no <- 1:nrow(ID)
Data1 <- do.call(rbind, X.list)
Z1 <- do.call(rbind, Z.list)
z1 <- 1-unlist(z.list)

#�f�[�^�̊m�F�Ɖ���
hist(y1, col="grey", breaks=30, main="�C�x���g����", xlab="�o�ߎ���")
table(z)   #�ł��؂萔


####���U�I���̉����ϐ��̔���####
##�����ϐ��̐ݒ�
#�p�����Ԃ̕ϗʌ��ʂ�����ϐ��Ƃ��Đݒ�
u_vec <- c()
index_id <- list()

for(i in 1:hh){
  num <- max(ID$t[ID$id==i])
  index_id[[i]] <- subset(1:nrow(ID), ID$id==i)
  u_vec <- c(u_vec, rep(f[i], num))
}

#�f�[�^�̌���
Data2 <- cbind(1, Data1[, -1], u_vec)
colnames(Data2) <- 1:ncol(Data2)

##���W�X�e�B�b�N��A���牞���ϐ��𔭐�
#�p�����[�^�̐ݒ�
theta00 <- runif(1, -1.2, -0.8)
theta01 <- c(runif(k1/2, 0, 0.8), runif(k1/2, -0.7, 0.9))
theta02 <- c(runif(k2, 0, 0.8), runif(k3, -0.8, 1.1))
theta03 <- runif(1, 0.5, 0.8)
theta0 <- c(theta00, theta01, theta02, theta03)

#���W�b�g�Ɗm���̌v�Z
logit <- Data2 %*% theta0
Pr0 <- exp(logit)/(1+exp(logit))

#�x���k�[�C�������牞���ϐ��𔭐�
y2 <- rbinom(length(Pr0), 1, Pr0)
round(cbind(y1, y2, Pr0, ID, Data2), 3)   #�����������f�[�^�̊m�F


####�}���R�t�A�������e�J�����@�Ōp�����ԂƗ��U�I���𓯎�����####
##���C�u�����n�U�[�h���f���̑ΐ��ޓx
loglike1 <- function(alpha, beta, v, y, X, Z, z){
  lambda <- exp(X %*% beta + Z %*% v)   #���`����
  LLi <- z*(log(lambda)+log(alpha)+(alpha-1)*log(y)) - lambda*y^alpha   #�ΐ��ޓx���v�Z
  LL <- sum(LLi)
  LL.val <- list(LL=LL, LLi=LLi)
  return(LL.val)
}

loglike2 <- function(theta, y, X){
  logit <- X %*% theta
  Pr <- exp(logit)/(1+exp(logit))
  Li <- y * log(Pr) + (1-y) * log(1-Pr)
  LL <- sum(Li)
  return(LL)
}

##�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- c()   #�ΐ��ޓx�̕ۑ��p

##���O���z�̐ݒ�
#�Œ���ʂ̎��O���z
betas <- rep(0, ncol(Data1))   #��A�W���̕��ς̎��O���z
sigma1 <- diag(0.01, ncol(Data1))   #��A�W���̎��O���z�̕��U

thetas <- rep(0, ncol(Data2))
sigma2 <- diag(0.01, ncol(Data2))

#�`��p�����[�^�̎��O���z
alpha_mu <- 0
alpha_sigma <- 2.5

#�ϗʌ��ʂ̎��O���z
Deltabar <- 0
Adelta <- 100   
tau1 <- 1   #�t�K���}���z�̌`��p�����[�^
tau2 <- 0.01   #�t�K���}���z�̃X�P�[���p�����[�^
beta.random <- rep(0, nrow=hh)   #�ϗʌ��ʂ̎��O���z�̕��ς�0�ɌŒ�

##�T���v�����O���ʂ̕ۑ��p�z��
BETA <- matrix(0, nrow=R/keep, ncol=ncol(Data1))
ALPHA <- matrix(0, nrow=R/keep, ncol=1)
RANDOM <- matrix(0, nrow=R/keep, ncol=hh)
THETA <- matrix(0, nrow=R/keep, ncol(Data2))
SIGMA <- matrix(0, nrow=R/keep, ncol=1)

##�����l�̐ݒ�
oldalpha <- runif(1, 0.6, 1.5)   #�`��p�����[�^
oldbetas <- c(runif(1, 0, 1.4), runif(k1/2, 0, 1.0), runif(k1/2, -1.0, 1.0), runif(k2+k3, -1.0, 1.0))   #�Œ���ʂ̃p�����[�^ 
oldtheta <- rep(0, ncol(Data2))
betas.random <- rnorm(hh, 0, 0.75)
cov.random <- 1.0


##�����_���E�H�[�N�̕��U��ݒ�
#���C�u�����n�U�[�h���f���̑ΐ��ޓx
llike <- function(theta, y, X, z){
  a <- exp(theta[1])
  beta <- theta[2:(ncol(X)+1)]
  
  lambda <- exp(X %*% beta)   #���`����
  LL <- sum(z*(log(lambda)+log(a)+(a-1)*log(y)) - lambda*y^a)   #�ΐ��ޓx���v�Z
  return(LL)
}


#���j���[�g���@�őΐ��ޓx���ő剻
for(i in 1:1000){
  print(i)
  #�����p�����[�^�̐ݒ�
  beta01 <- c(0, 1, runif(ncol(Data1)-1, -0.5, 0.5))
  beta02 <- beta01[-1]
  
  #���j���[�g���@�őΐ��ޓx���ő剻
  res1 <- try(optim(beta01, llike, y=y1, X=Data1, z=z1, method="BFGS", 
                   hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)
  
  res2 <- optim(beta02, loglike2, y=y2, X=Data1, method="BFGS", hessian=TRUE, 
                control=list(fnscale=-1))
  if(class(res1) == "try-error") {next} else {break}   #�G���[����
}

#���C�u�����f���̃����_���E�H�[�N�̕��U
rw_alpha <- sqrt(-diag(solve(res1$hessian))[1])
rw_beta <- 0.5*diag(-diag(solve(res1$hessian))[2:length(res1$par)])

#���W�b�g���f���̃����_���E�H�[�N�̕��U
rw_theta <- 0.5*diag(c(-diag(solve(res2$hessian)), -solve(res2$hessian)[2, 2]))


####MCMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##MH�T���v�����O�ŌŒ���ʂ��T���v�����O
  betad <- oldbetas
  betan <- betad + 0.5 * mvrnorm(1, rep(0, ncol(Data1)), rw_beta)
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew1 <- loglike1(oldalpha, betan, betas.random, y1, Data1, Z1, z1)$LL
  logold1 <- loglike1(oldalpha, betad, betas.random, y1, Data1, Z1, z1)$LL
  logpnew1 <- lndMvn(betan, betas, sigma1)
  logpold1 <- lndMvn(betad, betas, sigma1)
  
  #MH�T���v�����O
  alpha1 <- min(1, exp(lognew1 + logpnew1 - logold1 - logpold1))
  if(alpha == "NAN") alpha1 <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V�����Œ����beta���̑�
  if(u < alpha1){
    oldbetas <- betan
    logl1 <- lognew1
    
    #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
  } else {
    logl1 <- logold1
  }
  
  rate1 <- as.numeric((oldbetas!=betad)[1])
  
  
  ##MH�T���v�����O�Ō`��p�����[�^���T���v�����O
  alphad <- abs(oldalpha)
  alphan <- abs(alphad + rnorm(1, 0, 0.1))
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew2 <- loglike1(alphan, oldbetas, betas.random, y1, Data1, Z1, z1)$LL
  logold2 <- loglike1(alphad, oldbetas, betas.random, y1, Data1, Z1, z1)$LL
  logpnew2 <- -1/2 * alpha_sigma^(-1) * (log(alphan) - alpha_mu)^2
  logpold2 <- -1/2 * alpha_sigma^(-1) * (log(alphad) - alpha_mu)^2
  
  #MH�T���v�����O
  alpha2 <- min(1, exp(lognew2 + logpnew2 - logold2 - logpold2))
  if(alpha2 == "NAN") alpha2 <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V�����Œ����beta���̑�
  if(u < alpha2){
    oldalpha <- alphan
    
    #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
  } else {
    oldalpha <- alphad
  }
  
  ##MH�T���v�����O�Ńt���C���e�B�p�����[�^���T���v�����O
  betad.random <- betas.random
  betan.random <- betad.random + rnorm(hh, 0, 0.2)
  
  #���O���z�̌덷���v�Z
  er.new <- betan.random - beta.random
  er.old <- betad.random - beta.random
  inv.cov <- cov.random^-1
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew3 <- loglike1(oldalpha, oldbetas, betan.random, y1, Data1, Z1, z1)$LLi
  logold3 <- loglike1(oldalpha, oldbetas, betad.random, y1, Data1, Z1, z1)$LLi
  logpnew3 <- -0.5 * inv.cov * er.new^2
  logpold3 <- -0.5 * inv.cov * er.old^2
  
  #ID�ʂɑΐ��ޓx�̘a�����
  lognew.ind <- as.matrix(data.frame(logl=lognew3, id=ID$id) %>%
                            dplyr::group_by(id) %>%
                            dplyr::summarize(sum=sum(logl)))[, 2]
  
  logold.ind <- as.matrix(data.frame(logl=logold3, id=ID$id) %>%
                            dplyr::group_by(id) %>%
                            dplyr::summarize(sum=sum(logl)))[, 2]
  
  #MH�T���v�����O
  rand <- runif(hh)
  LLind.diff <- exp(lognew.ind + logpnew3 - logold.ind - logpold3)   #���p�����v�Z
  alpha3 <- ifelse(LLind.diff > 1, 1, LLind.diff)
  betas.random <- ifelse(alpha3 > rand, betan.random, betad.random)   #alpha��rand�������Ă�����̑�
  rate3 <- sum(betas.random==betad.random)/hh
  
  ##�t�K���}���z���番�U�������T���v�����O
  shape <- hh+tau1
  scale <- sum((betas.random - beta.random)^2)+tau2
  cov.random <- rinvgamma(1, shape, scale)
  
  
  ##���U�I�����f���̉�A�W�����T���v�����O
  #�ϗʌ��ʂ̃f�[�^�̐ݒ�
  w_vec <- rep(0, nrow(Data1))
  for(i in 1:hh){
    w_vec[index_id[[i]]] <- betas.random[i]
  }
  
  #�f�[�^�̌���
  Data2 <- cbind(Data1, w_vec)
  
  ##MH�T���v�����O�ŉ�A�W�����T���v�����O
  thetad <- oldtheta
  thetan <- thetad + 0.25 * mvrnorm(1, rep(0, ncol(Data2)), rw_theta)
  
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew4 <- loglike2(thetan, y2, Data2)
  logold4 <- loglike2(thetad, y2, Data2)
  logpnew4 <- lndMvn(thetan, thetas, sigma2)
  logpold4 <- lndMvn(thetad, thetas, sigma2)
  
  #MH�T���v�����O
  alpha4 <- min(1, exp(lognew4 + logpnew4 - logold4 - logpold4))
  if(alpha4 == "NAN") alpha4 <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V�����Œ����beta���̑�
  if(u < alpha4){
    oldtheta <- thetan
    logl2 <- lognew4
    
    #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
  } else {
    logl2 <- logold4
  }
  
  BETA <- matrix(0, nrow=R/keep, ncol=ncol(Data1))
  ALPHA <- matrix(0, nrow=R/keep, ncol=1)
  RANDOM <- matrix(0, nrow=R/keep, ncol=hh)
  THETA <- matrix(0, nrow=R/keep, ncol(Data2))
  SIGMA <- matrix(0, nrow=R/keep, ncol=1)
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbetas
    ALPHA[mkeep, ] <- oldalpha
    THETA[mkeep, ] <- oldtheta
    RANDOM[mkeep, ] <- betas.random
    SIGMA[mkeep, ] <- cov.random
    
    print(rp)
    print(round(logl1 + logl2, 2))
    print(round(rbind(c(oldalpha, oldbetas), c(alpha, b1)), 2))
    print(round(rbind(oldtheta, theta0), 2))
    print(round(c(sqrt(cov.random), v.par), 3))
    print(round(c(rate1, rate3), 3))
  }
}

