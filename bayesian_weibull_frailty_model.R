#####�J��Ԃ��̂��鐶�����ԃ��f��(���d�C�x���g���f��)#####
library(MASS)
library(survival)
library(bayesm)
library(MCMCpack)
library(frailtypack)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

##�f�[�^�̐ݒ�
hh <- 500
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

####�����ϐ��̐ݒ�####
##�p�����[�^�̐ݒ�
for(i in 1:10000){
  print(i)
  alpha <- runif(1, 0.8, 1.2)   #�`��p�����[�^
  b1 <- c(runif(1, 0, 1.2), runif(k1/2, 0, 0.7), runif(k1/2, -0.4, 0.9), runif(k2+k3, -0.5, 1.0))   #�Œ���ʂ̃p�����[�^
  u <- rnorm(hh, 0, 0.5)   #�t���C���e�B�[�̃p�����[�^
  
  #���C�u�����z����C�x���g���Ԃ𔭐�
  lambda <- exp(X %*% b1 + Z %*% u)   #���`����
  y <- rweibull(nrow(lambda), alpha, lambda)
  
  if(min(y) > 0.01) break
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
    y_vec <- c(y_ind[index1], dt-c_sum[length(index1)])
  } else if(max(c_sum) > 100 & length(index1)==0) {
    y_vec <- 100
    z <- 1
  } else {
    y_vec <- y_ind[index2]
    z[length(y_vec)] <- 1
  }

  #�ł��؂�ꂽ�ϐ����i�[
  y.list[[i]] <- y_vec[index2]
  ID.list[[i]] <- ID[id==i, ][index2, ]
  X.list[[i]] <- X[id==i, ][index2, ]
  Z.list[[i]] <- Z[id==i, ][index2, ]
  z.list[[i]] <- z[index2]
}

#���X�g���s�񂠂邢�̓x�N�g����
y <- unlist(y.list)
ID <- do.call(rbind, ID.list)
no <- 1:nrow(ID)
X <- do.call(rbind, X.list)
Z <- do.call(rbind, Z.list)
z <- unlist(z.list)

#�f�[�^�̊m�F�Ɖ���
round(data.frame(no, ID[, 2:3], y, z, x=X, z=Z[, 1:10]), 2)
hist(y, col="grey", breaks=30, main="�C�x���g����", xlab="�o�ߎ���")
table(z)   #�ł��؂萔


####�}���R�t�A�������e�J�����@�Ń��C�u�����L�t���C���e�B���f���𐄒�####
##���C�u�����n�U�[�h���f���̑ΐ��ޓx
loglike <- function(alpha, beta, v, y, X, Z, z){
  lambda <- exp(X %*% beta + Z %*% v)   #���`����
  LL <- sum(z*(log(lambda)+log(alpha)+(alpha-1)*log(y)) - lambda*y^alpha)   #�ΐ��ޓx���v�Z
  return(LL)
}

##�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- c()   #�ΐ��ޓx�̕ۑ��p

##���O���z�̐ݒ�
#�Œ���ʂ̎��O���z
betas <- rep(0, ncol(X))   #��A�W���̕��ς̎��O���z
sigma <- diag(rep(0.01), ncol(X))   #��A�W���̎��O���z�̕��U

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
BETA <- matrix(0, nrow=R/keep, ncol=ncol(X))
ALPHA <- matrix(0, nrow=R/keep, ncol=1)
RANDOM <- matrix(0, nrow=R/keep, ncol=hh)

##�����l�̐ݒ�
oldalpha <- runif(1, 0.6, 1.5)   #�`��p�����[�^
oldbetas <- c(runif(1, 0, 1.4), runif(k1/2, 0, 1.0), runif(k1/2, -1.0, 1.0), runif(k2+k3, -1.0, 1.0))   #�Œ���ʂ̃p�����[�^ 
betas.random <- rnorm(hh, 0, 0.75)

##�����_���E�H�[�N�̕��U
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
  beta0 <- c(0, 1, runif(ncol(X)-1, -0.5, 0.5))
  
  #���j���[�g���@�őΐ��ޓx���ő剻
  res <- try(optim(beta0, llike, y=y, X=X, z=z, method="BFGS", 
                   hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #�G���[����
}

rw_alpha <- sqrt(-diag(solve(res$hessian))[1])
rw_beta <- 0.5*diag(-diag(solve(res$hessian))[2:length(res$par)])


####MCMC�Ń��C�u�����L�t���C���e�B���f���𐄒�####
##�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O

##MH�T���v�����O�ŌŒ���ʂ��T���v�����O
betad <- oldbetas
betan <- betad + mvrnorm(1, rep(0, ncol(X)), rw_beta)

#�ΐ��ޓx�Ƒΐ����O���z���v�Z
lognew1 <- loglike(oldalpha, betan, betas.random, y, X, Z, z)
logold1 <- loglike(oldalpha, betad, betas.random, y, X, Z, z)
logpnew1 <- lndMvn(betan, betas, sigma)
logpold1 <- lndMvn(betad, betas, sigma)

#MH�T���v�����O
alpha1 <- min(1, exp(lognew1 + logpnew1 - logold1 - logpold1))
if(alpha == "NAN") alpha1 <- -1

#��l�����𔭐�
u <- runif(1)

#u < alpha�Ȃ�V�����Œ����beta���̑�
if(u < alpha1){
  oldbetas <- betan
  logl <- lognew1
  
  #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
} else {
  logl <- logold1
}

##MH�T���v�����O�Ō`��p�����[�^���T���v�����O