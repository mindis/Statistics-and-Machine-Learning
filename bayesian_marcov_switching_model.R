#####�x�C�W�A����ʏ�ԋ�ԃ}���R�t�؊������f��#####
library(MASS)
library(MSwM) 
library(matrixStats)
library(FAdist)
library(Rfast)
library(bayesm)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(18478)

####�f�[�^�̔���####
n <- 3000   #�ϑ����Ԑ�
s <- 4   #�؊�����
k1 <- 2   #�ϑ����f���̉����ϐ��̃p�����[�^��
k2 <- 2   #�ϑ����f���̐����ϐ��̃p�����[�^��

####�p�����[�^�̐ݒ�####
##�}���R�t�؊����s��̃p�����[�^��ݒ�
#�����m���̐ݒ�
p0 <- c(0.5, 0.2, 0.2, 0.1)

#���ڍs��̐ݒ�
pr1 <- c(0.7, 0.3, 0, 0)
pr2 <- c(0, 0, 0.6, 0.4)
pr3 <- c(0.75, 0.25, 0, 0)
pr4 <- c(0, 0, 0.8, 0.2)
Pr0 <- rbind(pr1, pr2, pr3, pr4)

##�ϑ����f���̃p�����[�^��ݒ�
#�W�������i�̃p�����[�^
tau01 <- c(8.0, 7.0, 9.0, 7.5)
tau02 <- c(0.5, 3.5, 1.0, 2.25)
tau0 <- cbind(tau01, tau02)
colnames(tau0) <- c("", "")


#����_���̃p�����[�^
beta00 <- c(2.1, 3.3, 1.5, 2.7)
beta01 <- c(1.0, 1.8, 0.8, 1.3)
beta02 <- c(0.5, 0.7, 0.3, 0.6)
beta0 <- rbind(beta00, beta01, beta02)
sigma0 <- 0.5


####�����ϐ��̔���####
#�f�[�^�̕ۑ��p�z��
Z0 <- matrix(0, nrow=n, ncol=s)
z0 <- rep(0, n)
y <- matrix(0, nrow=n, ncol=k1)
X <- cbind(1, matrix(0, nrow=n, ncol=k2))
X[, ncol(X)] <- rbinom(n, 1, 0.4)

#�����l�̐ݒ�
Z0[1, ] <- extraDistr::rmnom(1, 1, p0)
seg <- z0[1] <- Z0[1, ] %*% 1:s
y[1, 2] <- rbeta(1, tau01[seg], tau02[seg])
X[1, 2] <- y[1, 2]
y[1, 1] <- X[1, ] %*% beta0[, seg] + rnorm(1, 0, sigma0)

##���Ԃ��Ƃɒ����I�ɉ����ϐ��𔭐�������
for(i in 2:n){
  print(i)
  #�}���R�t���ڍs�񂩂���ݕϐ��̊����𐶐�
  obs_seg <- z0[i-1]
  Z0[i, ] <- extraDistr::rmnom(1, 1, Pr0[obs_seg, ])
  seg <- z0[i] <- Z0[i, ] %*% 1:s
  
  #���ݕϐ�����ϑ��f�[�^�𔭐�
  y[i, 2] <- rbeta(1, tau01[seg], tau02[seg])
  X[i, 2] <- y[i, 2]
  y[i, 1] <- X[i, ] %*% beta0[, seg] + rnorm(1, 0, sigma0)
}
table(z0)


####�}���R�t�A�������e�J�����@�Ń}���R�t�X�C�b�`���O���f���𐄒�####
##�A���S���Y���̐ݒ�
R <- 20000  
keep <- 4
iter <- 0
sbeta <- 1.5

##���O���z�̐ݒ�
#���`��A���f���̎��O���z
alpha0 <- rep(0, k2+1)
cov0 <- diag(0.01, k2+1)
s0 <- 0.01
v0 <- 0.01

#�x�[�^���z���f���̎��O���z
lambda0 <- rep(0, 2)
omega0 <- diag(0.01, 2)

#�}���R�t���ڍs��̎��O���z
phi01 <- rep(1, 2)
phi02 <- rep(1, s)

##�����l�̐ݒ�
#���`��A���f���̏����l
oldbeta <- matrix(solve(t(X) %*% X) %*% t(X) %*% y[, 1], nrow=k1+1, ncol=s)
oldbeta[1, ] <- c(1.5, 2.5, 1, 2)
oldsigma <- sd(rep(y[, 1], s) - as.numeric(X %*% oldbeta))

#�x�[�^���z�̏����l
oldtau <- matrix(beta.mle(y[, 2])$param, nrow=s, ncol=2, byrow=T)
oldtau[, 2] <- c(0.5, 2.0, 0.75, 1.5)
rw <- rbind(c(0.005, 0.0025), matrix(c(0.01, 0.01), nrow=s-1, ncol=2, byrow=T))


#���ڍs��̏����l
Pr1 <- c(0.7, 0.3, 0, 0)
Pr2 <- c(0, 0, 0.7, 0.3)
Pr3 <- c(0.7, 0.3, 0, 0)
Pr4 <- c(0, 0, 0.7, 0.3)
Pr <- rbind(Pr1, Pr2, Pr3, Pr4)
r0 <- r <- rep(0.25, s)

##�C���f�b�N�X�̐ݒ�
index_s <- list()
for(i in 1:s){
  index_s[[i]] <- which(Pr[i, ] > 0)
}

##�p�����[�^�̊i�[�p�z��
BETA <- matrix(0, nrow=R/keep, ncol=(k1+1)*s)
SIGMA <- rep(0, R/keep)
LAMBDA <- matrix(0, nrow=R/keep, ncol=k2*s)
Zi <- matrix(0, nrow=R/keep, ncol=n)
PR <- array(0, dim=c(s, s, R/keep))
storage.mode(Zi) <- "integer"


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##���Ԃ��Ƃɓ��I�Ȑ��ݕϐ�z�𐶐�(�V�X�e�����f���̃p�����[�^���T���v�����O)
  #�f�[�^�̊i�[�p�z��
  Z <- matrix(0, nrow=n, ncol=s)
  z <- rep(0, n)
  z_rate <- matrix(0, nrow=n, ncol=s)
  LLi1 <- matrix(0, nrow=n, ncol=s)
  LLi2 <- matrix(0, nrow=n, ncol=s)
  
  #�ϑ����f���̖ޓx���v�Z
  mu <- X %*% oldbeta   #���`��A���f���̕��ύ\��
  for(j in 1:s){
    LLi1[, j] <- dnorm(y[, 1], mu[, j], oldsigma)
    LLi2[, j] <- dbeta(y[, 2], oldtau[j, 1], oldtau[j, 2])
  }
  LLi <- LLi1 * LLi2
  
  #�������z���1���ڂ̐��ݕϐ�z�𐶐�
  z_rate[1, ] <- r*LLi[1, ] / sum(r*LLi[1, ])
  Z[1, ] <- extraDistr::rmnom(1, 1, z_rate[1, ])
  z[1] <- Z[1, ] %*% 1:s
  
  #2���ڈȍ~�̐��ݕϐ�z�𒀎��I�Ɋ��蓖�Ă�
  for(i in 2:n){
    #�������z�����ݕϐ�z�𐶐�
    index <- index_s[[z[i-1]]]
    z_rate[i, index] <- Pr[z[i-1], index]*LLi[i, index] / sum(Pr[z[i-1], index]*LLi[i, index])
    Z[i, ] <- extraDistr::rmnom(1, 1, z_rate[i, ])
    z[i] <- Z[i, ] %*% 1:s
  }
  
  ##�}���R�t���ڍs��̃p�����[�^���T���v�����O
  index_z <- list()
  for(j in 1:s){
    index_z[[j]] <- which(z[1:(n-1)]==j) + 1  
    index <- index_z[[j]]   #1���O�̐��ݕϐ��̊��������
    par <- colSums(Z[index, index_s[[j]]])   #�x�[�^���z�̃p�����[�^
    pr <- rbeta(1, par[1]+phi01[1], par[2]+phi01[2])   #�x�[�^���z���p�����[�^���T���v�����O
    Pr[j, index_s[[j]]] <- c(pr, 1-pr)
  }
  r <- as.numeric(extraDistr::rdirichlet(1, colSums(Z) + phi02))   #�f�B�N�������z��荬�������T���v�����O
  
  
  ##�ϑ����f���̃p�����[�^���T���v�����O
  mu <- rep(0, n)   #��A���f���̕��ύ\���̊i�[�p�z��
  
  for(j in 1:s){
    #���ݕϐ��̊����𒊏o
    index <- which(z==j)
    
    #���`��A���f���̉�A�W�����T���v�����O
    x <- X[index, ]
    XXV <- solve(t(x) %*% x + cov0)
    Xy <- t(x) %*% y[index, 1]
    beta_mu <- XXV %*% (Xy + cov0 %*% alpha0)
    oldbeta[, j] <- mvrnorm(1, beta_mu, oldsigma*XXV)   #���ϗʐ��K���z����A�W�����T���v�����O
    mu[index] <- x %*% oldbeta[, j]
    
    #�x�[�^���z�̃p�����[�^���T���v�����O
    taud <- oldtau[j, ]
    taun <- abs(taud + mvrnorm(1, rep(0, 2), diag(rw[j, ], 2)))
    
    #�ΐ��ޓx�Ƒΐ����O���z���v�Z
    lognew <- sum(dbeta(y[index, 2], taun[1], taun[2], log=TRUE))
    logold <- sum(dbeta(y[index, 2], taud[1], taud[2], log=TRUE))
    logpnew <- lndMvn(taun, lambda0, omega0)
    logpold <- lndMvn(taud, lambda0, omega0)
    
    #MH�T���v�����O
    alpha <- min(1, exp(lognew + logpnew - logold - logpold))
    if(alpha == "NAN") alpha <- -1
    
    #��l�����𔭐�
    u <- runif(1)
    
    #u < alpha�Ȃ�V�����p�����[�^���̑�
    if(u < alpha){
      oldtau[j, ] <- taun
      
      #�����łȂ��Ȃ�p�����[�^���X�V���Ȃ�
    } else {
      oldtau[j, ] <- taud
    }
  }
  
  #���`��A���f���̋��ʂ̕��U���T���v�����O
  er <- y[, 1] - mu
  s1 <- s0 + t(er) %*% er
  v1 <- v0 + n
  oldsigma <- sqrt(1/(rgamma(1, v1/2, s1/2)))   #�t�K���}���z����sigma^2���T���v�����O
  
  
  ##�T���v�����O���ʂ̕ۑ��ƕ\��
  if(rp%%keep==0){
  
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    BETA[mkeep, ] <- as.numeric(oldbeta)
    SIGMA[mkeep] <- oldsigma
    LAMBDA[mkeep, ] <- as.numeric(oldtau)
    Zi[mkeep, ] <- z
    PR[, , mkeep] <- Pr
    
    #�T���v�����O���ʂ̕\��
    print(rp)
    print(sum(log(LLi) * Z))
    print(round(cbind(Pr, Pr0), 3))
    print(round(cbind(oldbeta, beta0), 3))
    print(round(c(oldsigma, sigma0), 3))
    print(round(cbind(oldtau, tau0), 3))
  }
}

####�T���v�����O���ʂ̗v��Ɖ���####
burnin <- 2000/keep
RS <- R/keep

##�T���v�����O���ʂ̉���
matplot(BETA[, 1:3], type="l", xlab="�T���v�����O��", ylab="�p�����[�^����l")
matplot(BETA[, 4:6], type="l", xlab="�T���v�����O��", ylab="�p�����[�^����l")
matplot(BETA[, 7:9], type="l", xlab="�T���v�����O��", ylab="�p�����[�^����l")
matplot(BETA[, 10:12], type="l", xlab="�T���v�����O��", ylab="�p�����[�^����l")
plot(1:RS, SIGMA, type="l", xlab="�T���v�����O��", ylab="�p�����[�^����l")
matplot(LAMBDA[, 1:s], type="l", xlab="�T���v�����O��", ylab="�p�����[�^����l")
matplot(LAMBDA[, (s+1):ncol(LAMBDA)], type="l", xlab="�T���v�����O��", ylab="�p�����[�^����l")
matplot(t(PR[1, , ]), type="l", xlab="�T���v�����O��", ylab="���ڊm���̃T���v�����O����")
matplot(t(PR[2, , ]), type="l", xlab="�T���v�����O��", ylab="���ڊm���̃T���v�����O����")
matplot(t(PR[3, , ]), type="l", xlab="�T���v�����O��", ylab="���ڊm���̃T���v�����O����")
matplot(t(PR[4, , ]), type="l", xlab="�T���v�����O��", ylab="���ڊm���̃T���v�����O����")

##�T���v�����O���ʂ̗v��
#���`��A�̉�A�W���̐���l
round(cbind(matrix(colMeans(BETA[burnin:RS, ]), nrow=k1+1, ncol=s), beta0), 3)
round(matrix(apply(BETA[burnin:RS, ], 2, sd), nrow=k1+1, ncol=s), 3)

#���`��A�̕W���΍��̐���l
round(c(mean(SIGMA[burnin:RS]), sigma0), 3)
round(sd(SIGMA[burnin:RS]), 3)

#�x�[�^���z�̃p�����[�^�̐���l
round(cbind(matrix(colMeans(LAMBDA[burnin:RS, ]), nrow=s, ncol=2), tau0), 3)
round(matrix(apply(LAMBDA[burnin:RS, ], 2, sd), nrow=s, ncol=2), 3)

#�}���R�t���ڊm���̐���l
round(cbind(apply(PR[, , burnin:RS], c(1, 2), mean), Pr0), 3)
round(apply(PR[, , burnin:RS], c(1, 2), sd), 3)

##�\�����x���m�F
W <- t(apply(rbind(Zi, matrix(1:s, nrow=s, ncol=n)), 2, table)) - 1
W_rate <- W / rowSums(W)   #���ݕϐ��̊����m��
w <- apply(W_rate, 1, which.max)   #���ݕϐ��̊���
round(cbind(W_rate, w, z0), 3)   #�^�̐��ݕϐ��Ƃ̔�r
sum(w==z0) / n   #�\���̐�����

