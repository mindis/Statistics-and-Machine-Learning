#####�K�w�x�C�Y�������������W�b�g���f��#####
library(MASS)
library(mclust)
library(reshape2)
library(bayesm)
detach("package:bayesm", unload=TRUE)
library(ExtDist)
library(extraDistr)
library(matrixStats)
library(glmnet)
library(monomvn)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(534798)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 2000
pt0 <- rpois(hh, 3.8)
pt <- ifelse(pt0==0, 1, pt0)
hhpt <- sum(pt)
select <- 8

#ID�̐ݒ�
id <- rep(1:hh, pt)
time <- c()
for(i in 1:hh){
  time <- c(time, 1:pt[i])
}
ID <- data.frame(no=1:hhpt, id, time)

#ID�̃x�N�g����
id_v <- as.numeric(t(matrix(id, nrow=hhpt, ncol=select)))
time_v <- as.numeric(t(matrix(time, nrow=hhpt, ncol=select)))
ID_v <- data.frame(no=1:length(id_v), id=id_v, time=time_v)


####�����ϐ��̔���####
k <- 100   #�����ʐ�
cont <- 70   #�A���ϐ���
bin <- 30   #��l�ϐ���
Data <- array(0, dim=c(hh, k+1, select-1))

for(i in 1:(select-1)){
  p <- runif(1, 0.4, 0.6)
  x.cont <- matrix(rnorm(hh*cont, 0, 1), nrow=hh, ncol=cont)
  x.bin <- matrix(rbinom(hh*bin, 1, p), nrow=hh, ncol=bin)
  Data[, , i] <- cbind(1, x.cont, x.bin)
}

####�����ϐ��̔���#### 
##�ؕЂ��x�N�g���ϊ�
X <- matrix(diag(select), nrow=hhpt*select, ncol=select, byrow=T)[, -select]

for(i in 1:1000){
  print(i) 
  
  ##�p�����[�^�̐ݒ�
  #�K�w���f���̉�A�p�����[�^��ݒ�
  b01 <- runif(select-1, -0.8, 0.8)
  b02 <- matrix(runif(cont*(select-1), 0, 1.1), nrow=cont, ncol=select-1) * 
    matrix(rbinom(cont*(select-1), 1, 0.3), nrow=cont, ncol=select-1)
  b03 <- matrix(runif(bin*(select-1), -1.1, 1.3), nrow=bin, ncol=select-1) * 
    matrix(rbinom(bin*(select-1), 1, 0.3), nrow=bin, ncol=select-1)
  b0 <- rbind(b01, b02, b03)
  
  #�K�w���f���̕��U��ݒ�
  cov0 <- 0.35
  
  ##�����ϐ��̔���
  logit <- matrix(0, nrow=hhpt, ncol=select)
  tau0 <- matrix(0, nrow=hh, ncol=select-1)
  
  #���W�b�g�̒�`
  for(j in 1:(select-1)){
    tau0[, j] <- rnorm(hh, 0, cov0)
    logit[, j] <- (Data[, , j] %*% b0[, j] + tau0[, j])[ID$id, ]
  }
  
  #�������z���牞���ϐ��𔭐�
  Pr0 <- exp(logit) / rowSums(exp(logit))
  y <- t(apply(Pr0, 1, function(x) rmultinom(1, 1, x)))
  
  if(max(colMeans(y)) < 0.4 & min(colMeans(y)) > 0.05) break
}

#�����������f�[�^���m�F
colSums(y); round(colMeans(Pr0), 3)
rownames(b0) <- NULL

####�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y�������������W�b�g���f���𐄒�####
##�������W�b�g���f���̑ΐ��ޓx�֐���ݒ�
fr <- function(y, X, beta, select, n){
  
  #���W�b�g�Ɗm���̌v�Z
  logit <- matrix(X %*% beta, nrow=n, ncol=select, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=n, ncol=select)
  
  #�ΐ��ޓx���`
  LLi <- rowSums(y * log(Pr))
  LL <- sum(LLi)
  return(LL)
}

##�A���S���Y���̐ݒ�
R <- 20000
keep <- 4
sbeta <- 1.5
iter <- 0

##�����l�̐ݒ�
beta0 <- scale(colSums(y))
oldbeta <- mvrnorm(hh, beta0[-select], diag(0.2, select-1))
oldtheta <- matrix(0, nrow=k+1, ncol=select-1)
oldcov <- diag(0.1, select-1)
inv_cov <- solve(oldcov)

##�T���v�����O���ʂ̕ۑ��p�z��
THETA <- array(0, dim=c(k+1, select-1, R/keep))
BETA <- array(0, dim=c(hh, select-1, R/keep))
TAU <- array(0, dim=c(select-1, select-1, R/keep))

##�C���f�b�N�X�̐ݒ�
index_id <- list()
index_y <- list()
for(i in 1:hh){
  index_id[[i]] <- which(ID_v$id==i)
  index_y[[i]] <- which(ID$id==i)
}
lognew <- rep(0, hh)
logold <- rep(0, hh)
logpnew <- rep(0, hh)
logpold <- rep(0, hh)
lambda <- rep(1, select-1)
er_new <- matrix(0, nrow=hh, select-1)
er_old <- matrix(0, nrow=hh, select-1)


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�������W�b�g���f����ID�ʃp�����[�^���T���v�����O
  #�V�����p�����[�^���T���v�����O
  betad <- oldbeta
  betan <- betad + mvrnorm(hh, rep(0, select-1), diag(0.025 ,select-1))
  
  #�덷��ݒ�
  for(j in 1:(select-1)){
    mu <- Data[, , j] %*% oldtheta[, j]
    er_new[, j] <- betan[, j] - mu
    er_old[, j] <- betad[, j] - mu
  }
  
  #�K�w���f���̕��U�����U�s��𐄒�
  oldcov <- var(er_old)
  inv_cov <- solve(oldcov)
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  for(i in 1:hh){
    lognew[i] <- fr(y[index_y[[i]], ], X[index_id[[i]], ], betan[i, ], select, length(index_y[[i]]))
    logold[i] <- fr(y[index_y[[i]], ], X[index_id[[i]], ], betad[i, ], select, length(index_y[[i]]))
    logpnew[i] <- -0.5 * er_new[i, ] %*% inv_cov %*% er_new[i, ]
    logpold[i] <- -0.5 * er_old[i, ] %*% inv_cov %*% er_old[i, ]
  }
  
  #MH�T���v�����O�Ńp�����[�^�̍̑�������
  #���g���|���X�w�C�X�e�B���O�@�Ńp�����[�^�̍̑�������
  rand <- runif(hh)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew + logpnew - logold - logpold)   #�̑𗦂��v�Z
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=hh, ncol=select-1)
  oldbeta <- flag*betan + (1-flag)*betad   #alpha��rand�������Ă�����̑�
  
  
  ##�x�C�W�A��lasso�ŊK�w���f���̉�A�p�����[�^���T���v�����O
  for(j in 1:(select-1)){
    res <- blasso(X=Data[, -1, j], y=oldbeta[, j], beta=oldtheta[-1, j], lambda2=lambda[j], s2=diag(oldcov)[j], T=2)
    oldtheta[, j] <- c(res$mu[2], res$beta[2, ])
    lambda[j] <- res$lambda2[2]
  }
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- oldtheta
    BETA[, , mkeep] <- oldbeta
    TAU[, , mkeep] <- oldcov
    print(rp)
    print(sum(lognew))
    print(round(t(cbind(oldtheta, b0)[1:20, ]), 3))
  }
}

matplot(t(THETA[7, , ]), type="l")


