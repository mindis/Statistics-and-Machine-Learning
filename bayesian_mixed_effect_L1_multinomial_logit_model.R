#####�ϗʌ��ʐ������������W�b�g���f��#####
library(MASS)
library(mclust)
library(reshape2)
library(bayesm)
library(ExtDist)
library(extraDistr)
library(matrixStats)
library(glmnet)
library(monomvn)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)
gc(); gc()

#set.seed(534798)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 2000
pt <- ceiling(rgamma(hh, 3.7, 0.9))
hhpt <- sum(pt)
select <- 10

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
k <- 40   #�����ʐ�
data1 <- extraDistr::rdirichlet(hhpt, rep(0.2, k/2))[, -k/2]
data2 <- extraDistr::rdirichlet(hhpt, rep(0.2, k/2))[, -k/2]
Data <- cbind(1, data1, data2)

####�����ϐ��̔���#### 
##�ؕЂ��x�N�g���ϊ�
X <- matrix(diag(select), nrow=hhpt*select, ncol=select, byrow=T)[, -select]

for(i in 1:1000){
  print(i) 
  
  ##�p�����[�^�̐ݒ�
  #lasso���f���̉�A�p�����[�^��ݒ�
  b01 <- runif(select-1, -0.8, 0.8)
  b02 <- matrix(runif((k-2)*(select-1), -3.25, 4.0), nrow=k-2, ncol=select-1) * 
              matrix(rbinom((k-2)*(select-1), 1, 0.3), nrow=k-2, ncol=select-1)
  b0 <- rbind(b01, b02)
  
  #���U��ݒ�
  sigma0 <- 0.3
  
  ##�ϗʌ��ʂ̃p�����[�^��ݒ�
  cov0 <- diag(0.25, select-1)
  theta0 <- mvrnorm(hh, rep(0, select-1), diag(0.25, select-1))  
  
  ##�����ϐ��̔���
  logit <- matrix(0, nrow=hhpt, ncol=select)
  tau0 <- matrix(0, nrow=hhpt, ncol=select-1)
  
  #���W�b�g�̒�`
  for(j in 1:(select-1)){
    tau0[, j] <- rnorm(hhpt, 0, sigma0)
    logit[, j] <- (Data %*% b0[, j] + tau0[, j]) + theta0[ID$id, j]
  }
  
  #�������z���牞���ϐ��𔭐�
  Pr0 <- exp(logit) / rowSums(exp(logit))
  y <- t(apply(Pr0, 1, function(x) rmultinom(1, 1, x)))
  
  if(max(colMeans(y)) < 0.45 & min(colMeans(y)) > 0.05) break
}

#�����������f�[�^���m�F
colSums(y); round(colMeans(Pr0), 3)
rownames(b0) <- NULL
sum(b0!=0)

####�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y�������������W�b�g���f���𐄒�####
##�������W�b�g���f���̑ΐ��ޓx�֐���ݒ�
fr <- function(y, X, beta, gamma, select, hhpt){
  
  #���W�b�g�Ɗm���̌v�Z
  logit <- cbind(beta + gamma, 0)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hhpt, ncol=select)
  
  #�ΐ��ޓx���`
  LLi <- rowSums(y * log(Pr))
  return(LLi)
}


##�A���S���Y���̐ݒ�
R <- 40000
keep <- 4
sbeta <- 1.5
iter <- 0

##���O���z�̐ݒ�
nu <- select
V <- solve(nu * diag(select-1))


##�����l�̐ݒ�
#�x�C�W�A��lasso��A�̏����l
beta0 <- scale(colSums(y))
oldgamma <- mvrnorm(hhpt, beta0[-select], diag(0.2, select-1))
oldtheta <- matrix(0, nrow=k-2+1, ncol=select-1)
oldsigma <- diag(0.1, select-1)
inv_sigma <- solve(oldsigma)

#�ϗʌ��ʂ̏����l
oldcov <- diag(0.1, select-1)
inv_cov <- solve(oldcov)
oldbeta <- mvrnorm(hh, rep(0, select-1), oldcov)

##�T���v�����O���ʂ̕ۑ��p�z��
THETA <- array(0, dim=c(k-2+1, select-1, R/keep))
BETA <- array(0, dim=c(hh, select-1, R/keep))
COV <- array(0, dim=c(select-1, select-1, R/keep))
SIGMA <- matrix(0, nrow=R/keep, ncol=select-1)

##�C���f�b�N�X�̐ݒ�
index_id <- list()
index_y <- list()
for(i in 1:hh){
  index_id[[i]] <- which(ID_v$id==i)
  index_y[[i]] <- which(ID$id==i)
}

lognew1 <- rep(0, hh)
logold1 <- rep(0, hh)
logpnew1 <- rep(0, hh)
logpold1 <- rep(0, hh)
lambda <- rep(1, select-1)
er_new <- matrix(0, nrow=hh, select-1)
er_old <- matrix(0, nrow=hh, select-1)

#ID�̐ݒ�
id_vec <- c()
for(i in 1:hh){
  id_vec <- c(id_vec, rep(i, sum(ID$id==i)*select))
}


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){

  ##�ϗʌ��ʂ��T���v�����O
  #�V�����p�����[�^���T���v�����O
  betad <- oldbeta 
  betan <- betad + mvrnorm(hh, rep(0, select-1), diag(0.025, select-1))
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew0 <- fr(y, X, betan[ID$id, ], oldgamma, select, hhpt)
  logold0 <- fr(y, X, betad[ID$id, ], oldgamma, select, hhpt)
  logpnew1 <- -0.5 * rowSums(betan %*% inv_cov * betan)
  logpold1 <- -0.5 * rowSums(betad %*% inv_cov * betad)
  
  for(i in 1:hh){
    lognew1[i] <- sum(lognew0[index_y[[i]]])
    logold1[i] <- sum(logold0[index_y[[i]]])
  }
  
  #���g���|���X�w�C�X�e�B���O�@�Ńp�����[�^�̍̑�������
  rand <- runif(hh)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew1 + logpnew1 - logold1 - logpold1)   #�̑𗦂��v�Z
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=hh, ncol=select-1)
  oldbeta <- flag*betan + (1-flag)*betad   #alpha��rand�������Ă�����̑�
  
  ##�t�E�B�V���[�g���z����K�w���f���̕��U�����U�s����T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^
  R_par <- V + t(oldbeta) %*% oldbeta
  Sn <- nu + hh
  
  #�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
  oldcov <- rwishart(Sn, solve(R_par))$IW
  inv_cov <- solve(oldcov)
  
  
  ##�������W�b�g���f���̃T���v�����Ƃ̃p�����[�^���T���v�����O
  #�V�����p�����[�^���T���v�����O
  gammad <- oldgamma
  gamman <- gammad + mvrnorm(hhpt, rep(0, select-1), diag(0.05 ,select-1))
  
  #�덷��ݒ�
  mu <- Data %*% oldtheta
  er_new <- gamman - mu
  er_old <- gammad - mu
  
  #lasso�̕��U�����U�s��𐄒�
  oldsigma <- diag(diag(var(er_old)))
  inv_sigma <- solve(oldsigma)
  
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew2 <- fr(y, X, oldbeta[ID$id, ], gamman, select, hhpt)
  logold2 <- fr(y, X, oldbeta[ID$id, ], gammad, select, hhpt)
  logpnew2 <- -0.5 * rowSums(er_new %*% inv_sigma * er_new)
  logpold2 <- -0.5 * rowSums(er_old %*% inv_sigma * er_old)
  
  #���g���|���X�w�C�X�e�B���O�@�Ńp�����[�^�̍̑�������
  rand <- runif(hhpt)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew2 + logpnew2 - logold2 - logpold2)   #�̑𗦂��v�Z
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff

  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=hhpt, ncol=select-1)
  oldgamma <- flag*gamman + (1-flag)*gammad   #alpha��rand�������Ă�����̑�
  
  
  ##�x�C�W�A��lasso�ŊK�w���f���̉�A�p�����[�^���T���v�����O
  for(j in 1:(select-1)){
    res <- blasso(X=Data[, -1], y=oldgamma[, j], beta=oldtheta[-1, j], lambda2=lambda[j], s2=diag(oldsigma)[j], T=2)
    oldtheta[, j] <- c(res$mu[2], res$beta[2, ])
    lambda[j] <- res$lambda2[2]
  }

  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- oldtheta
    BETA[, , mkeep] <- oldbeta
    COV[, , mkeep] <- oldcov
    SIGMA[mkeep, ] <- diag(oldsigma)
    print(rp)
    print(sum(lognew1))
    print(round(lambda, 3))
    print(round(t(cbind(oldtheta, b0)[1:15, ]), 3))
    print(round(rbind(diag(oldcov), diag(cov0)), 3))
  }
}

matplot(t(THETA[4, , ]), type="l")
matplot(t(BETA[10, , ]), type="l")


