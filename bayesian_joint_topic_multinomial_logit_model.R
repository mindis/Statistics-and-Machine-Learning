#####�Ή��g�s�b�N���f��#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
detach("package:gtools", unload=TRUE)
detach("package:bayesm", unload=TRUE)
library(extraDistr)
library(monomvn)
library(glmnet)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

set.seed(8079)

####�f�[�^�̔���####
#set.seed(423943)
#�����f�[�^�̐ݒ�
k <- 10   #�g�s�b�N��
d <- 7000   #������
v <- 300   #��b��
w <-  rpois(d, rgamma(d, 30, 0.5))   #1����������̒P�ꐔ
a <- 100   #�⏕�ϐ���
x <- rpois(d, 20)
select <- 7   #�I������

#�p�����[�^�̐ݒ�
alpha0 <- rep(0.2, k) #�����̃f�B���N�����O���z�̃p�����[�^
alpha1 <- rep(0.15, v)   #�P��̃f�B���N�����O���z�̃p�����[�^
alpha2 <- rep(0.25, a)   #�⏕���̃f�B�N�������O���z�̃p�����[�^

#�f�B���N�������̔���
theta0 <- theta <- extraDistr::rdirichlet(d, alpha0)   #�����̃g�s�b�N���z���f�B���N���������甭��
phi0 <- phi <- extraDistr::rdirichlet(k, alpha1)   #�P��̃g�s�b�N���z���f�B���N���������甭��
omega0 <- omega <- extraDistr::rdirichlet(k, alpha2) #�⏕�f�[�^�̃g�s�b�N���z���f�B�N�����������甭��

#�������z�̗�������f�[�^�𔭐�
WX <- matrix(0, nrow=d, ncol=v)
AX <- matrix(0, nrow=d, ncol=a)
Z1 <- list()
Z2 <- list()

for(i in 1:d){
  print(i)
  
  #�����̃g�s�b�N���z�𔭐�
  z1 <- t(rmultinom(w[i], 1, theta[i, ]))   #�����̃g�s�b�N���z�𔭐�
  
  #�����̃g�s�b�N���z����P��𔭐�
  zn <- z1 %*% c(1:k)   #0,1�𐔒l�ɒu��������
  zdn <- cbind(zn, z1)   #apply�֐��Ŏg����悤�ɍs��ɂ��Ă���
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi[x[1], ])))   #�����̃g�s�b�N����P��𐶐�
  wdn <- colSums(wn)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  WX[i, ] <- wdn  
  
  #�����̃g�s�b�N���z����⏕�ϐ��𔭐�
  z2 <- t(rmultinom(x[i], 1, theta[i, ]))
  zx <- z2 %*% 1:k
  zax <- cbind(zx, z2)
  an <- t(apply(zax, 1, function(x) rmultinom(1, 1, omega[x[1], ])))
  adn <- colSums(an)
  AX[i, ] <- adn
  
  #�����g�s�b�N����ѕ⏕���g�s�b�N���i�[
  Z1[[i]] <- z1
  Z2[[i]] <- z2
}

#�f�[�^�s��𐮐��^�s��ɕύX
storage.mode(WX) <- "integer"
storage.mode(AX) <- "integer"

####�����ϐ��̔���####
#�����ϐ��̊i�[�p�z��
y <- matrix(0, nrow=d, ncol=select)
Pr <- matrix(0, nrow=d, ncol=select)
Pr0 <- matrix(0, nrow=d, ncol=select)

##�Ó��ȉ����ϐ�����������܂Ŕ���������
for(j in 1:5000){
  ##�p�����[�^�̐ݒ�
  #�g�s�b�N���f���̃p�����[�^
  sparse1 <- matrix(rbinom((select-1)*k, 1, 0.4), nrow=k, ncol=select-1)   #�p�����[�^�̃X�p�[�X�s��
  b00 <- runif(select-1, -0.5, 0.5)
  b01 <- (matrix(runif((select-1)*k, -4.0, 4.0), nrow=k, ncol=select-1)) * sparse1
  b02 <- (b01 + mvrnorm(k, rep(0, select-1), diag(0.2, select-1))) * sparse1
  b0 <- rbind(b00, b01, b02)
  rownames(b0) <- NULL
  
  ##�������ƂɊm���Ɖ����ϐ��𔭐�
  for(i in 1:d){
    #���݃g�s�b�Nz���W�v
    z1 <- colSums(Z1[[i]])
    z1 <- ifelse(z1==0, 1, z1)
    z2 <- colSums(Z2[[i]])
    z2 <- ifelse(z2==0, 1, z2)
    
    #�������z���牞���ϐ��𔭐�
    logit <- c(c(1, log(z1), log(z2)) %*% b0, 0)
    Pr[i, ] <- exp(logit) / sum(exp(logit))
    y[i, ] <- t(rmultinom(1, 1, Pr[i, ]))
  }
  
  t1 <- sum(apply(Pr, 1, which.max)==y %*% 1:select)/d
  print(round(c(t1, min(colSums(y))), 3))
  
  if(t1 > 0.85 &  min(colSums(y)) > 250) break
}



####�}���R�t�A�������e�J�����@�Ō����g�s�b�N���f��+L1���������W�b�g���f���𐄒�####
####�g�s�b�N���f���̂��߂̃f�[�^�Ɗ֐��̏���####
##���ꂼ��̕������̒P��̏o������ѕ⏕���̏o�����x�N�g���ɕ��ׂ�
##�f�[�^����pID���쐬
ID1_list <- list()
wd_list <- list()
ID2_list <- list()
ad_list <- list()

#���l���Ƃɋ��lID����ђP��ID���쐬
for(i in 1:nrow(WX)){
  print(i)
  
  #�P���ID�x�N�g�����쐬
  ID1_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
  num2 <- subset(num1, num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
  
  #�⏕����ID�x�N�g�����쐬
  ID2_list[[i]] <- rep(i, x[i])
  num1 <- (AX[i, ] > 0) * (1:a)
  num2 <- subset(num1, num1 > 0)
  A1 <- AX[i, (AX[i, ] > 0)]
  number <- rep(num2, A1)
  ad_list[[i]] <- number
}

#���X�g���x�N�g���ɕϊ�
ID1_d <- unlist(ID1_list)
ID2_d <- unlist(ID2_list)
wd <- unlist(wd_list)
ad <- unlist(ad_list)

##�C���f�b�N�X���쐬
doc1_list <- list()
word_list <- list()
doc2_list <- list()
aux_list <- list()

for(i in 1:length(unique(ID1_d))) {doc1_list[[i]] <- which(ID1_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- which(wd==i)}
for(i in 1:length(unique(ID2_d))) {doc2_list[[i]] <- which(ID2_d==i)}
for(i in 1:length(unique(ad))) {aux_list[[i]] <- which(ad==i)}
gc(); gc()


####�}���R�t�A�������e�J�����@�őΉ��g�s�b�N���f���𐄒�####
##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #���S�W���̊i�[�p
  for(kk in 1:k){
    #���S�W�����v�Z
    Bi <- rep(theta[, kk], w) * phi[kk, c(wd)]   #�ޓx
    Bur[, kk] <- Bi   
  }
  
  Br <- Bur / rowSums(Bur)   #���S���̌v�Z
  r <- colSums(Br) / sum(Br)   #�������̌v�Z
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##�������W�b�g���f���̑ΐ��ޓx
fr <- function(beta, y, x, hh, select){
  
  #���W�b�g�Ɗm���̌v�Z
  logit <- t(x %*% t(beta))
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #�ΐ��ޓx��ݒ�
  LLi <- rowSums(y*log(Pr)) 
  return(LLi)
}


##�A���S���Y���̐ݒ�
R <- 20000   #�T���v�����O��
keep <- 4   #4���1��̊����ŃT���v�����O���ʂ��i�[
iter <- 0

##�w�K�f�[�^�ƃe�X�g�f�[�^�ɕ���
index_test <- sample(1:d, 1000)
n1 <- d-length(index_test)
n2 <- length(index_test)
y_train <- y[-index_test, ]
y_test <- y[index_test, ]
y_vec <- y_test %*% 1:select


##���O���z�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z
alpha01 <- rep(1.0, k)
beta0 <- rep(0.5, v)
gamma0 <- rep(0.5, a)
alpha01m <- matrix(alpha01, nrow=d, ncol=k, byrow=T)
beta0m <- matrix(beta0, nrow=v, ncol=k)
gamma0m <- matrix(gamma0, nrow=a, ncol=k)


##�p�����[�^�̏����l
#�g�s�b�N���f���̏����l
theta.ini <- runif(k, 0.5, 2)
phi.ini <- runif(v, 0.5, 1)
omega.ini <- runif(a, 0.5, 1)
theta <- rdirichlet(d, theta.ini)   #�����g�s�b�N�̃p�����[�^�̏����l
phi <- rdirichlet(k, phi.ini)   #�P��g�s�b�N�̃p�����[�^�̏����l
omega <- rdirichlet(k, omega.ini)   #�⏕���g�s�b�N�̃p�����[�^�̏����l

#�������W�b�g���f���̏����l
par <- 2*k+1
beta0 <- scale(colSums(y_train))
oldbeta <- mvrnorm(n1, beta0[-select], diag(0.2, select-1))
oldtheta <- matrix(0, nrow=par, ncol=select-1)
oldcov <- diag(0.1, select-1)
inv_cov <- solve(oldcov)
lambda <- c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01)


##�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
OMEGA <- array(0, dim=c(k, a, R/keep))
W_SEG <- matrix(0, nrow=R/(keep*10), ncol=sum(w))
A_SEG <- matrix(0, nrow=R/(keep*10), ncol=sum(x))
storage.mode(W_SEG) <- "integer"
storage.mode(A_SEG) <- "integer"
gc(); gc()

##MCMC����p�z��
wsum0 <- matrix(0, nrow=d, ncol=k)
vf0 <- matrix(0, nrow=v, ncol=k)
af0 <- matrix(0, nrow=a, ncol=k)
asum0 <- matrix(0, nrow=d, ncol=k)
x_diag <- diag(select)[, -select]


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�P��g�s�b�N���T���v�����O
  #�P�ꂲ�ƂɃg�s�b�N�̏o�������v�Z
  word_rate <- burden_fr(theta, phi, wd, w, k)$Br
  
  #�������z����P��g�s�b�N���T���v�����O
  Zi1 <- rmnom(nrow(word_rate), 1, word_rate) 
  
  ##�����g�s�b�N�̃p�����[�^���X�V
  #�f�B�N�������z����theta���T���v�����O
  for(i in 1:d){
    wsum0[i, ] <- colSums(Zi1[doc1_list[[i]], ]) 
  }
  wsum <- wsum0 + alpha01m   #�f�B�N�������z�̃p�����[�^
  theta <- extraDistr::rdirichlet(d, wsum)   #�f�B�N�������z����theta���T���v�����O
  
  ##�P��g�s�b�N�̃p�����[�^���X�V
  #�f�B�N�������z����phi���T���v�����O
  for(i in 1:v){
    vf0[i, ] <- colSums(Zi1[word_list[[i]], ])
  }
  vf <- t(vf0 + beta0m)   #�f�B�N�������z�̃p�����[�^
  phi <- extraDistr::rdirichlet(k, vf)   #�f�B�N�������z����phi���T���v�����O

  
  ##�⏕���g�s�b�N���T���v�����O
  #�����g�s�b�N��p���ĕ⏕��񂲂ƂɃg�s�b�N�̏o�������v�Z
  aux_rate <- burden_fr(theta, omega, ad, x, k)$Br
  
  #�������z����⏕���g�s�b�N���T���v�����O
  Zi2 <- rmnom(nrow(aux_rate), 1, aux_rate) 
  
  
  ##�⏕���g�s�b�N�̃p�����[�^���X�V
  #�f�B�N�������z����omega���T���v�����O
  for(i in 1:a){
    af0[i, ] <- colSums(Zi2[aux_list[[i]], ])
  }
  af <- t(af0 + gamma0m)   #�f�B�N�������z�̃p�����[�^
  omega <- extraDistr::rdirichlet(k, af)   #�f�B�N�������z����omega���T���v�����O
  
  
  ####�x�C�W�A��L1���������W�b�g���f����beta���T���v�����O####
  ##�f�[�^�̐ݒ�
  #�����g�s�b�N�̐����ϐ��̐ݒ�
  wsum0[wsum0==0] <- 1
  Data01 <- log(wsum0)
  
  #�⏕���g�s�b�N�̐����ϐ��̐ݒ�
  for(i in 1:d){
    asum0[i, ] <- colSums(Zi2[doc2_list[[i]], ])
  }
  asum0[asum0==0] <- 1
  Data02 <- log(asum0)
  
  #�w�K�f�[�^�ƃe�X�g�f�[�^�ɕ���
  Data <- cbind(1, Data01, Data02)
  Data1 <- Data[-index_test, ]
  Data2 <- Data[index_test, ]

  
  ##�������W�b�g���f���̃p�����[�^���T���v�����O
  #�V�����p�����[�^���T���v�����O
  betad <- oldbeta
  betan <- betad + mvrnorm(n1, rep(0, select-1), diag(0.25, select-1))
  
  #�덷��ݒ�
  er_new <- betan - mu
  er_old <- betad - mu

  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew <- fr(betan, y_train, x_diag, n1, select)
  logold <- fr(betad, y_train, x_diag, n1, select)
  logpnew <- -0.5 * rowSums(er_new %*% inv_cov * er_new)
  logpold <- -0.5 * rowSums(er_old %*% inv_cov * er_old)  
  
  #���g���|���X�w�C�X�e�B���O�@�Ńp�����[�^�̍̑�������
  rand <- runif(n1)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew + logpnew - logold - logpold)   #�̑𗦂��v�Z
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=n1, ncol=select-1)
  oldbeta <- flag*betan + (1-flag)*betad   #alpha��rand�������Ă�����̑�
  
  ##lasso�ŊK�w���f���̉�A�p�����[�^���T���v�����O
  for(j in 1:(select-1)){
    
    #lambda�����ȏ�Ȃ狭���I�ɏ�����lambda�ɂ���
    if(lambda[j] > 0.5){
      lambda[j] <- 0.001
    }
    
    #�x�C�W�A��lasso�ŊK�w���f���̉�A�W�����X�V
    res <- blasso(X=Data1[, -1], y=oldbeta[, j], beta=oldtheta[-1, j], lambda2=lambda[j], s2=diag(oldcov)[j], 
                  normalize=TRUE, T=2)
    oldtheta[, j] <- c(res$mu[2], res$beta[2, ])
    lambda[j] <- res$lambda2[2]
  }
  
  ##�K�w���f���̕��U�����U�s��𐄒�
  mu <- Data1 %*% oldtheta
  oldcov <- var(oldbeta - mu)
  inv_cov <- solve(oldcov)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    mkeep1 <- rp/keep
    THETA[, , mkeep1] <- theta
    PHI[, , mkeep1] <- phi
    OMEGA[, , mkeep1] <- omega
    #if(rp%%(keep*10)==0){
    #  mkeep2 <- rp/(keep*5)
    #  W_SEG[mkeep2, ] <- word_z
    #  A_SEG[mkeep2, ] <- aux_z
    #}
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(round(c(sum(lognew), mean(alpha)), 3))
    print(round(lambda, 4))
    print(round(cbind(theta[1:10, ], theta0[1:10, ]), 3))
    #print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
    #print(round(cbind(omega[, 1:10], omegat[, 1:10]), 3))
    
    #�\�����z���v�Z
    logit <- Data2 %*% cbind(oldtheta, 0)
    Pr <- exp(logit) / rowSums(exp(logit))
    print(mean(y[index_test, ] %*% 1:select==apply(Pr, 1, which.max)))
  }
}


round(cbind(oldtheta, b0), 3)
round(Pr, 3)
round(cbind(Pr, y[index_test, ]), 3)