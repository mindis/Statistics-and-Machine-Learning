#####�g�s�b�N���f�����܂񂾃��W�X�e�B�b�N��A���f��#####
library(MASS)
library(lda)
library(bayesm)
library(RMeCab)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(317209)

####�f�[�^�̔���####
#�f�[�^�̐ݒ�
k <- 5   #�g�s�b�N��
d <- 1000   #���[�U�[��
v <- 50   #�A�C�e����
w <- rpois(d, 50)   #�w����

####�g�s�b�N���f������A�C�e���w���𐶐�####
##�g�s�b�N���f���̐ݒ�
#�p�����[�^�̐ݒ�
alpha0 <- runif(k, 0.1, 0.8)   #���[�U�[�̃f�B���N�����O���z�̃p�����[�^
alpha1 <- rep(0.25, v)   #�A�C�e���̃f�B���N�����O���z�̃p�����[�^

#�f�B���N�������̔���
theta0 <- rdirichlet(d, alpha0)   #���[�U�[�̃g�s�b�N���z���f�B���N���������甭��
phi0 <- rdirichlet(k, alpha1)   #�A�C�e���̃g�s�b�N���z���f�B���N���������甭��

#�������z����f�[�^�𔭐�
WX <- matrix(0, nrow=d, ncol=v)
Z <- matrix(0, nrow=d, ncol=k)
Rate <- matrix(0, nrow=d, ncol=k)
ZS <- list()

for(i in 1:d){
  #�g�s�b�N�𐶐�
  z <- t(rmultinom(w[i], 1, theta0[i, ]))   #���[�U�[�̃A�C�e���̃g�s�b�N�𐶐�
  Z[i, ] <- colSums(z)
  Rate[i, ] <- Z[i, ] / w[i]
  zn <- z %*% c(1:k)   #0,1�𐔒l�ɒu��������
  zdn <- cbind(zn, z)   #apply�֐��Ŏg����悤�ɍs��ɂ��Ă���
 
  #�g�s�b�N���牞���ϐ��𐶐�
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi0[x[1], ])))   #���[�U�[�̃g�s�b�N����A�C�e���𐶐�
  
  wdn <- colSums(wn)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  WX[i, ] <- wdn
  ZS[[i]] <- cbind(rep(i, w[i]), zdn[, 1])
  print(i)
}

#���X�g���s������ɕύX
ZS <- do.call(rbind, ZS)

#�g�s�b�N�̒P���W�v
z_table <- table(ZS[, 2])
z_ratet <- z_table/sum(z_table)


####�����������g�s�b�N���牞���ϐ��𔭐�####
##�⏕�ϐ��𔭐�
#�A���ϐ��̔���
k1 <- 3
X1 <- matrix(rnorm(d*k1, 0, 1), nrow=d, ncol=k1)

#��l�ϐ��̔���
k2 <- 3
X2 <- matrix(0, d, k2)
for(i in 1:k2){
  X2[, i] <- rbinom(d, 1, runif(1, 0.3, 0.7))
}

#�g�s�b�N�̕s�v�ϐ����폜
ColSums(Rate1)
Rate1 <- Rate[, -which.min(colSums(Rate))]

#�f�[�^������
X <- cbind(1, X1, X2)
XZ <- cbind(X, Rate1)


##���W�X�e�B�b�N��A���f�����牞���ϐ��𔭐�
#�p�����[�^�̐ݒ�
b0 <- runif(1, -1.2, -0.5)
b1 <- c(runif(k1, 0, 0.9), runif(k2, -0.9, 1.1))
b2 <- runif(k-1, -2.5, 4.0)
b <- c(b0, b1, b2)

#�m���̌v�Z
logit <- XZ %*% b   #���W�b�g
Pr <- exp(logit)/(1+exp(logit))

#�񍀕��z���牞���ϐ��𔭐�
y <- rbinom(d, 1, Pr)

loglt <- as.numeric(logLik(glm(y ~ ., data.frame(XZ[, -1]), family=binomial)))

res1 <- glm(y ~ ., data.frame(XZ[, -1]), family=binomial)
res2 <- glm(y ~ ., data.frame(XZ[, -c(1, 8:11)]), family=binomial)
summary(res1)
summary(res2)


####�}���R�t�A�������e�J�����@�Ńg�s�b�N���f��+���W�X�e�B�b�N��A���f���𐄒�####
##�A���S���Y���̐ݒ�
R <- 10000   #�T���v�����O�� 
keep <- 2
iter <- 0

##���O���z�̃p�����[�^�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z�̃p�����[�^
alpha01 <- alpha0   #�����̃f�B���N�����O���z�̃p�����[�^
beta01 <- alpha1[1]   #�P��̃f�B���N�����O���z�̃p�����[�^

#���W�X�e�b�N��A���f���̎��O���z�̃p�����[�^
betas <- rep(0, ncol(XZ))  #��A�W���̏����l
betaBi <- rep(0, ncol(XZ))   #��A�W���̎��O���z�̕���
rootBi <- 0.01*diag(ncol(XZ))   #���O���z�̐��x

##�p�����[�^�̏����l�̐ݒ�
theta.ini <- runif(k, 0.3, 1.5)
theta <- rdirichlet(d, theta.ini)   #���[�U�[�g�s�b�N�̃p�����[�^�̏����l
phi.ini <- runif(v, 0.5, 1)
phi <- rdirichlet(k, phi.ini)   #�A�C�e���g�s�b�N�̃p�����[�^�̏����l

##�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
W.SEG <- matrix(0, nrow=sum(w), R/keep)
RATE <- matrix(0, nrow=R/keep, ncol=k)
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XZ))


####MCMC����̂��߂̃f�[�^�̏���####
#ID���쐬
d.id <- rep(1:d, rep(v, d))
w.id <- rep(1:v, d)
ID <- data.frame(d.id, w.id)

#�C���f�b�N�X���쐬
index_w <- list()
index_v <- list()

for(i in 1:d) {index_w[[i]] <- subset(1:length(d.id), d.id==i)}
for(i in 1:v) {index_v[[i]] <- subset(1:length(w.id), w.id==i)}


##�g�s�b�N�����̏����l�𐶐�
Zx <- matrix(0, nrow=nrow(ID), ncol=k)

#���[�U�[���ƂɃA�C�e���̃g�s�b�N�𐶐�
for(i in 1:d){
  theta.m <- matrix(theta[i, ], nrow=k, ncol=v)   #theta���s��`���ɕύX 
  z.rate <- t(phi * theta.m) / matrix(rowSums(t(phi * theta.m)), nrow=v, ncol=k)   #���������v�Z
  Zx[index_w[[i]], ] <- t(apply(cbind(WX[i, ], z.rate), 1, function(x) rmultinom(1, x[1], x[-1])))   #�P��̃g�s�b�N����������
}

#�g�s�b�N�������̏����l���v�Z
#�S�̂ł̃g�s�b�N����
k_sum <- colSums(Zx) 

#���[�U�[���Ƃ̃g�s�b�N����
kw_sum <- matrix(0, nrow=d, ncol=k)
for(i in 1:d) {kw_sum[i, ] <- colSums(Zx[index_w[[i]], ])}

#�A�C�e���̃g�s�b�N����
kv_sum <- matrix(0, nrow=v, ncol=k)
for(i in 1:v) {kv_sum[i, ] <- colSums(Zx[index_v[[i]], ])}


#�g�s�b�N�������x�N�g���`���ɕύX
seg_vec <- unlist(apply(Zx, 1, function(x) rep(1:k, x)))   

#�s��`���ɕύX
seg_mx <- matrix(0, nrow=length(seg_vec), ncol=k)
for(i in 1:nrow(seg_mx)) {seg_mx[i, seg_vec[i]] <- 1}


#�g�s�b�N�����x�N�g����ID���쐬
id_vec11 <- rep(1:d, w)
id_vec12 <- c()
for(i in 1:d) {id_vec12 <- c(id_vec12, rep(1:v, rowSums(Zx[index_w[[i]], ])))}

Z1 <- matrix(0, nrow=length(id_vec11), ncol=k)


#�M�u�X�T���v�����O�p�̃C���f�b�N�X���쐬
index_user <- list()
for(i in 1:d) {index_user[[i]] <- subset(1:length(id_vec11), id_vec11==i)}


##���W�X�e�B�b�N��A���f���̑ΐ��ޓx���`
loglike <- function(b, X, Y){
  
  #�ޓx���`���č��v����
  logit <- X %*% b 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}


####�}���R�t�A�������e�J�����@�Ńg�s�b�N���f��+���W�X�e�B�b�N��A���f���𐄒�####
for(rp in 1:R){
  
  ##�g�s�b�N���T���v�����O
  for(i in 1:d){
    ##�P��̃g�s�b�N���T���v�����O
    for(wd in 1:length(index_user[[i]])){
      index <- index_user[[i]][wd]   #�P��̃C���f�b�N�X
      
      #�g�s�b�N��������P��̃g�s�b�N����菜��
      mx <- seg_mx[index, ]
      k1 <- k_sum - mx
      kw <- kw_sum[i, ] - mx
      kv <- kv_sum[id_vec12[index], ] - mx
    id_vec12
    kv_sum
    
      #�P��̃g�s�b�N�����m�����v�Z
      z_sums <- (kw + alpha01) * (kv + beta01) / (k1 + beta01*v)
      z_rate <- z_sums / sum(z_sums)
      
      #�g�s�b�N���T���v�����O
      Z1 <- t(rmultinom(1, 1, z_rate))
      
      #�f�[�^���X�V
      k_sum <- k1 + Z1
      kw_sum[i, ] <- kw + Z1
      kv_sum[id_vec12[index], ] <- kv + Z1
      seg_mx[index, ] <- Z1
      
    }
  }
  
  ##���W�X�e�B�b�N��A���f���̃p�����[�^���T���v�����O
  #�g�s�b�N���W�v
  ZL <- matrix(0, nrow=d, ncol=k)
  
  for(i in 1:d){
    z_ind <- colSums(seg_mx[index_user[[i]], ])
    ZL[i, ] <- z_ind / sum(z_ind)
  }
  XZ <- cbind(X, ZL[, -k])   #�f�[�^������
  
  #beta�̃T���v�����O
  if(rp %in% c(1, 10, 100, 1000)){
    res <- glm(y ~ ., data=data.frame(XZ[, -1]), family=binomial)
    invBi <- diag(summary(res)[[12]][, 2]^2)
  }
  
  #�����_���E�H�[�N�T���v�����O
  betad <- betas
  betan <- betad + 0.1 * mvrnorm(1, rep(0, length(betad)), invBi)   #�V����beta�������_���E�H�[�N�ŃT���v�����O
  
  #�ΐ��ޓx�Ƒΐ����O���z�̌v�Z
  lognew <- loglike(betan, XZ, y)
  logold <- loglike(betad, XZ, y)
  logpnew <- lndMvn(betan, betaBi, rootBi)
  logpold <- lndMvn(betad, betaBi, rootBi)
  
  #MH�T���v�����O
  alpha <- min(1, exp(lognew + logpnew - logold - logpold))
  if(alpha == "NAN") alpha <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V����beta���̑�
  if(u < alpha){
    betas <- betan
    logl <- lognew
    
    #�����łȂ��Ȃ�beta���X�V���Ȃ�
  } else {
    logl <- logold
  } 
  
  ##�T���v�����O���ʂ�ۑ�
  mkeep <- rp/keep
  if(rp%%keep==0){
    
    #�������̌v�Z
    rate <- colSums(seg_mx)/nrow(seg_mx)
    
    #�T���v�����O���ʂ�ۑ�
    W.SEG[, mkeep] <- seg_mx %*% 1:k 
    RATE[mkeep, ] <- rate
    
    #�T���v�����O�󋵂�\��
    print(rp)
    print(round(rbind(rate1=rate, ratet=z_ratet), 3))
    print(round(c(logl, loglt), 2))
    print(round(rbind(betas, b), 3))
  }
}


