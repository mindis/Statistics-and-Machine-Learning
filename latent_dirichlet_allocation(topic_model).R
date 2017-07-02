#####�g�s�b�N���f��(Latent Dirichlet Allocation Model)#####
library(MASS)
library(lda)
library(RMeCab)
library(gtools)
library(reshape2)
library(plyr)
library(ggplot2)


####�f�[�^�̔���####
#set.seed(423943)
#�f�[�^�̐ݒ�
k <- 5   #�g�s�b�N��
d <- 200   #������
v <- 100   #��b��
w <- 200   #1����������̒P�ꐔ 

#�p�����[�^�̐ݒ�
alpha0 <- runif(k, 0.1, 1.25)   #�����̃f�B���N�����O���z�̃p�����[�^
alpha1 <- rep(0.5, v)   #�P��̃f�B���N�����O���z�̃p�����[�^

#�f�B���N�������̔���
theta0 <- rdirichlet(d, alpha0)   #�����̃g�s�b�N���z���f�B���N���������甭��
phi0 <- rdirichlet(k, alpha1)   #�P��̃g�s�b�N���z���f�B���N���������甭��

#�������z�̗�������f�[�^�𔭐�
WX <- matrix(0, d, v)
ZS <- list()

for(i in 1:d){
  z <- t(rmultinom(w, 1, theta0[i, ]))   #�����̃g�s�b�N���z�𔭐�
  
  zn <- z %*% c(1:k)   #0,1�𐔒l�ɒu��������
  zdn <- cbind(zn, z)   #apply�֐��Ŏg����悤�ɍs��ɂ��Ă���
  
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi0[x[1], ])))   #�����̃g�s�b�N����P��𐶐�
  wdn <- colSums(wn)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  WX[i, ] <- wdn  
  ZS[[i]] <- zdn[, 1]
  print(i)
}
barplot(colSums(z), names.arg=c("seg1", "seg2", "seg3", "seg4", "seg5"))
round(colSums(WX)/sum(WX), 3)   #�P��̏o���p�x


####�}���R�t�A�������e�J�����@�Ńg�s�b�N���f���𐄒�####
##�}���R�t�A�������e�J�����@�̐ݒ�
#�A���S���Y���̐ݒ�
R <- 5000   #�T���v�����O��
keep <- 2   #4���1��̊����ŃT���v�����O���ʂ𗘗p
iter <- 0

#�n�C�p�[�p�����[�^�̎��O���z�̐ݒ�
alpha <- rep(1.5, k)
beta <- rep(100/v, v)

#�p�����[�^�̏����l
theta.ini <- runif(k, 0.5, 2)
phi.ini <- runif(v, 0.5, 1)
theta <-    rdirichlet(d, theta.ini)   #�����g�s�b�N�̃p�����[�^�̏����l
phi <- rdirichlet(k, phi.ini)   #�P��g�s�b�N�̃p�����[�^�̏����l


#�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
W.SEG <- matrix(0, nrow=d*w, ncol=R/(keep*2))


####�M�u�X�T���v�����O�Ńg�s�b�N���f���𐄒�####
id1 <- rep(1:d, v)
id2 <- rep(1:v, rep(d, v))

##�P��p�x�s����x�N�g���ɕϊ�
d.vec <- c() 
v.vec <- c()
for(i in 1:v){
  for(j in 1:d){
    d.vec <- c(d.vec, rep(j, WX[j, i]))
    v.vec <- c(v.vec, rep(i, WX[j, i]))
  }
}

index.vec <- cbind(v.vec, d.vec)
id.w <- cbind(id2, id1)

#���o������݊m���̍s���擾
index.z <- apply(index.vec, 1, function(x) subset(1:nrow(id.w), id.w[, 1]==x[1] & id.w[, 2]==x[2]))

#�g�s�b�Nz���T���v�����O���邽�߂̔z��
theta.seg <- array(0, dim=c(d, v, k))
phi.seg <- array(0, dim=c(d, v, k))
burden.seg <- array(0, dim=c(d, v, k))
zpt <- array(0, dim=c(d, v, k))

#�p�����[�^���T���v�����O���邽�߂̔z��
dir.theta <- matrix(0, nrow=d, ncol=k) 
dir.phi <- matrix(0, nrow=k, ncol=v)

####�M�u�X�T���v�����O####
##�P��̃g�s�b�Nz���T���v�����O
for(rp in 1:R){
  
  burden.sum <- matrix(0, nrow=d, ncol=v) 
  #���݊m���̃p�����[�^���v�Z
  for(i in 1:k){
    theta.seg[, , i] <- matrix(theta[, i], nrow=d, ncol=v)
    phi.seg[, , i] <- matrix(phi[i, ], nrow=d, ncol=v, byrow=T)
    burden.seg[, , i] <- theta.seg[, , i] * phi.seg[, , i]
    burden.sum <- burden.sum + burden.seg[, , i]
  }

  #���݊m�����v�Z
  for(i in 1:k){
    zpt[, , i] <- burden.seg[, , i] / burden.sum
  }

  ZP <- matrix(zpt, nrow=v*d, ncol=k)   #���݊m����P�ꏇ�ɔz��
  ZW <- ZP[index.z, ]

  #�P�ꂲ�ƂɃg�s�b�Nz���T���v�����O
  Z <- t(apply(ZW, 1, function(x) rmultinom(1, 1, x)))
  
  ##�����g�s�b�N���ztheta���T���v�����O
  #�f�B���N�����z�̃p�����[�^������
  for(i in 1:d){
   dir.theta[i, ] <- colSums(Z[d.vec==i, ]) + alpha 
  }
  
  theta <- t(apply(dir.theta, 1, function(x) rdirichlet(1, x)))   #�f�B���N�����z����theta���T���v�����O

  ##�P��g�s�b�N���zphi���T���v�����O
  #�f�B���N�����z�̃p�����[�^������
  for(i in 1:v){
   dir.phi[, i] <- colSums(Z[v.vec==i, ]) + beta[i]
  }
  
  phi <- t(apply(dir.phi, 1, function(x) rdirichlet(1, x)))   #�f�B���N�����z����phi���T���v�����O
  
  ##�n�C�p�[�p�����[�^�̍X�V
  
  
  ##�T���v�����O���ʂ̕ۑ�
  #�T���v�����O��ۑ�����񐔂Ȃ�beta����������
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
  }
  if(rp%%(keep*2)==0){
    W.SEG[, mkeep/2] <- Z %*% 1:k
    print(rp)
  }
}


####���茋�ʂƓK���x####
burnin <- 1500
Rkeep <- R/keep

##�T���v�����O���ʂ̃v���b�g
matplot(t(THETA[1:3, 1, ]), type="l", lty=1, ylab="value")
matplot(t(THETA[4:6, 1, ]), type="l", lty=1, ylab="value")
matplot(t(PHI[1, 1:3, ]), type="l", lty=1, ylab="value")
matplot(t(PHI[2, 1:3, ]), type="l", lty=1, ylab="value")

#theta�̎��㕽�ςƐ^��theta�̔�r
round(THETA1 <- cbind(theta0, rowMeans(THETA[, 1, burnin:Rkeep]), rowMeans(THETA[, 2, burnin:Rkeep]), 
      rowMeans(THETA[, 3, burnin:Rkeep]), rowMeans(THETA[, 4, burnin:Rkeep]), 
      rowMeans(THETA[, 5, burnin:Rkeep])), 3)

#phi�̎��㕽�ςƐ^��phi�̔�r
round(PHI1 <- rbind(phi0, rowMeans(PHI[1, , burnin:Rkeep]), rowMeans(PHI[2, , burnin:Rkeep]),
      rowMeans(PHI[3, , burnin:Rkeep]), rowMeans(PHI[4, , burnin:Rkeep]), 
      rowMeans(PHI[5, , burnin:Rkeep])), 3)

#�g�s�b�N�̎��㕪�z
segment <- apply(W.SEG[, (burnin/2):(Rkeep/2)], 1, table)

