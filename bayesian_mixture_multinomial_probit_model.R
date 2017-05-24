#####�x�C�W�A�����������v���r�b�g���f��#####
library(MASS)
library(bayesm)
library(condMVNorm)
library(MCMCpack)
library(gtools)
library(MNP)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####���ϗʐ��K�����𔭐�������֐�####
##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

#���ϗʉ�A���f���̑��֍s����쐬
##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

##���֍s�񂩂番�U�����U�s����쐬����֐����`
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #�ŗL�l�����ŋ����I�ɐ���l�s��ɏC������
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}


####�f�[�^�̔���####
#set.seed(8437)
##�f�[�^�̐ݒ�
seg <- 3   #�Z�O�����g��
hh <- 500   #�Z�O�����g���Ƃ̍w���Ґ�
H <- seg*hh
choise <- 5   #�I���\��
st <- 5   #��u�����h
k <- 5   #��A�W���̐�
g <- rep(1:seg, rep(hh, seg))

##�����ϐ��̔���
#�ʏ퉿�i�̔���
PRICE <- matrix(runif(H*choise, 0.7, 1), nrow=H, ncol=choise)   

#�f�B�X�J�E���g���̔���
DISC <- matrix(runif(H*choise, 0, 0.3), nrow=H, ncol=choise)

#���ʒ�̔���
DISP <- matrix(0, nrow=H, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(H, 1, r)
}

#���ʃL�����y�[���̔���
CAMP <- matrix(0, nrow=H, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(H, 1, r)
}

##�Z�O�����g���Ƃ̕��U�����U�s��̐ݒ�
Cov <- list()
for(i in 1:seg){
  corM <- corrM(col=choise-1, lower=-0.55, upper=0.60)   #���֍s����쐬
  Sigma <- covmatrix(col=choise-1, corM=corM, lower=1, upper=4)   #���U�����U�s��
  Cov[[i]] <- Sigma$covariance
}


##�Z�O�����g�̐ݒ�
z.seg <- rep(1:seg, rep(hh, seg))

##�p�����[�^�̐ݒ�
beta1 <- c(-9.0, -7.0, -5.0)   #���i�̃p�����[�^
beta2 <- c(8.2, 5.5, 9.2)   #�������̃p�����[�^
beta3 <- c(0.7, 2.3, 3.2)   #���ʒ�̃p�����[�^
beta4 <- c(0.5, 1.2, 3.3)   #�L�����y�[���̃p�����[�^

#�u�����h1�`4�̑��΃x�[�X�̔���
beta0 <- matrix(0, nrow=seg, ncol=choise-1)
for(i in 1:seg){
  beta0[i, ] <- c(runif(1, -1.2, 2.8), runif(1, -1.5, 3.6), runif(1, -1.0, 4.3), runif(1, -1.5, 4.0))   
}

#��A�W��������
betat <- matrix(0, nrow=seg, ncol=ncol(beta0)+k-1)
for(i in 1:seg){
  betat[i, ] <- c(beta0[i, ], beta1[i], beta2[i], beta3[i], beta4[i])
}
round(betat, 2)


##���Ό��p�𔭐������A�I�����ꂽ�u�����h������
#��u�����h�Ƃ̑��ΐ����ϐ�
PRICE.r <- PRICE[, -5] - PRICE[, 5]
DISC.r <- DISC[, -5] - DISC[, 5]
DISP.r <- DISP[, -5] - DISP[, 5]
CAMP.r <- CAMP[, -5] - CAMP[, 5]

#���Ό��p�̕��ύ\���𔭐�
U.mean <- matrix(0, nrow=H, ncol=choise-1)
for(s in 1:seg){
  index <- subset(1:length(g), g==s)
  for(b in 1:ncol(beta0)){
    U.mean[index, b] <- beta0[s, b] + PRICE.r[g==s, b]*beta1[s] + DISC.r[g==s, b]*beta2[s] + 
                        DISP.r[g==s, b]*beta3[s] + CAMP.r[g==s, b]*beta4[s]
  }
}

#�덷�\�������������Ό��p
U <- matrix(0, nrow=H, ncol=choise-1)
for(s in 1:seg){
  index <- subset(1:length(g), g==s)
  U[index, ] <- t(apply(U.mean[g==s, ], 1, function(x) mvrnorm(1, x, Cov[[s]])))
}

#���p�ő剻�����Ɋ�Â��I���u�����h������
Y <- apply(U, 1, function(x) ifelse(max(x) < 0, 5, which.max(x)))

#�w����0�A1�s��ɕύX
BUY <- matrix(0, nrow=H, ncol=choise)
for(i in 1:H){
  BUY[i, Y[i]] <- 1
}

table(Y)   #�I���u�����h�̏W�v
round(data.frame(Y, U=U), 1)   #���p�ƑI���u�����h���r


####�}���R�t�A�������e�J�����@�ő����v���r�b�g���f���𐄒�####
####MCMC����̂��߂̐��菀��####

##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##MCMC�A���S���Y���̐ݒ�
R <- 20000
keep <- 2

##��A���f���𐄒肷�邽�߂ɐ����ϐ����x�N�g���`���ɕύX�ݒ�
#�ؕЂ̐ݒ�
p <- c(1, rep(0, choise-1))
bp <- matrix(p, nrow=H*choise, ncol=choise-1, byrow=T)
BP <- subset(bp, rowSums(bp) > 0)

#�����ϐ��̐ݒ�
PRICE.v <- as.numeric(t(PRICE.r))
DISC.v <- as.numeric(t(DISC.r))
DISP.v <- as.numeric(t(DISP.r))
CAMP.v <- as.numeric(t(CAMP.r))

X <- data.frame(BP=BP, PRICE.v, DISC.v, DISP.v, CAMP.v)   #�f�[�^�̌���
XM <- as.matrix(X)   #�f�[�^�`�����s��ɕϊ�

#ID�̃C���f�b�N�X��ݒ�
index <- rep(1:H, rep(choise-1, H))   

#X���x�N�g�������Ă���
XV <- matrix(0, nrow=H, ncol=(choise-1)*ncol(XM))
for(i in 1:H){
  XV[i, ] <- as.numeric(XM[index==i, ])
}

##���O���z�̐ݒ�
a <- rep(10, seg)   #�f�B�N�������O���z

##�T���v�����O���ʂ̕ۑ��p�z��
Util <- array(0, dim=c(hh, choise-1, R/keep))
BETA <- matrix(0, nrow=R/keep, length(beta0)+k-1)
SIGMA <- array(0, dim=c(choise-1, choise-1, R/keep))
Z <- matrix(0, nrow=R/keep, ncol=hh)
THETA <- matrix(0, nrow=R/keep, ncol=seg)

##�����p�����[�^�̐ݒ�
#�����v���r�b�g���f���ŏ����l��ݒ�
data1 <- list(p=choise, y=Y[g==1], X=XM[1:(length(g[g==1])*(choise-1)), ])
mcmc1 <- list(R=10000, keep=1)

out <- rmnpGibbs(Data=data1, Mcmc=mcmc1)   #�����v���r�b�g���f���𐄒�
betaf <- colMeans(out$betadraw[5000:nrow(out$betadraw), ])   
sigmaf <- matrix(colMeans(out$sigmadraw[5000:nrow(out$sigmadraw), ]), nrow=choise-1, ncol=choise-1)

#�Z�O�����g�ʂɏ����l��ݒ�
#beta�̏����l
betaf.M <- matrix(betaf, nrow=seg, ncol=length(betaf), byrow=T)
betaold <- betaf.M + matrix(runif(seg*length(betaf), -0.5, 0.5), nrow=seg, ncol=length(betaf))

#sigma�̏����l
sigmaold <- list()
for(i in 1:seg){
  sigmaold[[i]] <- sigmaf
}

#theta�̏����l
theta <- c(0.3, 0.3, 0.4)


####�}���R�t�A�������e�J�����@�ō��������v���r�b�g���f���𐄒�####
##�Z�O�����g�ʂɑ����v���r�b�g���f���𐄒�
##����̂��߂̏���
#�����ϐ��̃Z�O�����g�̎w���ϐ����쐬
gx <- matrix(0, nrow=H, ncol=choise-1)
for(s in 1:(choise-1)){
  gx[, s] <- g
}
gx <- as.numeric(t(gx))

burnin <- 3000   #�o�[���C������
R <- 10000

#���茋�ʂ̊i�[�p�z��
beta.s <- matrix(0, nrow=seg, ncol=length(betaf))
sigma.s <- list() 
BETA.S <- list()
SIGMA.S <- list()
  
##�����v���r�b�g���f�����M�u�X�T���v�����O�ŃZ�O�����g�ʂɐ���
for(s in 1:seg){
#�f�[�^�̐ݒ�
  data1 <- list(p=choise, y=Y[g==s], X=XM[gx==s, ])
  mcmc1 <- list(R=10000, keep=1)
  
  #�����v���r�b�g���f���𐄒�
  out <- rmnpGibbs(Data=data1, Mcmc=mcmc1)   #�����v���r�b�g���f���𐄒�
  beta.s[s, ] <- colMeans(out$betadraw[burnin:nrow(out$betadraw), ])   
  sigma.s[[s]] <- matrix(colMeans(out$sigmadraw[burnin:nrow(out$sigmadraw), ]), nrow=choise-1, ncol=choise-1)
  BETA.S[[s]] <- out$betadraw
  SIGMA.S[[s]] <- out$sigmadraw
}

####���������v���r�b�g���f���̐��茋�ʂƗv��####
##�p�����[�^�̎��ʐ����m��
BETA.SI <- list()
SIGMA.SI <- list()

#�����U�s���(1, 1)�v�f��1�ɌŒ肷�鐧��
for(s in 1:seg){
  BETA.SI[[s]] <- BETA.S[[s]] / matrix(SIGMA.S[[s]][, 1], nrow=R, ncol=ncol(beta.s))
  SIGMA.SI[[s]] <- SIGMA.S[[s]] / matrix(SIGMA.S[[s]][, 1], nrow=R, ncol=(choise-1)^2)
}

##�p�����[�^�̗v��
round(colMeans(BETA.SI[[1]][burnin:R, ]), 3)

  
SIGMA.S[[1]]

matrix(SIGMA.S[[1]])[, 1] / 

round(beta.s, 3)
round(betat, 3)
sigma.s
Cov

##���ݕϐ�Z�̃T���v�����O
#�����v���r�b�g���f���̊m���̌v�Z
Pr <- list()
for(s in 1:seg){
  Pr[[s]] <- t(apply(XV, 1, function(x) abs(mnpProb(betaold[s, ] / sigmaold[[s]][1, 1], 
                                                    sigmaold[[s]] / sigmaold[[s]][1, 1], 
                                                    matrix(x, nrow=choise-1, ncol=ncol(XM)), r=25))))
}

#�����ϐ�Y�ƑΉ�����Z�O�����g���Ƃ̊m�����v�Z
Pr.y <- matrix(0, nrow=H, ncol=seg)
for(s in 1:seg){0
  Pr.y[, s] <- rowSums(Pr[[s]] * BUY)
}

#���ݕϐ�Z�𔭐�
z1 <- Pr.y * matrix(theta, nrow=H, ncol=length(theta), byrow=T)
zp <- z1 / rowSums(z1)
z <- t(apply(zp, 1, function(x) rmultinom(1, 1, x)))
Z <- apply(z, 1, which.max)

##�M�u�X�T���v�����O�ő����v���r�b�g���f���𐄒�
#�ŗL�l�����ŋ����I�ɐ���l�s��ɏC������
for(i in 1:seg){
  UDU <- eigen(sigmaold[[s]])
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  sigmaold[[s]] <- vec %*% diag(D) %*% t(vec)
}

for(s in 1:seg){
  index.z <- subset(1:length(Z), Z==s)
  XZ <- XM[index %in% index.z, ]
  
  #�f�[�^�̐ݒ�
  Data1 <- list(p=choise, y=Y[index.z], X=XZ)
  Mcmc1 <- list(beta0=betaold[s, ], sigma0=sigmaold[[s]], R=1000, keep=1)
  
  #�����v���r�b�g���f�����M�u�X�T���v�����O
  out=rmnpGibbs(Data=Data1,Mcmc=Mcmc1)
  betaold[s, ] <- colMeans(out$betadraw[500:1000, ])
  sigmaold[[s]] <- matrix(colMeans(out$sigmadraw[500:1000, ]), nrow=choise-1, ncol=choise-1)
}

##theta���M�u�X�T���v�����O
theta <- table(Z) / sum(table(Z))
round(Pr.y[1:20, ], 3)
theta
betaold
round(betat, 2)


Pr <- list()
for(s in 1:seg){
  Pr[[s]] <- t(apply(XV, 1, function(x) abs(mnpProb(betat[s, ] / Cov[[s]][1, 1], 
                                                    Cov[[s]] / Cov[[s]][1, 1], 
                                                    matrix(x, nrow=choise-1, ncol=ncol(XM)), r=25))))
}

#�����ϐ�Y�ƑΉ�����Z�O�����g���Ƃ̊m�����v�Z
Pr.y <- matrix(0, nrow=H, ncol=seg)
for(s in 1:seg){
  Pr.y[, s] <- rowSums(Pr[[s]] * BUY)
}

cbind(Y, round(Pr.y, 2))

cbind(Y, round(Pr[[1]], 2), round(Pr[[2]], 2))

betat
betaold