#####�\���v���r�b�g���f��#####
library(MASS)
library(bayesm)
library(condMVNorm)
library(MCMCpack)
library(glmm)
library(lme4)
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
hh <- 1000   #�v���C���[��
member <- 10   #�I���\��
k <- 5   #��A�W���̐�

##ID�̐ݒ�
id <- rep(1:hh, rep(member-1, hh))
c <- rep(1:(member-1), hh)
ID <- data.frame(no=1:(hh*(member-1)), id=id, c=c)
id_r <- matrix(1:(hh*(member-1)), nrow=hh, ncol=member-1, byrow=T)


##�����ϐ��̔���
cont <- 2
X1.cont <- matrix(rnorm(hh*member*2, 0, 1), nrow=hh, ncol=2*member)

bin <- 2
X1.bin <- matrix(0, nrow=hh, ncol=2*member)
for(i in 1:(2*member)){
  X1.bin[, i]  <- rbinom(hh, 1, runif(1, 0.35, 0.6))
}



##�����ϐ����x�N�g���`���̃t�H�[�}�b�g�ɕύX
#ID��ݒ�
id <- rep(1:hh, rep(member-1, hh))
m <- rep(1:(member-1), hh)
ID <- data.frame(no=1:length(id), id=id, m=m)

#�ؕЂ̐ݒ�
p <- c(1, rep(0, member-1))
Pop <- matrix(p, nrow=hh*member, ncol=member-1, byrow=T)
POP <- subset(Pop, rowSums(Pop) > 0)

#�\���ϐ��̐ݒ�
c <- 3 
g1 <- rep(c(rep(1, (member-1)/c), rep(0, 2*(member-1)/c)), hh)
g2 <- rep(c(rep(0, (member-1)/c), rep(1, (member-1)/c), rep(0, (member-1)/c)), hh)
g3 <- rep(c(rep(0, 2*(member-1)/c), rep(1, (member-1)/c)), hh)
G <- cbind(g1, g2)


#�����ϐ��̐ݒ�
#���Ό��p�ɕύX
X1r.cont1 <- X1.cont[, 1:(member-1)] - X1.cont[, member]
X1r.cont2 <- X1.cont[, (member+1):(2*member-1)] - X1.cont[, (2*member)]
X1r.bin1 <- X1.bin[, 1:(member-1)] - X1.bin[, member]
X1r.bin2 <- X1.bin[, (member+1):(2*member-1)] - X1.bin[, (2*member)]

#�x�N�g���`���ɕύX
X1v.cont1 <- as.numeric(t(X1r.cont1))
X1v.cont2 <- as.numeric(t(X1r.cont2))
X1v.bin1 <- as.numeric(t(X1r.bin1))
X1v.bin2 <- as.numeric(t(X1r.bin2))

##�f�[�^������
X <- data.frame(pop=POP, G, c1=X1v.cont1, c2=X1v.cont2, b1=X1v.bin1, b2=X1v.bin2)
XM <- as.matrix(X)

##���U�����U�s��̐ݒ�
corM <- corrM(col=member-1, lower=-0.6, upper=0.70)   #���֍s����쐬
Sigma <- covmatrix(col=member-1, corM=corM, lower=1, upper=1)   #���U�����U�s��
Cov <- Sigma$covariance


##�����f�[�^�𔭐�
##�p�����[�^�̐ݒ�
##�Ó��ȉ����ϐ�����������܂ŌJ��Ԃ�
for(i in 1:10000){
  print(i)
  
  #��A�W���̃p����-�^�̐ݒ�
  b0 <- runif(member-1, 0.1, 1.7)
  b1 <- runif(c-1, 0.2, 1.3)
  b2 <- runif(cont, 0, 1.1)
  b3 <- runif(bin, -0.9, 1.2)
  b <- c(b0, b1, b2, b3)
  beta.t <- b
  
  ##���Ό��p�𔭐�������
  err <- mvrnorm(hh, rep(0, member-1), Cov)   #�덷�\��
  U.mean <- matrix(XM %*% b, nrow=hh, ncol=member-1, byrow=T)   #���Ό��p�̕��ύ\��
  U <- U.mean + err   #�덷�\�����������\��
  
  ##���p�ő剻�����Ɋ�Â��I�������o�[������
  y <- apply(cbind(U, 0), 1, which.max)
  
  #������o�[���K���Ȑl���ɑI�΂��܂Ń��[�v������
  if(sum(y==member) > 15 & sum(y==member) < 100) {break} else {next}
}

#�I�������o�[��0�A1�s��ɕύX
Y <- matrix(0, hh, member)
for(i in 1:hh){
  Y[i, y[i]] <- 1
}

table(y)   #�I���u�����h�̏W�v
round(cbind(Y, U, U.mean), 2)   #���p�ƑI���u�����h���r


####�}���R�t�A�������e�J�����@�ő����v���r�b�g���f���𐄒�####
####MCMC����̂��߂̐��菀��####

##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##���ϗʐ��K���z�̏����t�����Ғl�ƕ��U���v�Z����֐�
cdMVN <- function(mu, Cov, dependent, U){
  
  #���U�����U�s��̃u���b�N�s����`
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent]
  Cov21 <- Cov[-dependent, dependent]
  Cov22 <- Cov[-dependent, -dependent]
  
  #�����t�����U�Ə����t�����ς��v�Z
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mu[, dependent] + t(CDinv %*% t(U[, -dependent] - mu[, -dependent]))   #�����t�����ς��v�Z
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #�����t�����U���v�Z
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- array(0, dim=c(R/keep))   #�ΐ��ޓx�̕ۑ��p


##�����ϐ��𑽎����z��
X.array <- array(0, dim=c(member-1, ncol(X), hh))
for(i in 1:hh){
  X.array[, , i] <- XM[ID[, 2]==i, ]
}
YX.array <- array(0, dim=c(member-1, ncol(X)+1, hh))


#����v���Z�X�̊i�[�z��
UM <- matrix(0, nrow=hh, ncol=member-1)
util.M <- matrix(0, nrow=hh, ncol=member-1)   

##���O���z�̐ݒ�
nu <- member   #�t�E�B�V���[�g���z�̎��R�x
V <- solve((1/10)*diag(member-1))    #�t�E�B�V���[�g���z�̃p�����[�^
Deltabar <- rep(0, ncol(X))  #��A�W���̕��ς̎��O���z
Adelta <- solve(100 * diag(rep(1, ncol(X))))   #��A�W���̎��O���z�̕��U

##�T���v�����O���ʂ̕ۑ��p�z��
Util <- array(0, dim=c(hh, member-1, R/keep))
BETA <- matrix(0, nrow=R/keep, length(b))
SIGMA <- matrix(0, nrow=R/keep, ncol=(member-1)^2)

##�����l�̐ݒ�
#��A�W���̏����l
oldbeta <- c(runif(member-1, 0, 2), runif(c-1, 0, 1.5), runif(2, 0, 2.0), runif(2, -2, 2))   


#���U�����U�s��̏����l
corM.f <- corrM(col=member-1, lower=0, upper=0)   #���֍s����쐬
Sigma.f <- covmatrix(col=member-1, corM=corM.f, lower=1, upper=1)   #���U�����U�s��
oldcov <- Sigma.f$covariance

#���p�̕��ύ\���̏����l
old.utilm <- matrix(XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)

#���p�̏����l
old.util <- old.utilm + mvrnorm(nrow(old.utilm), rep(0, member-1), oldcov)


####�}���R�t�A�������e�J�����@�ő����v���r�b�g���f���𐄒�####
for(rp in 1:R){
  
  ##�I�����ʂƐ����I�Ȑ��݌��p�𔭐�������
  #�����t�����Ғl�Ə����t�����U���v�Z
  S <- rep(0, member-1)
  
  for(j in 1:(member-1)){
    MVR <- cdMVN(mu=old.utilm, Cov=oldcov, dependent=j, U=old.util)   #�����t�����z���v�Z
    UM[, j] <- MVR$CDmu   #�����t�����Ғl�����o��
    S[j] <- sqrt(MVR$CDvar)    #�����t�����U�����o��
    
    #���ݕϐ��𔭐�������
    #�ؒf�̈�̐ݒ�
    max.u <- apply(cbind(old.util[, -j], 0), 1, max)
    max.u <- ifelse(y==member, 0, max.u)
    
    #�ؒf���K���z�����ݕϐ��𔭐�
    old.util[, j] <- ifelse(y==j, rtnorm(mu=UM[, j], sigma=S[j], a=max.u, b=100), 
                            rtnorm(mu=UM[, j], sigma=S[j], a=-100, b=max.u))
    old.util[, j] <- ifelse(is.infinite(old.util[, j]), ifelse(y==j, max.u + runif(1), max.u - runif(1)), old.util[, j])
  }
  util.v <- as.numeric(t(old.util))
  
  ##beta�̕��z�̃p�����[�^�̌v�Z��mcmc�T���v�����O
  #z.vec��X.vec���������đ������z��ɕύX
  YX.bind <- cbind(util.v, XM)
  for(i in 1:hh){
    YX.array[, , i] <- YX.bind[id_r[i, ], ]
  }
  
  ##��A���f���̃M�u�X�T���v�����O��beta��sigma�𐄒�
  #beta�̃M�u�X�T���v�����O
  invcov <- solve(oldcov)
  xvx.vec <- rowSums(apply(X.array, 3, function(x) t(x) %*% invcov %*% x))
  XVX <- matrix(xvx.vec, nrow=ncol(X), ncol=ncol(X), byrow=T)
  XVY <- rowSums(apply(YX.array, 3, function(x) t(x[, -1]) %*% invcov %*% x[, 1]))
  
  #beta�̕��z�̕��U�����U�s��̃p�����[�^
  inv_XVX <- solve(XVX + Adelta)
  
  #beta�̕��z�̕��σp�����[�^
  B <- inv_XVX %*% (XVY + Adelta %*% Deltabar)   #beta�̕���
  b1 <- as.numeric(B)
  
  #���ϗʐ��K���z�����A�W�����T���v�����O
  oldbeta <- mvrnorm(1, b1, inv_XVX)
  
  ##Cov�̕��z�̃p�����[�^�̌v�Z��mcmc�T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^���v�Z
  R.error <- matrix(util.v - XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)
  IW.R <- V + matrix(rowSums(apply(R.error, 1, function(x) x %*% t(x))), nrow=member-1, ncol=member-1, byrow=T)
  
  #�t�E�B�V���[�g���z�̎��R�x���v�Z
  Sn <- nu + hh
  
  #�t�E�B�V���[�g���z����Cov���T���v�����O
  Cov_hat <- rwishart(Sn, solve(IW.R))$IW
  oldcov <- cov2cor(Cov_hat)
  
  ##���݌��p�ƃp�����[�^���X�V
  #���݌��p�Ɛ��݌��p�̕��ς��X�V
  old.utilm <- matrix(XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)
  
  ##�T���v�����O���ʂ�ۑ�
  if(rp%%keep==0){
    print(rp)
    mkeep <- rp/keep
    Util[, , mkeep] <- old.util
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    print(round(cbind(oldcov, Cov), 2))
    print(round(rbind(oldbeta, beta.t), 2))
  }
}


####�֐��Ő���####
Data1 <- list(p=member, y=y, X=XM)
Mcmc1 <- list(R=10000, keep=2)

#�����v���r�b�g���f���𐄒�
out <- rmnpGibbs(Data=Data1,Mcmc=Mcmc1)
BETA.out <- out$betadraw
SIGMA.out <- out$sigmadraw

####���茋�ʂ̗v��ƓK���x�̊m�F####
burnin <- 5000/keep   #�o�[���C������

##�T���v�����O���ʂ�����
#��A�W���̃v���b�g
matplot(BETA[, 1:3], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 4:6], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 7:9], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 10:12], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 13:15], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")

#���U�����U�s��̉���
matplot(SIGMA[, 1:4], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 5:8], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 9:12], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 13:16], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")


##����l�̎��㕽�ς̔�r
#beta�̗v�񓝌v��
round(colMeans(BETA.out[burnin:nrow(BETA.out), ] / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(�֐�����)�̎��㕽��
round(colMeans(BETA[burnin:nrow(BETA), ]), 3)   #beta�̎��㕽��
round(beta.t, 3)   #�^�̒l
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5�����ʓ_
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95�����ʓ_
round(apply(BETA[burnin:nrow(BETA), ], 2, sd), 2)   #����W���΍�

#sigma�̗v�񓝌v��
round(colMeans(SIGMA.out[burnin:nrow(SIGMA.out), ]  / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(�֐�����)�̎��㕽��
round(colMeans(SIGMA[burnin:nrow(SIGMA), ]), 3)   #beta�̎��㕽��
round(as.numeric(Cov), 3)   #�^�̒l
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.05)), 2)   #5�����ʓ_
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.95)), 2)   #95�����ʓ_
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, sd), 2)   #����W���΍�