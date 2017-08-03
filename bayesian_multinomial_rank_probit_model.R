#####����-�����N�����v���r�b�g���f��#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(gtools)
library(MNP)
library(reshape2)
library(dplyr)
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
hh <- 1500   #�T���v����
member <- 10   #�I���\�����o�[��
k <- 3   #3�ʂ܂őI��

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
X <- data.frame(pop=POP, c1=X1v.cont1, c2=X1v.cont2, b1=X1v.bin1, b2=X1v.bin2)
XM <- as.matrix(X)


####�����ϐ��𔭐�####
##�����N�f�[�^�𔭐�������
##���U�����U�s��̐ݒ�
corM <- corrM(col=member-1, lower=-0.6, upper=0.70)   #���֍s����쐬
Sigma <- covmatrix(col=member-1, corM=corM, lower=1, upper=1)   #���U�����U�s��
Cov <- Sigma$covariance

##�Ó��ȃ����L���O����������܂ŌJ��Ԃ�
for(i in 1:10000){
  print(i)
  
  #��A�W���̃p����-�^�̐ݒ�
  a0 <- runif(member-1, 0.2, 3.1)
  a1 <- runif(cont, 0, 1.1)
  a2 <- runif(bin, -0.9, 1.2)
  a <- c(a0, a1, a2)
  alpha.t <- a
  
  ##���Ό��p�𔭐�������
  err1 <- mvrnorm(hh, rep(0, member-1), Cov)   #�덷�\��
  U1.mean <- matrix(XM %*% a, nrow=hh, ncol=member-1, byrow=T)   #���Ό��p�̕��ύ\��
  U1 <- U1.mean + err1   #�덷�\�����������\��
  
  ##���p�ő剻�����Ɋ�Â����Ώ��ʂ�����
  Rank.full <- t(apply(cbind(U1, 0), 1, function(x) order(x, decreasing=TRUE)))
  Rank <- Rank.full[, 1:3]
  
  #������o�[���K���Ȑl���ɑI�΂��܂Ń��[�v������
  if(sum(Rank.full[, 1:k]==member) > 30 & sum(Rank.full[, 1:k]==member) < 150) {break} else {next}
}

#�����������f�[�^�̊m�F
Rank.full
apply(Rank.full, 2, table) #���ʂ��Ƃ̏W�v


##�����f�[�^�𔭐�
##�p�����[�^�̐ݒ�
##�Ó��ȉ����ϐ�����������܂ŌJ��Ԃ�
for(i in 1:10000){
  print(i)
  
  #��A�W���̃p����-�^�̐ݒ�
  b <- a + runif(length(a), -0.4, 0.4)
  beta.t <- b
  
  ##���Ό��p�𔭐�������
  err2 <- mvrnorm(hh, rep(0, member-1), Cov)   #�덷�\��
  U2.mean <- matrix(XM %*% b, nrow=hh, ncol=member-1, byrow=T)   #���Ό��p�̕��ύ\��
  U2 <- U2.mean + err2   #�덷�\�����������\��
  
  ##���p�ő剻�����Ɋ�Â��I�������o�[������
  y <- apply(cbind(U2, 0), 1, which.max)
  
  #������o�[���K���Ȑl���ɑI�΂��܂Ń��[�v������
  if(sum(y==member) > 15 & sum(y==member) < 100) {break} else {next}
}

#�I�������o�[��0�A1�s��ɕϊ�
Y <- matrix(0, nrow=hh, ncol=member)
for(i in 1:hh){
  Y[i, y[i]] <- 1
}

table(y)   #�I�������o�[�̒P���W�v


####�}���R�t�A�������e�J�����@�ő���-�����N�v���r�b�g���f���𐄒�####
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

##�f�[�^�̐ݒ�
#�����ϐ��̑������z��
X.array <- array(0, dim=c(member-1, ncol(XM), hh))
for(i in 1:hh){
  X.array[, , i] <- XM[ID$id==i, ]
}
YX1.array <- array(0, dim=c(member-1, ncol(XM)+1, hh))
YX2.array <- array(0, dim=c(member-1, ncol(XM)+1, hh))

#ID�̐ݒ�
id_r <- matrix(1:nrow(XM), nrow=hh, ncol=member-1, byrow=T)

##����v���Z�X�̊i�[�p�z��
UM1 <- matrix(0, nrow=hh, ncol=member-1)
UM2 <- matrix(0, nrow=hh, ncol=member-1)
util.M1 <- matrix(0, nrow=hh, ncol=member-1)
util.M2 <- matrix(0, nrow=hh, ncol=member-1)

##���O���z�̐ݒ�
nu <- member   #�t�E�B�V���[�g���z�̎��R�x
V <- solve(0.1*diag(member-1))    #�t�E�B�V���[�g���z�̃p�����[�^
Deltabar <- rep(0, ncol(XM))  #��A�W���̕��ς̎��O���z
Adelta <- solve(100 * diag(rep(1, ncol(XM)))) #��A�W���̎��O���z�̕��U

##�T���v�����O���ʂ̕ۑ��p�z��
Util1 <- array(0, dim=c(hh, member-1, R/keep))
Util2 <- array(0, dim=c(hh, member-1, R/keep))
ALPHA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
THETA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
SIGMA <- matrix(0, nrow=R/keep, ncol=(member-1)^2)


##�����l�̐ݒ�
oldalpha <- c(runif(member-1, 0, 3), runif(ncol(XM)-(member-1), 0, 1.0))  
oldbeta <- oldalpha


#���U�����U�s��̏����l
corM.f <- corrM(col=member-1, lower=0, upper=0)   #���֍s����쐬
Sigma.f <- covmatrix(col=member-1, corM=corM.f, lower=1, upper=1)   #���U�����U�s��
oldcov1 <- Sigma.f$covariance
oldcov2 <- oldcov1


#���p�̕��ύ\���̏����l
old.utilm1 <- matrix(XM %*% oldalpha, nrow=hh, ncol=member-1, byrow=T)
old.utilm2 <- matrix(XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)

#���p�̏����l
old.util1 <- old.utilm1 + mvrnorm(nrow(old.utilm1), rep(0, member-1), oldcov1)
old.util2 <- old.utilm2 + mvrnorm(nrow(old.utilm2), rep(0, member-1), oldcov2)


####�}���R�t�A�������e�J�����@�ő���-�����N�v���r�b�g���f���𐄒�####
for(rp in 1:R){
  
  ##���ʑI�����ʂƐ����I�Ȑ��݌��p�𔭐�������
  #���ϗʐ��K���z�̏����t�����Ғl�Ə����t�����U���v�Z
  S1 <- rep(0, member-1)
  
  for(j in 1:(member-1)){
    MVR1 <- cdMVN(mu=old.utilm1, Cov=oldcov1, dependent=j, U=old.util1)   #�����t�����Ғl�Ə����t�����U���v�Z
    UM1[, j] <- MVR1$CDmu   #�����t�����Ғl�����o��
    S1[j] <- sqrt(MVR1$CDvar)   #�����t�����U�����o��

    
    #���ݕϐ��𔭐�������
    #�ؒf�̈�̐ݒ�
    rank.u <- t(apply(cbind(old.util1[, -j], 0), 1, function(x) sort(x, decreasing=TRUE)))[, 1:3]
    rank.u <- ifelse(Rank==member, 0, rank.u)
  
    #�ؒf���K���z�����ݕϐ��𔭐�
    old.util1[, j] <- ifelse(Rank[, 1]==j, rtnorm(mu=UM1[, j], S1[j], a=rank.u[, 1], b=100), 
                            ifelse(Rank[, 2]==j, rtnorm(mu=UM1[, j], S1[j], a=rank.u[, 2], b=rank.u[, 1]),
                                   ifelse(Rank[, 3]==j, rtnorm(mu=UM1[, j], S1[j], a=rank.u[, 3], b=rank.u[, 2]),
                                          rtnorm(mu=UM1[, j], sigma=S1[j], a=-100, b=rank.u[, 3]))))
  }
  
  util1.v <- as.numeric(t(old.util1))   #�������������݌��p���x�N�g���ɕϊ�
  
  
  ##beta�̕��z�̃p�����[�^�̌v�Z��MCMC�T���v�����O
  #z.vec��X.vec���������đ������z��ɕύX
  YX1.bind <- cbind(util1.v, XM)
  for(i in 1:hh){
    YX1.array[, , i] <- YX1.bind[id_r[i, ], ]
  }
  
  ##�M�u�X�T���v�����O��lalpha��sigma�𐄒�
  #alpha�̃M�u�X�T���v�����O
  invcov1 <- solve(oldcov1)
  xvx1.vec <- rowSums(apply(X.array, 3, function(x) t(x) %*% invcov1 %*% x))
  XVX1 <- matrix(xvx1.vec, nrow=ncol(XM), ncol=ncol(XM), byrow=T)
  XVY1 <- rowSums(apply(YX1.array, 3, function(x) t(x[, -1]) %*% invcov1 %*% x[, 1]))
  
  #alpha�̕��z�̕��U�����U�s��̃p�����[�^
  inv_XVX1 <- solve(XVX1 + Adelta)
  
  #alpha�̕��z�̕��σp�����[�^
  A <- inv_XVX1 %*% (XVY1 + Adelta %*% Deltabar)   #alpha�̕��� 
  a1 <- as.numeric(A)

  #���ϗʐ��K���z�����A�W�����T���v�����O
  oldalpha <- mvrnorm(1, a1, inv_XVX1)
  
  
  ##Cov�̕��z�̃p�����[�^�̌v�Z��mcmc�T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^���v�Z
  R.error1_M <- matrix(util1.v - XM %*% oldalpha, nrow=hh, ncol=member-1, byrow=T)
  R.error1 <- matrix(rowSums(apply(R.error1_M, 1, function(x) x %*% t(x))), nrow=member-1, ncol=member-1, byrow=T)
  IW.R1 <- V + R.error1
  
  #�t�E�B�V���[�g���z�̎��R�x���v�Z
  Sn1 <- nu + hh
  
  #�t�E�B�V���[�g���z����Cov���T���v�����O
  Cov_hat1 <- rwishart(Sn1, solve(IW.R1))$IW
  oldcov1 <- cov2cor(Cov_hat1)
  
  
  ##�I�����ʂƐ����I�Ȑ��݌��p�𔭐�������
  #�����t�����Ғl�Ə����t�����U���v�Z
  S2 <- rep(0, member-1)
  
  for(j in 1:(member-1)){
    MVR2 <- cdMVN(mu=old.utilm2, Cov=oldcov2, dependent=j, U=old.util2)   #�����t�����z���v�Z
    UM2[, j] <- MVR2$CDmu   #�����t�����Ғl�����o��
    S2[j] <- sqrt(MVR2$CDvar)    #�����t�����U�����o��
    
    #���ݕϐ��𔭐�������
    #�ؒf�̈�̐ݒ�
    max.u <- apply(cbind(old.util2[, -j], 0), 1, max)
    max.u <- ifelse(y==member, 0, max.u)
    
    #�ؒf���K���z�����ݕϐ��𔭐�
    old.util2[, j] <- ifelse(y==j, rtnorm(mu=UM2[, j], sigma=S2[j], a=max.u, b=100), 
                            rtnorm(mu=UM2[, j], sigma=S2[j], a=-100, b=max.u))
    old.util2[, j] <- ifelse(is.infinite(old.util2[, j]), ifelse(Y==j, max.u + runif(1), max.u - runif(1)), old.util2[, j])
  }
  
  util2.v <- as.numeric(t(old.util2))   #�������������݌��p���x�N�g���ɕϊ�
  
  
  ##beta�̕��z�̃p�����[�^�̌v�Z��mcmc�T���v�����O
  #z.vec��X.vec���������đ������z��ɕύX
  YX2.bind <- cbind(util2.v, XM)
  for(i in 1:hh){
    YX2.array[, , i] <- YX2.bind[id_r[i, ], ]
  }
  
  ##��A���f���̃M�u�X�T���v�����O��beta��sigma�𐄒�
  #beta�̃M�u�X�T���v�����O
  invcov2 <- solve(oldcov2)
  xvx2.vec <- rowSums(apply(X.array, 3, function(x) t(x) %*% invcov2 %*% x))
  XVX2 <- matrix(xvx2.vec, nrow=ncol(XM), ncol=ncol(XM), byrow=T)
  XVY2 <- rowSums(apply(YX2.array, 3, function(x) t(x[, -1]) %*% invcov2 %*% x[, 1]))
  
  #alpha�̕��z�̕��U�����U�s��̃p�����[�^
  inv_XVX2 <- solve(XVX2 + Adelta)
  inv_XVX3 <- solve(XVX1/2 + XVX2/2 + Adelta)
  
  #beta�̕��z�̕��σp�����[�^
  B <- inv_XVX2 %*% (XVY2 + Adelta %*% Deltabar)
  C <- inv_XVX2 %*% (XVY2/2 + XVY1/2 + Adelta %*% Deltabar)
  b1 <- as.numeric(B)
  c1 <- as.numeric(C)
  
  #���ϗʐ��K���z�����A�W�����T���v�����O
  oldbeta <- mvrnorm(1, b1, inv_XVX2)   #beta�̃T���v�����O
  oldtheta <- mvrnorm(1, c1, inv_XVX3)   #alpha��beta�̌�����A�W��theta�̃T���v�����O
  
  
  ##Cov�̕��z�̃p�����[�^�̌v�Z��mcmc�T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^���v�Z
  R.error2_M <- matrix(util2.v - XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)
  R.error2 <- matrix(rowSums(apply(R.error2_M, 1, function(x) x %*% t(x))), nrow=member-1, ncol=member-1, byrow=T)
  IW.R2 <- IW.R1 + R.error2 
  
  #�t�E�B�V���[�g���z�̎��R�x���v�Z
  Sn2 <- nu + hh
  
  #�t�E�B�V���[�g���z����Cov���T���v�����O
  Cov_hat2 <- rwishart(Sn2, solve(IW.R2))$IW
  oldcov2 <- cov2cor(Cov_hat2)
  
  
  ##���݌��p�ƃp�����[�^���X�V
  old.utilm1 <- matrix(XM %*% oldalpha, nrow=hh, ncol=member-1, byrow=T)
  old.utilm2 <- matrix(XM %*% oldbeta, nrow=hh, ncol=member-1, byrow=T)
  
  
  ##�T���v�����O���ʂ�ۑ�
  if(rp%%keep==0){
    print(rp)
    mkeep <- rp/keep
    Util1[, , mkeep] <- old.util1
    Util2[, , mkeep] <- old.util2
    ALPHA[mkeep, ] <- oldalpha
    BETA[mkeep, ] <- oldbeta
    THETA[mkeep, ] <- oldtheta
    SIGMA[mkeep, ] <- as.numeric(oldcov2)
    
    print(round(cbind(oldcov2, Cov), 2))
    print(round(rbind(oldalpha, alpha.t), 2))
    print(round(rbind(oldbeta, beta.t), 2))
    print(round(oldtheta, 2))
  }
}

####���茋�ʂƗv��####
##�T���v�����O���ʂ��v���b�g
#��A�W���̃v���b�g
matplot(ALPHA[, 1:4], type="l", ylab="alpha1�̉�A�W��", xlab="�T���v�����O��")
matplot(ALPHA[, 5:9], type="l", ylab="alpha2�̉�A�W��", xlab="�T���v�����O��")
matplot(ALPHA[, 10:13], type="l", ylab="alpha3�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA[, 1:4], type="l", ylab="beta1�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA[, 5:9], type="l", ylab="beta2�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA[, 10:13], type="l", ylab="beta3�̉�A�W��", xlab="�T���v�����O��")
matplot(THETA[, 1:4], type="l", ylab="theta1�̉�A�W��", xlab="�T���v�����O��")
matplot(THETA[, 5:9], type="l", ylab="theta2�̉�A�W��", xlab="�T���v�����O��")
matplot(THETA[, 10:13], type="l", ylab="theta3�̉�A�W��", xlab="�T���v�����O��")


#���U�����U�s��̃v���b�g
matplot(SIGMA[, 1:9], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 10:18], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 19:27], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 28:36], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 37:45], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 46:54], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 55:63], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 64:72], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 73:81], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")