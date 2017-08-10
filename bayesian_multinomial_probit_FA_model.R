#####�����v���r�b�g���f��#####
library(MASS)
library(bayesm)
library(pyhch)
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
hh <- 2500   #�v���C���[��
choise <- 10   #�I���\��
st <- choise   #��u�����h
k <- 5   #��A�W���̐�

##ID�̐ݒ�
id <- rep(1:hh, rep(choise-1, hh))
c <- rep(1:(choise-1), hh)
ID <- data.frame(no=1:(hh*(choise-1)), id=id, c=c)
id_r <- matrix(1:(hh*(choise-1)), nrow=hh, ncol=choise-1, byrow=T)


##�����ϐ��̔���
#�ʏ퉿�i�̔���
PRICE <- matrix(runif(hh*choise, 0.7, 1), nrow=hh, ncol=choise)   

#�f�B�X�J�E���g���̔���
DISC <- matrix(runif(hh*choise, 0, 0.3), nrow=hh, ncol=choise)

#���ʒ�̔���
DISP <- matrix(0, nrow=hh, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hh, 1, r)
}

#���ʃL�����y�[���̔���
CAMP <- matrix(0, nrow=hh, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hh, 1, r)
}

##���U�����U�s��̐ݒ�
corM <- corrM(col=choise-1, lower=-0.55, upper=0.75)   #���֍s����쐬
Sigma <- covmatrix(col=choise-1, corM=corM, lower=1, upper=1)   #���U�����U�s��
Cov <- Sigma$covariance

##�p�����[�^�̐ݒ�
beta1 <- -6.3   #���i�̃p�����[�^
beta2 <- 7.2   #�������̃p�����[�^
beta3 <- 2.0   #���ʒ�̃p�����[�^
beta4 <- 1.8   #�L�����y�[���̃p�����[�^
beta0 <- runif(choise-1, -0.4, 2.1)   #�u�����h1�`4�̑��΃x�[�X�̔���
betat <- c(beta0, beta1, beta2, beta3, beta4)

#��u�����h�Ƃ̑��ΐ����ϐ�
PRICE.r <- PRICE[, -5] - PRICE[, 5]
DISC.r <- DISC[, -5] - DISC[, 5]
DISP.r <- DISP[, -5] - DISP[, 5]
CAMP.r <- CAMP[, -5] - CAMP[, 5]

##��A���f���𐄒肷�邽�߂ɐ����ϐ����x�N�g���`���ɕύX�ݒ�
#�ؕЂ̐ݒ�
p <- c(1, rep(0, choise-1))
bp <- matrix(p, nrow=hh*choise, ncol=choise-1, byrow=T)
BP <- subset(bp, rowSums(bp) > 0)


#�����ϐ��̐ݒ�
PRICE.v <- as.numeric(t(PRICE.r))
DISC.v <- as.numeric(t(DISC.r))
DISP.v <- as.numeric(t(DISP.r))
CAMP.v <- as.numeric(t(CAMP.r))

X <- data.frame(BP=BP, PRICE.v, DISC.v, DISP.v, CAMP.v)   #�f�[�^�̌���
XM <- as.matrix(X)

##���Ό��p�𔭐������A�I�����ꂽ�u�����h������
U.mean <- matrix(XM %*% betat, nrow=hh, ncol=choise-1, byrow=T)   #���Ό��p�̕��ύ\��
U <- U.mean + mvrnorm(hh, rep(0, choise-1), Cov)   #�덷�\�������������p

#���p�ő剻�����Ɋ�Â��I���u�����h������
Y <- apply(U, 1, function(x) ifelse(max(x) < 0, choise, which.max(x)))


#�w����0�A1�s��ɕύX
BUY <- matrix(0, hh, choise)
for(i in 1:hh){
  BUY[i, Y[i]] <- 1
}

table(Y)   #�I���u�����h�̏W�v
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
factors <- 3
llike <- array(0, dim=c(R/keep))   #�ΐ��ޓx�̕ۑ��p


##�����ϐ��𑽎����z��
X.array <- array(0, dim=c(choise-1, ncol(X), hh))
for(i in 1:hh){
  X.array[, , i] <- XM[ID[, 2]==i, ]
}
YX.array <- array(0, dim=c(choise-1, ncol(X)+1, hh))


#����v���Z�X�̊i�[�z��
UM <- matrix(0, nrow=hh, ncol=choise-1)
util.M <- matrix(0, nrow=hh, ncol=choise-1)   

##���O���z�̐ݒ�
nu <- choise   #�t�E�B�V���[�g���z�̎��R�x
V <- solve((1/10)*diag(choise-1))    #�t�E�B�V���[�g���z�̃p�����[�^
Deltabar <- rep(0, ncol(X))  #��A�W���̕��ς̎��O���z
Adelta <- solve(100 * diag(rep(1, ncol(X))))   #��A�W���̎��O���z�̕��U

inv.facov <- solve(diag(100, factors))
alpha_d <- 1
beta_d <- 100


##�T���v�����O���ʂ̕ۑ��p�z��
Util <- array(0, dim=c(hh, choise-1, R/keep))
BETA <- matrix(0, nrow=R/keep, length(beta0)+k-1)
SIGMA <- matrix(0, nrow=R/keep, ncol=(choise-1)^2)
FA.A <- matrix(0, nrow=R/keep, ncol=(choise-1)*factors)
FA.D <- matrix(0, nrow=R/keep, ncol=choise-1)
FA.F <- array(0, dim=c(hh, factors, R/keep))


##�����l�̐ݒ�
#��A�W���̏����l
oldbeta <- c(runif(choise-1, 0, 3), -3.0, 3.0, runif(2, 0, 2))   

#���U�����U�s��̏����l
corM.f <- corrM(col=choise-1, lower=0, upper=0)   #���֍s����쐬
Sigma.f <- covmatrix(col=choise-1, corM=corM.f, lower=1, upper=1)   #���U�����U�s��
oldcov <- Sigma.f$covariance

#���p�̕��ύ\���̏����l
old.utilm <- matrix(XM %*% oldbeta, nrow=hh, ncol=choise-1, byrow=T)

#���p�̏����l
old.util <- old.utilm + mvrnorm(nrow(old.utilm), rep(0, choise-1), oldcov)

#���q���חʂƓƎ����q�̏����l
A <- matrix(runif((choise-1)*factors, -1, 1), nrow=choise-1, ncol=factors)
D <- diag(runif(choise-1, 0, 0.5))


####�}���R�t�A�������e�J�����@�ő����v���r�b�g���f���𐄒�####
for(rp in 1:R){
  
  ##�I�����ʂƐ����I�Ȑ��݌��p�𔭐�������
  #�����t�����Ғl�Ə����t�����U���v�Z
  S <- rep(0, choise-1)
  
  for(j in 1:(choise-1)){
    MVR <- cdMVN(mu=old.utilm, Cov=oldcov, dependent=j, U=old.util)   #�����t�����z���v�Z
    UM[, j] <- MVR$CDmu   #�����t�����Ғl�����o��
    S[j] <- sqrt(MVR$CDvar)    #�����t�����U�����o��
    
    #���ݕϐ��𔭐�������
    #�ؒf�̈�̐ݒ�
    max.u <- apply(cbind(old.util[, -j], 0), 1, max)
    max.u <- ifelse(Y==choise, 0, max.u)
    
    #�ؒf���K���z�����ݕϐ��𔭐�
    old.util[, j] <- ifelse(Y==j, rtnorm(mu=UM[, j], sigma=S[j], a=max.u, b=100), 
                            rtnorm(mu=UM[, j], sigma=S[j], a=-100, b=max.u))
    old.util[, j] <- ifelse(is.infinite(old.util[, j]), ifelse(Y==j, max.u + runif(1), max.u - runif(1)), old.util[, j])
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
  R.error <- matrix(util.v - XM %*% oldbeta, nrow=hh, ncol=choise-1, byrow=T)
  IW.R <- V + matrix(rowSums(apply(R.error, 1, function(x) x %*% t(x))), nrow=choise-1, ncol=choise-1)
  
  #�t�E�B�V���[�g���z�̎��R�x���v�Z
  Sn <- nu + hh
  
  #�t�E�B�V���[�g���z����Cov���T���v�����O
  Cov_hat <- rwishart(Sn, solve(IW.R))$IW
  oldcov <- cov2cor(Cov_hat)
  
  ##���݌��p�ƃp�����[�^���X�V
  #���݌��p�Ɛ��݌��p�̕��ς��X�V
  old.utilm <- matrix(XM %*% oldbeta, nrow=hh, ncol=choise-1, byrow=T)
  Z <- old.util - old.utilm 
  Z <- scale(Z)
  
  ##���݌��p�̌덷��������q���̓��f���𐄒�
  #���ϗʐ��K���z������ݕϐ�f(���ʈ��q)���T���v�����O
  ADA <- t(A) %*% solve(A %*% t(A) + D)
  F_mean <- Z %*% t(ADA)   #���ʈ��q�̕���
  F_var <- diag(factors) - ADA %*% A    #���ʈ��q�̕��U�����U�s��
  Fi <- t(apply(F_mean, 1, function(x) mvrnorm(1, x, F_var)))   #���ϗʐ��K���z���狤�ʈ��q���T���v�����O 
  
  
  #�K���}���z����Ǝ����qd���T���v�����O
  Z.error <- Z - Fi %*% t(A)
  Zv.R <- matrix(rowSums(apply(Z.error, 1, function(x) x %*% t(x))), nrow=choise-1, ncol=choise-1)
  
  gamma_alpha <- (hh + alpha_d)/2   #alpha���v�Z
  gamma_beta <- (diag(Zv.R) + beta_d)/2   #beta���v�Z
  D <- diag(rgamma(length(gamma_beta), gamma_alpha, gamma_beta))   #�K���}���z����Ǝ����q���T���v�����O
  
  #���ϗʐ��K���z������q���ח�A���T���v�����O
  FF <- t(Fi) %*% Fi
  FZ <- t(Fi) %*% Z
  d_sigma <- list()
  
  for(i in 1:(choise-1)){
    d_sigma[[i]]  <- 1/diag(D)[i] * inv.facov
    A_mu <- solve(d_sigma[[i]] + FF) %*% FZ[, i]
    A_cov <- solve(inv.facov + diag(D)[i]*FF) 
    A[i, ] <- mvrnorm(1, A_mu, A_cov)
  }
  A[1, 1] <- 0
  A[2, 2] <- 0
  A[3, 3] <- 0
  
  ##�T���v�����O���ʂ�ۑ�
  if(rp%%keep==0){
    print(rp)

    mkeep <- rp/keep
    Util[, , mkeep] <- old.util
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    FA.A[mkeep, ] <- as.numeric(A)
    FA.D[mkeep, ] <- diag(D)
    FA.F[, , mkeep] <- Fi
    print(round(cbind(oldcov, Cov), 2))
    print(round(rbind(oldbeta, betat), 2))
    print(round(t(A), 2))
  }
}


####�֐��Ő���####
Data1 <- list(p=choise, y=Y, X=XM)
Mcmc1 <- list(R=10000, keep=2)

#�����v���r�b�g���f���𐄒�
out <- rmnpGibbs(Data=Data1,Mcmc=Mcmc1)
BETA.out <- out$betadraw
SIGMA.out <- out$sigmadraw

####���茋�ʂ̗v��ƓK���x�̊m�F####
burnin <- 10000/keep   #�o�[���C������

##�T���v�����O���ʂ�����
#��A�W���̃v���b�g
matplot(BETA[, 1:4], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 5:8], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")

#���U�����U�s��̉���
matplot(SIGMA[, 1:4], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 5:9], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 10:13], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 14:18], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 19:22], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 23:27], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 28:31], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 32:36], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 37:40], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 41:45], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 46:49], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 50:54], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 55:58], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 59:63], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 64:67], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 68:72], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 73:76], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 77:81], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")

#���q���חʂ̉���
matplot(FA.A[, 1:4], type="l", main="���q���חʂ̃T���v�����O����", ylab="�p�����[�^����l")
matplot(FA.A[, 5:9], type="l", main="���q���חʂ̃T���v�����O����", ylab="�p�����[�^����l")
matplot(FA.A[, 10:13], type="l", main="���q���חʂ̃T���v�����O����", ylab="�p�����[�^����l")
matplot(FA.A[, 14:18], type="l", main="���q���חʂ̃T���v�����O����", ylab="�p�����[�^����l")
matplot(FA.A[, 19:22], type="l", main="���q���חʂ̃T���v�����O����", ylab="�p�����[�^����l")
matplot(FA.A[, 23:27], type="l", main="���q���חʂ̃T���v�����O����", ylab="�p�����[�^����l")


##����l�̎��㕽�ς̔�r
#beta�̗v�񓝌v��
round(colMeans(BETA.out[burnin:nrow(BETA.out), ] / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(�֐�����)�̎��㕽��
round(colMeans(BETA[burnin:nrow(BETA), ]), 3)   #beta�̎��㕽��
round(betat, 3)   #�^�̒l
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5�����ʓ_
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95�����ʓ_
round(apply(BETA[burnin:nrow(BETA), ], 2, sd), 2)   #����W���΍�

#sigma�̗v�񓝌v��
round(colMeans(SIGMA.out[burnin:nrow(SIGMA.out), ]  / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(�֐�����)�̎��㕽��
round(colMeans(SIGMA[burnin:nrow(SIGMA), ]), 3)   #beta�̎��㕽��
round(as.numeric(Cov), 3)   #�^�̒l
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.05)), 2)   #5�����ʓ_
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.95)), 2)   #95�����ʓ_
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, sd), 2) #����W���΍�