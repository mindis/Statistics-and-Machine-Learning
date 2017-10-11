#####���ւ̂��鑽�����W�b�g���f��#####
library(MASS)
library(mlogit)
library(MCMCpack)
library(bayesm)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�C�ӂ̕��U�����U�s����쐬������֐�####
##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho <- ifelse(abs(rho) < 0.1, 0, rho)
  rho[upper.tri(rho)] <- 0
  
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  X.Sigma <- eigen(Sigma)
  Lambda <- diag(X.Sigma$values)
  P <- X.Sigma$vector
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda)
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
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
hh <- 3000   #�T���v����
choise <- 5   #�I���\��
st <- 5   #��u�����h
k <- 5   #�����ϐ��̐�

##�����ϐ��̔���
#�ʏ퉿�i�̔���
PRICE <- matrix(runif(hh*choise, 0.6, 1), nrow=hh, ncol=choise) - 1   

#�f�B�X�J�E���g���̔���
DISC <- matrix(runif(hh*choise, 0, 0.4), nrow=hh, ncol=choise) 

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

#�J�e�S���[���C�����e�B
ROYL <- matrix(rnorm(hh), nrow=hh, ncol=1)

##�����ϐ����x�N�g����
#�ؕЂ̃x�N�g����
BP <- matrix(as.numeric(diag(choise)), nrow=hh*choise, ncol=choise, byrow=T)[, -choise]

#�u�����h���C�����e�B�̃x�N�g����
Royl_vec <- matrix(0, nrow=hh*choise, ncol=choise)
for(i in 1:hh){
  r <- ((i-1)*choise+1):((i-1)*choise+choise)
  Royl_vec[r, ] <- diag(ROYL[i], choise)
}
Royl_vec <- Royl_vec[, -choise]

#���̑��̕ϐ��̃x�N�g����
Price_vec <- as.numeric(t(PRICE))
Disc_vec <- as.numeric(t(DISC))
Disp_vec <- as.numeric(t(DISP))
Camp_vec <- as.numeric(t(CAMP))

#�f�[�^������
round(X <- data.frame(BP=BP, PRICE=Price_vec, DISC=Disc_vec, DISP=Disp_vec, CAMP=Camp_vec, ROYL=Royl_vec), 2)
XM <- as.matrix(X)

##ID�̐ݒ�
id <- rep(1:hh, rep(choise, hh))
c <- rep(1:choise, hh)
ID <- data.frame(no=1:(hh*choise), id=id, c=c)


####�����ϐ��𔭐�####
##���U�����U�s��̐ݒ�
corM <- corrM(col=choise, lower=-0.55, upper=0.9, eigen_lower=0.01, eigen_upper=0.2)   #���֍s����쐬
Sigma <- covmatrix(col=choise, corM=corM, lower=1, upper=1)   #���U�����U�s��
Cov <- Sigma$covariance

##�p�����[�^�̐ݒ�
beta1 <- -5.5   #���i�̃p�����[�^
beta2 <- 5.7   #�������̃p�����[�^
beta3 <- 2.0   #���ʒ�̃p�����[�^
beta4 <- 1.8   #�L�����y�[���̃p�����[�^
beta5 <- c(1.1, 0.6, -0.5, 0.3)   #�J�e�S���[���C�����e�B�̃p�����[�^
beta0 <- c(0.5, 1.1, 1.4, 2.2)   #�u�����h1�`4�̑��΃x�[�X�̔���
betat <- c(beta0, beta1, beta2, beta3, beta4, beta5)

##���W�b�g�Ɗm���̌v�Z
logit0 <- matrix(XM %*% betat, nrow=hh, ncol=choise, byrow=T)
logit <- logit0 + mvrnorm(hh, rep(0, choise), Cov)
Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=choise)

##�J�e�S���J�����z����I�����𔭐�
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
colMeans(Y); colSums(Y)


####�}���R�t�A�������e�J�����@�ő��֍\���̂��鍬���^���W�X�e�B�b�N��A���f���𐄒�####
##�������W�b�g���f���̑ΐ��ޓx�֐�
LLike <- function(beta, theta, X, Y, hh, choise){
  logit <- matrix(X %*% beta, nrow=hh, ncol=choise, byrow=T) + theta
  
  d <- rowSums(exp(logit))
  LLl <- rowSums(Y * logit) - log(d)
  LL <- sum(LLl)
  LL_val <- list(LL=LL, LLl=LLl)
  return(LL_val)
}

##�����l�̐ݒ�p�̑������W�b�g�̑ΐ��ޓx
loglike <- function(theta, X, Y, hh, choise){
  matrix(XM %*% runif(ncol(XM)), nrow=hh, ncol=choise, byrow=T)
  
  logit <- matrix(XM %*% theta, nrow=hh, ncol=choise, byrow=T)
  
  d <- rowSums(exp(logit))
  LLl <- rowSums(Y * logit) - log(d)
  LL <- sum(LLl)
  return(LL)
}


##MCMC�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- c()   #�ΐ��ޓx�̕ۑ��p

##�����ϐ��𑽎����z��
X.array <- array(0, dim=c(choise, ncol(XM), hh))
for(i in 1:hh){
  X.array[, , i] <- XM[ID$id==i, ]
}
YX.array <- array(0, dim=c(choise, ncol(XM)+1, hh))
id_r <- matrix(1:(hh*choise), nrow=hh, ncol=choise, byrow=T)


##���O���z�̐ݒ�
betas <- rep(0, ncol(XM))
rootBi <- 0.01 * diag(ncol(XM))

nu <- ncol(XM)-2   #�t�E�B�V���[�g���z�̎��R�x
V <- nu*diag(choise)   #�t�E�B�V���[�g���z�̃p�����[�^
Deltabar <- rep(0, ncol(XM))  #��A�W���̕��ς̎��O���z
Adelta <- 0.01 * diag(rep(1, ncol(XM)))   #��A�W���̎��O���z�̕��U

##�T���v�����O���ʂ̕ۑ��p
THETA <- array(0, dim=c(hh, choise, R/keep))
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
SIGMA <- matrix(0, nrow=R/keep, ncol=choise^2)

##�����l�̐ݒ�
#��A�W���̏����l
beta00 <- rep(0, ncol(XM))
res <- optim(beta00, loglike, gr=NULL, XM, Y, hh, choise, method="BFGS", hessian=TRUE, control=list(fnscale=-1))
oldbeta <- res$par
rw <- diag(diag(-solve(res$hessian)))

#���U�����U�s��̏����l
oldcov <- diag(choise)
oldcov0 <- oldcov
cov_inv <- solve(oldcov)

#���p�̕��ύ\���̏����l
theta_mu <- matrix(XM %*% oldbeta, nrow=hh, ncol=choise, byrow=T)
oldtheta <- mvrnorm(hh, rep(0, choise), oldcov)


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){

  ##MH�@�ŉ�A�W��beta���T���v�����O
  betad <- oldbeta
  betan <- betad + 0.25 * mvrnorm(1, rep(0, length(oldbeta)), rw)
  
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew1 <- LLike(beta=betan, theta=oldtheta, X=XM, Y=Y, hh=hh, choise=choise)$LL
  logold1 <- LLike(beta=betad, theta=oldtheta, X=XM, Y=Y, hh=hh, choise=choise)$LL
  logpnew1 <- lndMvn(betan, betas, rootBi)
  logpold1 <- lndMvn(betad, betas, rootBi)
  
  #�T���v�����O���邩�ǂ����̌���
  #MH�T���v�����O
  alpha1 <- min(1, exp(lognew1 + logpnew1 - logold1 - logpold1))
  if(alpha1 == "NAN") alpha1 <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V����beta���̑�
  if(u < alpha1){
    oldbeta <- betan
    logl <- lognew1
    
    #�����łȂ��Ȃ�beta���X�V���Ȃ�
  } else {
    logl <- logold1
  }
  
  ##MH�T���v�����O��theta���T���v�����O
  thetad <- oldtheta
  thetan <- oldtheta + 0.5 * mvrnorm(hh, rep(0, choise), diag(choise))
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew2 <- LLike(beta=oldbeta, theta=thetan, X=XM, Y=Y, hh=hh, choise=choise)$LLl
  logold2 <- LLike(beta=oldbeta, theta=thetad, X=XM, Y=Y, hh=hh, choise=choise)$LLl
  logpnew2 <- apply(thetan, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  logpold2 <- apply(thetad, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  
  
  #�T���v�����O���̑����邩�ǂ���������
  rand <- runif(hh)   #��l�������痐���𔭐�
  LLind.diff <- exp(lognew2 + logpnew2 - logold2 - logpold2)   #�̑𗦂��v�Z
  alpha2 <- ifelse(LLind.diff > 1, 1, LLind.diff)
  alpha2 <- matrix(alpha2, nrow=hh, ncol=choise)
  
  #alpha�Ɋ�Â�beta���̑�
  oldtheta.r <- ifelse(alpha2 > rand, thetan, oldtheta)   #alpha��rand�������Ă�����̑�
  adopt <- sum(oldtheta[, 1]!=oldtheta.r[, 1])/hh   #�̑�
  oldtheta <- oldtheta.r   #�p�����[�^���X�V
  
  ##�K�w���f���̕��U�����U�p�����[�^���T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^���v�Z
  R.error <- oldtheta
  IW.R <- V + t(R.error) %*% R.error
  
  #�t�E�B�V���[�g���z�̎��R�x���v�Z
  Sn <- nu + hh
  
  #�t�E�B�V���[�g���z����Cov���T���v�����O
  Cov_hat <- rwishart(Sn, solve(IW.R))$IW
  
  #���ʐ��̂��߃p�����[�^�ɐ����������
  oldcov <- cov2cor(Cov_hat)
  cov_inv <- solve(oldcov)
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep   
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- oldcov
    
    print(rp)
    print(round(c(logl, res$value), 2))
    print(round(rbind(oldbeta, res$par, betat), 3))
    print(round(cbind(oldcov, Cov), 3))
    print(round(adopt, 3))
  }
}

matplot(BETA[, 7:12], type="l")
matplot(SIGMA[, 1:5], type="l")