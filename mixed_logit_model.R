####���ւ̂��鑽�����W�b�g���f��#####
library(MASS)
library(mlogit)
library(MCMCpack)
library(bayesm)
library(caret)
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
hh <- 1000   #�T���v����
choise <- 5   #�I���\��
st <- 5   #��u�����h
k <- 5   #�����ϐ��̐�

##�����ϐ��̔���
#�ʏ퉿�i�̔���
PRICE <- matrix(runif(hh*choise, 0.6, 1), nrow=hh, ncol=choise)   

#�f�B�X�J�E���g���̔���
DISC <- matrix(runif(hh*choise, 0, 0.5), nrow=hh, ncol=choise)

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
ROYL <- matrix(runif(hh, 0, 1), nrow=hh, ncol=1)

##���U�����U�s��̐ݒ�
MU <- matrix(c(1, 1, 0, 1, 1), nrow=hh, ncol=choise, byrow=T)
Cov <- matrix(c(3.5, 0, 0, 1), nrow=2, ncol=2)
RMV <- mvrnorm(hh, rep(0, 2), Cov)   #���U����
MU <- cbind(RMV[, 1], RMV[, 1], 0, RMV[, 2], RMV[, 2])


##�p�����[�^�̐ݒ�
beta1 <- -5.8   #���i�̃p�����[�^
beta2 <- 5.5   #�������̃p�����[�^
beta3 <- 2.0   #���ʒ�̃p�����[�^
beta4 <- 1.8   #�L�����y�[���̃p�����[�^
betat <- c(beta1, beta2, beta3, beta4)

##�K�w���f���̐ݒ�
b1 <- c(1.1, 0.6, -0.7, -0.3)   #�J�e�S���[���C�����e�B�̃p�����[�^
b0 <- c(0.5, 0.8, 1.2, 2.0)   #�u�����h1�`4�̑��΃x�[�X�̔���
beta0 <- matrix(b0, nrow=hh, ncol=choise-1, byrow=T)

##���p�𔭐������A�I�����ꂽ�u�����h������
#���ϗʐ��K���z���烍�W�b�g�𔭐�
logit <- matrix(0, nrow=hh, ncol=st)
for(i in 1:(st-1)){
  logit[, i] <- beta0[, i] + beta1*PRICE[, i]  + beta2*DISC[, i] + beta3*DISP[, i] + beta4*CAMP[, i] + MU[, i]
}
#��ϐ��̃��W�b�g���v�Z
logit[, st] <- beta1*PRICE[, st]  + beta2*DISC[, st] + beta3*DISP[, st] + beta4*CAMP[, st] + MU[, st]

##�������������W�b�g����I���u�����h������
#�u�����h�I���m�����v�Z
Pr <- exp(logit)/rowSums(exp(logit))
colMeans(Pr); apply(Pr, 2, summary)

round(cbind(Pr, Pr1), 3)

#�I���u�����h�𔭐�
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
colMeans(Y); apply(Y, 2, table)

round(cbind(Y %*% 1:choise, Pr), 3)

####�}���R�t�A�������e�J�����@�ő��֍\���̂��鍬���^���W�X�e�B�b�N��A���f���𐄒�####
##��A���f���𐄒肷�邽�߂ɐ����ϐ����x�N�g���`���ɕύX�ݒ�
#ID�̐ݒ�
ID <- rep(1:hh, rep(choise, hh))

#�ؕЂ̐ݒ�
p <- c(1, rep(0, choise))
bp <- matrix(p, nrow=hh*(choise+1), ncol=choise, byrow=T)
BP <- subset(bp, rowSums(bp) > 0)
BP <- BP[, -st] 

#�J�e�S�����C�����e�B�̐ݒ�
ROYL.v <- matrix(0, nrow=hh*choise, ncol=choise)
for(i in 1:hh){
  ROYL.v[ID==i, ] <- diag(c(rep(ROYL[i, ], choise-1), 0))
}
ROYL.v <- ROYL.v[, -st]

#�����ϐ��̐ݒ�
PRICE.v <- as.numeric(t(PRICE))
DISC.v <- as.numeric(t(DISC))
DISP.v <- as.numeric(t(DISP))
CAMP.v <- as.numeric(t(CAMP))

round(X <- data.frame(BP, PRICE=PRICE.v, DISC=DISC.v, DISP=DISP.v, CAMP=CAMP.v), 2)   #�f�[�^�̌���
XM <- as.matrix(X)

#ID�̐ݒ�
brand <- rep(1:choise, hh)
id <- rep(1:hh, rep(choise, hh))
ID <- data.frame(id, brand)

#Z�̐ݒ�
z <- cbind(c(1, 1, 0, 0, 0), c(rep(0, choise)), c(0, 0, 0, 1, 1))
Z <- kronecker(diag(hh), z)
index.z <- subset(1:ncol(Z), colSums(Z)==0)
Z <- Z[, -index.z]

##MCMC�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 2
llike <- c()   #�ΐ��ޓx�̕ۑ��p

##���O���z�̐ݒ�
#�Œ���ʂ̎��O���z
betas.fix <- rep(0, ncol(XM))  #��A�W���̕��ς̎��O���z
sigma.fix <- diag(rep(0.01, ncol(XM)))   #��A�W���̎��O���z�̕��U

#�ϗʌ��ʂ̎��O���z
Deltabar <- rep(0, hh*choise)
Adelta <- 0.01*diag(2)
nu <- 2   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(rep(1, 2))
beta.random <- matrix(0, nrow=hh, ncol=2)   #�ϗʌ��ʂ̎��O���z�̕��ς�0�ɌŒ�


##�T���v�����O���ʂ̕ۑ��p
Util <- array(0, dim=c(hh, choise-1, R/keep))
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
BETA0 <- matrix(0, R/keep, ncol=2)
SIGMA <- matrix(0, nrow=R/keep, ncol=2^2)

##�����l�̐ݒ�
#��A�W���̏����l
oldbeta.f <- c(c(runif(choise-1, 0, 3)), -3.0, 3.0, runif(2, 0, 2))  

#�ϗʌ��ʂ̏����l
cov.random <- diag(runif(2, 0.1, 1))
oldbeta.r <- matrix(mvrnorm(hh, rep(0, 2), cov.random), nrow=hh, ncol=2, byrow=T)   #�ϗʌ��ʂ̏����l

#�K�w���f���̏����l
beta.random <- matrix(0, nrow=hh, ncol=2)

####�}���R�t�A�������e�J�����@�Ő���####
##mixed logit���f���̑ΐ��ޓx
LLike <- function(beta, b, X, Z, Y, hh, choise){
  d <- rowSums(exp(matrix(X %*% beta + Z %*% b, nrow=hh, ncol=choise, byrow=T)))
  LLl <- rowSums(Y * matrix(X %*% beta + Z %*% b, nrow=hh, ncol=choise, byrow=T)) - log(d)
  LL <- sum(LLl)
  LL.val<- list(LLl=LLl, LL=LL)
  return(LL.val)
}


for(rp in 1:R){
  ##MH�T���v�����O�ŌŒ����beta�̃T���v�����O
  oldbeta.rv <- as.numeric(t(oldbeta.r))
  betad.f <- oldbeta.f
  betan.f <- betad.f + rnorm(length(betad.f), 0, 0.05)   #�����_���E�H�[�N�T���v�����O
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew.f <- LLike(beta=betan.f, b=oldbeta.rv, X=XM, Z=Z, Y=Y, hh=hh, choise=choise)$LL
  logold.f <- LLike(beta=betad.f, b=oldbeta.rv, X=XM, Z=Z, Y=Y, hh=hh, choise=choise)$LL
  logpnew.f <- lndMvn(betan.f, betas.fix, sigma.fix)
  logpold.f <- lndMvn(betad.f, betas.fix, sigma.fix)
  
  #MH�T���v�����O
  alpha.f <- min(1, exp(lognew.f + logpnew.f - logold.f - logpold.f))
  if(alpha.f == "NAN") alpha.f <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V�����Œ����beta���̑�
  if(u < alpha.f){
    oldbeta.f <- betan.f
    logl.f <- lognew.f
    
    #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
  } else {
    logl.f <- logold.f
  }
  
  
  ##MH�T���v�����O�Ōl�ʂɕϗʌ���beta���T���v�����O
  betad.random <- oldbeta.r 
  rw <- t(1.25 * chol(cov.random) %*% t(matrix(rnorm(hh*(2)), nrow=hh, ncol=2)))
  
  betan.random <- betad.random + rw
  betad.r <- as.numeric(t(betad.random))
  betan.r <- as.numeric(t(betan.random))
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew.r <- LLike(beta=oldbeta.f, b=betan.r, X=XM, Z=Z, Y=Y, hh=hh, choise=choise)$LLl
  logold.r <- LLike(beta=oldbeta.f, b=betad.r, X=XM, Z=Z, Y=Y, hh=hh, choise=choise)$LLl
  logpnew.r <- apply((betan.random - beta.random), 1, function(x) -0.5 * x %*% solve(cov.random) %*% x)
  logpold.r <- apply((betad.random - beta.random), 1, function(x) -0.5 * x %*% solve(cov.random) %*% x)
  
  #MH�T���v�����O
  rand <- matrix(runif(hh), nrow=hh, ncol=2)
  LLind.diff <- exp(lognew.r + logpnew.r - logold.r - logpold.r)   #���p�����v�Z
  alpha <- matrix(ifelse(LLind.diff > 1, 1, LLind.diff), nrow=hh, ncol=2)      
  
  oldbeta.r <- ifelse(alpha > rand, oldbeta.r <- betan.random, oldbeta.r <- betad.random)   #alpha��rand�������Ă�����̑�
  logl <- ifelse(alpha[, 1] > rand[, 1], logl <- lognew.r, logl <- logold.r)
  
 
  ##���K���z����beta0���T���v�����O
  #beta0 <- colMeans(oldbeta.r)
  #ohm <- cov.random/hh
  #beta0.mv <-  mvrnorm(1, beta0, ohm)
  #beta0.mv + t(chol(cov.random/hh))*rnorm(choise-1)
  #beta.random <- matrix(beta0.mv, nrow=hh, ncol=choise-1, byrow=T)
  
  ##�t�E�B�V���[�g���z����sigma���T���v�����O
  V <- var(oldbeta.r)
  VK <- 2 * diag(2) + hh * V
  nu1 <- hh + nu - 1 
  
  cov.random <- rwishart(nu1, solve(VK))$IW   #�t�E�B�V���[�g���z���番�U�����U�s��𔭐�

  
  ##�M�u�X�T���v�����O��Delta���T���v�����O
  #M <- matrix(c(1, 0), hh, 2, byrow=T)   #���z�I��0�̐����ϐ����쐬
  #DeltaM <- matrix(Deltabar, 2, 2, byrow=T)   #���z�I��0�̉�A�W���̎��O���z���쐬
  
  #out <- rmultireg(oldbeta.r, M, DeltaM, Adelta, nu, V)   #���ϗʉ�A���f���̃M�u�X�T���v���[
  #cov.random <- out$Sigma
  
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta.f
    BETA0[mkeep, ] <- beta.random[1, ]
    SIGMA[mkeep, ] <- as.numeric(cov.random)
    Util[, , mkeep] <- oldbeta.r
    
    print(sum(logl))
    print(rp)
    print(round(mean(alpha), 3)); print(round(alpha.f, 3))
    print(round(rbind(oldbeta.f, c(b0, betat)), 3))
    print(round(cbind(cov.random, Cov), 3))
  }
}

burnin <- 5000

matplot(BETA[, 1:4], type="l")
matplot(BETA[, 5:8], type="l")
matplot(SIGMA[, c(1, 4)], type="l")

colMeans(SIGMA[, c(1, 4)]) - pi^2/6


