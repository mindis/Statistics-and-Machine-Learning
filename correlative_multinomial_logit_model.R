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
hh <- 1500   #�T���v����
choise <- 5   #�I���\��
st <- 5   #��u�����h
k <- 5   #�����ϐ��̐�

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

#�J�e�S���[���C�����e�B
ROYL <- matrix(runif(hh), nrow=hh, ncol=1)

##���U�����U�s��̐ݒ�
corM <- corrM(col=choise, lower=-0.7, upper=0.8)   #���֍s����쐬
Sigma <- covmatrix(col=choise, corM=corM, lower=1, upper=1.5)   #���U�����U�s��
Cov <- Sigma$covariance

##�p�����[�^�̐ݒ�
beta1 <- -6.5   #���i�̃p�����[�^
beta2 <- 6.3   #�������̃p�����[�^
beta3 <- 2.0   #���ʒ�̃p�����[�^
beta4 <- 1.8   #�L�����y�[���̃p�����[�^
beta5 <- c(1.1, 0.6, -0.5, 0.3)   #�J�e�S���[���C�����e�B�̃p�����[�^
beta0 <- c(0.5, 1.1, 1.4, 2.2)   #�u�����h1�`4�̑��΃x�[�X�̔���
betat <- c(beta0, beta1, beta2, beta3, beta4, beta5)

##���p�𔭐������A�I�����ꂽ�u�����h������
#���ϗʐ��K���z���烍�W�b�g�𔭐�
logit.l <- matrix(0, nrow=hh, ncol=st)
for(i in 1:(st-1)){
  logit.l[, i] <- beta0[i] + beta1*PRICE[, i]  + beta2*DISC[, i] + beta3*DISP[, i] + 
                beta4*CAMP[, i] + beta5[i]*ROYL
}
#��ϐ��̃��W�b�g���v�Z
logit.l[, st] <- beta1*PRICE[, st]  + beta2*DISC[, st] + beta3*DISP[, st] + beta4*CAMP[, st]

##���ϗʐ��K�����𔭐�
logit <- logit.l + mvrnorm(n=hh, rep(0, choise), Cov)

##�������������W�b�g����I���u�����h������
#�u�����h�I���m�����v�Z
Pr <- exp(logit)/rowSums(exp(logit))
colMeans(Pr); apply(Pr, 2, summary)

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

round(X <- data.frame(BP=BP, PRICE=PRICE.v, DISC=DISC.v, DISP=DISP.v, CAMP=CAMP.v, ROYL=ROYL.v), 2)   #�f�[�^�̌���
XM <- as.matrix(X)

#ID�̐ݒ�
brand <- rep(1:choise, hh)
id <- rep(1:hh, rep(choise, hh))
ID <- data.frame(id, brand)

##MCMC�A���S���Y���̐ݒ�
R <- 15000
sbeta <- 1.5
keep <- 2
llike <- c()   #�ΐ��ޓx�̕ۑ��p

##���O���z�̐ݒ�
nu <- ncol(X)-2   #�t�E�B�V���[�g���z�̎��R�x
V <- nu*diag(choise)   #�t�E�B�V���[�g���z�̃p�����[�^
Deltabar <- rep(0, ncol(X))  #��A�W���̕��ς̎��O���z
Adelta <- 100 * diag(rep(1, ncol(X)))   #��A�W���̎��O���z�̕��U

##�T���v�����O���ʂ̕ۑ��p
Util <- array(0, dim=c(hh, choise, R/keep))
BETA <- matrix(0, nrow=R/keep, ncol=choise)
SIGMA <- array(0, dim=c(choise, choise, R/keep))

##�����l�̐ݒ�
#��A�W���̏����l
oldbeta <- c(runif(choise-1, 0, 3), -3.0, 3.0, runif(2, 0, 2), runif(choise-1, -2, 3))  

#���U�����U�s��̏����l
corM.f <- corrM(col=choise, lower=-0.6, upper=0.6)   #���֍s����쐬
Sigma.f <- covmatrix(col=choise, corM=corM.f, lower=1, upper=5)   #���U�����U�s��
oldcov <- Sigma.f$covariance


#���p�̕��ύ\���̏����l
b <- oldbeta[(choise+1):length(oldbeta)]

round(cbind(Y, Pr, exp(old.util)/rowSums(exp(old.util)), 
exp(matrix(u, nrow=hh, ncol=choise, byrow=T))/rowSums(exp(matrix(u, nrow=hh, ncol=choise, byrow=T)))), 2)


####�}���R�t�A�������e�J�����@�Ő���####
##�f�[�^�g��@�ő��ϗʐ��K���z������݌��p�𔭐�������
u <- XM %*% oldbeta   #�x�N�g���`���̌��p�̕��ύ\��
U <- matrix(u, nrow=hh, ncol=choise, byrow=T)   #�s��`���̌��p�̕��ύ\��
util.M <- t(apply(U, 1, function(x) mvrnorm(1, x, oldcov)))   #���ϗʐ��K�����Ő��݌��p�𔭐�
U.old <- as.numeric(t(util.M))
I <- diag(hh)

for(rp in 1:50000){
##��A���f���̃M�u�X�T���v�����O��beta��sigma�𐄒�
  u.vec <- as.numeric(t(util.M))   #���݌��p���x�N�g���ɕύX
  
  ##beta�̃M�u�X�T���v�����O
  oldcovi <- solve(oldcov)
  SIGMA.B <- kronecker(I, oldcovi)
  
  #��A�W���̕��ύ\��
  B <- solve(t(XM) %*% SIGMA.B %*% XM) %*% t(XM) %*% SIGMA.B %*% u.vec   #��A�W���̍ŏ���搄���
  XVX <- t(XM) %*% SIGMA.B %*% XM
  BETA.M <- solve(XVX + solve(Adelta)) %*% (XVX %*% B + solve(Adelta) %*% Deltabar)
  
  #��A�W���̕��U�����U�s��
  BETA.SIG <- solve(XVX + solve(Adelta))
  
  #���ϗʐ��K���z�����A�W�����T���v�����O
  oldbeta <- mvrnorm(1, as.numeric(BETA.M), BETA.SIG)
  
  ##sigma�̃M�u�X�T���v�����O
  #�t�E�B�V���[�g���z�̎��R�x���v�Z
  Sn <- nu + hh

  #�t�E�B�V���[�g���z�̃p�����[�^���v�Z
  #���덷���v�Z���Ęa�����
  Vi <- solve(V) 
  EE <- matrix(0, nrow=choise, ncol=choise)
  redi <- u.vec - XM %*% oldbeta
  
  for(i in 1:hh){
    r <- (i-1)*(choise)
    error <- redi[(r+1):(r+choise)] %*% t(redi[(r+1):(r+choise)])
    EE <- EE + error
  }
  R <- solve(EE + Vi)
  
  #�t�E�B�V���[�g�����𔭐�
  oldcov <- rwishart(Sn, R)$IW
 
  u <- XM %*% oldbeta   #�x�N�g���`���̌��p�̕��ύ\��
  U <- matrix(u, nrow=hh, ncol=choise, byrow=T)   #�s��`���̌��p�̕��ύ\��
  util.new <- t(apply(U, 1, function(x) mvrnorm(1, x, oldcov)))   #���ϗʐ��K�����Ő��݌��p�𔭐�
  util.old <- util.M
  
  dnew <- exp(util.new[, 1]) + exp(util.new[, 2]) + exp(util.new[, 3]) + 
          exp(util.new[, 4]) + exp(util.new[, 5])
  dold <- exp(util.old[, 1]) + exp(util.old[, 2]) + exp(util.old[, 3]) + 
          exp(util.old[, 4]) + exp(util.old[, 5])
  
  LLind.new <- Y[, 1]*util.new[, 1] + Y[, 2]*util.new[, 2] + Y[, 3]*util.new[, 3] + 
               Y[, 4]*util.new[, 4] + Y[, 5]*util.new[, 5] - log(dnew)
  LLind.old <- Y[, 1]*util.old[, 1] + Y[, 2]*util.old[, 2] + Y[, 3]*util.old[, 3] + 
               Y[, 4]*util.old[, 4] + Y[, 5]*util.old[, 5] - log(dold)                                                                        

  rand <- matrix(runif(hh), nrow=hh, ncol=choise)
  LLind.diff <- exp(LLind.new - LLind.old)
  alpha <- matrix(ifelse(LLind.diff > 1, 1, LLind.diff), nrow=hh, ncol=choise)

  util.M <- ifelse(alpha > rand, util.M <- util.new, util.M <- util.old)
  logl <- ifelse(alpha > rand, logl <- LLind.new, logl <- LLind.old)
  LL <- sum(logl)
  
  print(LL)
  print(rp)
  print(rbind(round(oldbeta, 3), round(betat, 3)))
  print(cbind(round(cov2cor(oldcov), 3), round(cov2cor(Cov), 3)))
}