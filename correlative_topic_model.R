#####���փg�s�b�N���f��#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
detach("package:gtools", unload=TRUE)
library(bayesm)
library(ExtDist)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(842573)

####���ϗʐ��K���z�̗����𔭐�������֐����`####
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
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
#set.seed(423943)
#�����f�[�^�̐ݒ�
k <- 8   #�g�s�b�N��
d <- 2000   #������
v <- 300   #��b��
w <- rpois(d, rgamma(d, 160, 1.0))   #1����������̒P�ꐔ

#ID�̐ݒ�
word_id <- rep(1:d, w)

##�p�����[�^�̐ݒ�
alpha1 <- rep(0.25, v)   #�P��̃f�B���N�����O���z�̃p�����[�^

#���ϗʐ��K���z�̃p�����[�^��ݒ�
alpha0 <- rep(0, k-1)   #�����̂̃p�����[�^
tau0 <- corrM(k-1, -0.6, 0.9, 0.01, 0.2)
sigma0 <- covmatrix(k-1, tau0, 2.5, 3.0)$covariance

#���ϗʐ��K���z����̗����𑽍����W�b�g�ϊ����ĕ����g�s�b�N��ݒ�
mv <- cbind(mvrnorm(d, alpha0, sigma0), 0)
theta0 <- exp(mv) / rowSums(exp(mv))

#�f�B���N�������̔���
phi0 <- extraDistr::rdirichlet(k, alpha1)   #�P��̃g�s�b�N�����f�B���N�����z���甭��

##�������z����g�s�b�N����ђP��f�[�^�𔭐�
WX <- matrix(0, nrow=d, ncol=v)
Z1 <- list()

#�������ƂɃg�s�b�N�ƒP��𒀎�����
for(i in 1:d){
  print(i)

  #�����̃g�s�b�N���z�𔭐�
  z1 <- t(rmultinom(w[i], 1, theta0[i, ]))   #�����̃g�s�b�N���z�𔭐�
  
  #�����̃g�s�b�N���z����P��𔭐�
  zn <- z1 %*% c(1:k)   #0,1�𐔒l�ɒu��������
  zdn <- cbind(zn, z1)   #apply�֐��Ŏg����悤�ɍs��ɂ��Ă���
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi0[x[1], ])))   #�����̃g�s�b�N����P��𐶐�
  wdn <- colSums(wn)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  WX[i, ] <- wdn  
  
  #�����������g�s�b�N���i�[
  Z1[[i]] <- zdn[, 1]
}
storage.mode(WX) <- "integer"   #�f�[�^�s��𐮐��^�s��ɕύX


####�g�s�b�N���f���̂��߂̃f�[�^�Ɗ֐��̏���####
##���ꂼ��̕������̒P��̏o������ѕ⏕���̏o�����x�N�g���ɕ��ׂ�
##�f�[�^����pID���쐬
ID_list <- list()
wd_list <- list()

#���l���Ƃɋ��lID����ђP��ID���쐬
for(i in 1:nrow(WX)){
  print(i)
  
  #�P���ID�x�N�g�����쐬
  ID_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
  num2 <- subset(num1, num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
}

#���X�g���x�N�g���ɕϊ�
ID_d <- unlist(ID_list)
wd <- unlist(wd_list)

##�C���f�b�N�X���쐬
doc_list <- list()
word_list <- list()
for(i in 1:length(unique(ID_d))) {doc_list[[i]] <- subset(1:length(ID_d), ID_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- subset(1:length(wd), wd==i)}
gc(); gc()

x %*% etad[4, ]
t(x %*% t(etad))[4, ]

####�}���R�t�A�������e�J�����@�őΉ��g�s�b�N���f���𐄒�####
##�������W�b�g���f���̑ΐ��ޓx�֐�
loglike <- function(beta, y, X, N, select){
  
  #���W�b�g�Ɗm���̌v�Z
  logit <- t(X %*% t(beta))
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=N, ncol=select)
  
  #�ΐ��ޓx���`
  LLi <- rowSums(y * log(Pr))
  LL <- sum(LLi)
  val <- list(LLi=LLi, LL=LL)
  return(val)
}

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

##�A���S���Y���̐ݒ�
R <- 10000   #�T���v�����O��
keep <- 2   #2���1��̊����ŃT���v�����O���ʂ��i�[
iter <- 0

##���O���z�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z
alpha01 <- rep(0, k-1)
nu <- k+1
V <- nu * diag(k-1)
alpha02 <- rep(0.5, v)
beta0m <- matrix(1, nrow=v, ncol=k)

##�p�����[�^�̏����l
#���ϗʐ��K���z����g�s�b�N���z�𔭐�
oldalpha <- rep(0, k-1)
oldcov <- diag(k-1)
cov_inv <- solve(oldcov)
oldeta <- cbind(mvrnorm(d, oldalpha, oldcov), 1)
theta <- exp(mv) / rowSums(exp(mv))

#�f�B�N�������z����P�ꕪ�z�𔭐�
phi.ini <- runif(v, 0.5, 1)
phi <- extraDistr::rdirichlet(k, phi.ini)   #�P��g�s�b�N�̃p�����[�^�̏����l


##�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(d, k, R/keep))
SIGMA <- array(0, dim=c(k-1, k-1, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
Z_SEG <- matrix(0, nrow=length(wd), ncol=k)
storage.mode(Z_SEG) <- "integer"
gc(); gc()


##MCMC����p�z��
wsum0 <- matrix(0, nrow=d, ncol=k)
vf0 <- matrix(0, nrow=v, ncol=k)
x <- diag(k)[, -k]
lognew <- rep(0, d)
logold <- rep(0, d)
logpnew <- rep(0, d)
logpold <- rep(0, d)


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){

  ##�P�ꂲ�ƂɃg�s�b�N���T���v�����O
  #�P�ꂲ�ƂɃg�s�b�N�̏o�������v�Z
  word_rate <- burden_fr(theta, phi, wd, w, k)$Br
  
  #�������z����P��g�s�b�N���T���v�����O
  Zi1 <- rmnom(nrow(word_rate), 1, word_rate)
  word_z <- as.numeric(Zi1 %*% 1:k)
  
  ##���g���|���X�w�C�X�e�B���O�@�ŒP��g�s�b�N�̃p�����[�^���X�V
  #�������W�b�g���f���̉����ϐ��𐶐�
  for(i in 1:d){wsum0[i, ] <- colSums(Zi1[doc_list[[i]], ])}
  
  #���ւ̂��鑽�����W�b�g���f���̃p�����[�^�̃T���v�����O
  #�V�����p�����[�^���T���v�����O
  etad <- oldeta[, -k]
  etan <- etad + matrix(rnorm(d*(k-1), 0, 0.1), nrow=d, ncol=k-1)
  
  #���O���z�̌덷���v�Z
  er_new <- etan - matrix(0, nrow=d, ncol=k-1, byrow=T)
  er_old <- etad - matrix(0, nrow=d, ncol=k-1, byrow=T)

  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew <- loglike(etan, wsum0, x, d, k)$LLi
  logold <- loglike(etad, wsum0, x, d, k)$LLi
  logpnew <- apply(er_new, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  logpold <- apply(er_old, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  
  #���g���|���X�w�C�X�e�B���O�@�Ńp�����[�^�̍̑�������
  rand <- runif(d)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew + logpnew - logold - logpold)   #�̑𗦂��v�Z
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=d, ncol=k-1)
  oldeta <- flag*etan + (1-flag)*etad   #alpha��rand�������Ă�����̑�
  
  #�������W�b�g���f������theta���X�V
  oldeta0 <- cbind(oldeta, 0)
  theta <- exp(oldeta0) / rowSums(exp(oldeta0))
  
  ##�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
  V_par <- d + nu
  R_par <- solve(V) + t(oldeta) %*% oldeta
  oldcov <- rwishart(V_par, solve(R_par))$IW
  cov_inv <- solve(oldcov)

  ##�f�B�N�������z����phi���T���v�����O
  for(i in 1:v){vf0[i, ] <- colSums(Zi1[word_list[[i]], ])}
  vf <- vf0 + beta0m
  phi <- extraDistr::rdirichlet(k, t(vf))
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep1 <- rp/keep
    THETA[, , mkeep1] <- theta
    PHI[, , mkeep1] <- phi
    SIGMA[, , mkeep1] <- cov2cor(oldcov)
    
    #�g�s�b�N�����̓T���v�����O���Ԃ̔����𒴂�����i�[����
    if(rp >= R/2){
      Z_SEG <- Z_SEG + Zi1
    }
  
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(sum(lognew))
    print(round(cbind(cov2cor(oldcov), cov2cor(sigma0)), 3))
    print(round(rbind(theta[1:5, ], theta0[1:5, ]), 3))
  }
}

####�T���v�����O���ʂ̉����Ɨv��####
burnin <- 2000   #�o�[���C������

##�T���v�����O���ʂ̉���
#�����̃g�s�b�N���z�̃T���v�����O����
matplot(t(THETA[1, , ]), type="l", ylab="�p�����[�^", main="����1�̃g�s�b�N���z�̃T���v�����O����")
matplot(t(THETA[2, , ]), type="l", ylab="�p�����[�^", main="����2�̃g�s�b�N���z�̃T���v�����O����")
matplot(t(THETA[3, , ]), type="l", ylab="�p�����[�^", main="����3�̃g�s�b�N���z�̃T���v�����O����")
matplot(t(THETA[4, , ]), type="l", ylab="�p�����[�^", main="����4�̃g�s�b�N���z�̃T���v�����O����")

#���U�����U�s��̃T���v�����O����
matplot(t(SIGMA[1,  , ]), type="l", ylab="�p�����[�^")
matplot(t(SIGMA[2,  , ]), type="l", ylab="�p�����[�^")
matplot(t(SIGMA[3,  , ]), type="l", ylab="�p�����[�^")
matplot(t(SIGMA[4,  , ]), type="l", ylab="�p�����[�^")
matplot(t(SIGMA[5,  , ]), type="l", ylab="�p�����[�^")
matplot(t(SIGMA[6,  , ]), type="l", ylab="�p�����[�^")
matplot(t(SIGMA[7,  , ]), type="l", ylab="�p�����[�^")

#�P��̏o���m���̃T���v�����O����
matplot(t(PHI[1, 1:10, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N1�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[2, 11:20, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N2�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[3, 21:30, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N3�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[4, 31:40, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N4�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[5, 41:50, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N5�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[6, 51:60, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N6�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[7, 61:70, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N7�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[8, 71:80, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N8�̒P��̏o�����̃T���v�����O����")


##�T���v�����O���ʂ̗v�񐄒��
#�g�s�b�N���z�̎��㐄���
topic_mu <- apply(THETA[, , burnin:(R/keep)], c(1, 2), mean)   #�g�s�b�N���z�̎��㕽��
round(cbind(topic_mu, theta0), 3)
round(topic_sd <- apply(THETA[, , burnin:(R/keep)], c(1, 2), sd), 3)   #�g�s�b�N���z�̎���W���΍�

#�P��o���m���̎��㐄���
word_mu <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #�P��̏o�����̎��㕽��
round(rbind(word_mu, phi0)[, 1:50], 3)

#���U�����U�s��̎��㐄���
sigma_mu <- apply(SIGMA[, , burnin:(R/keep)], c(1, 2), mean)   #���U�����U�s��̏o�����̎��㕽��
round(rbind(sigma_mu, sigma0), 3)



