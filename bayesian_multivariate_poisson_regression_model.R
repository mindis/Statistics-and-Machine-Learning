#####���ϗʃ|�A�\����A���f��#####
library(MASS)
library(nlme)
library(glmm)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(98437)

####�C�ӂ̕��U�����U�s����쐬������֐�####
##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
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
N <- 4000   #�T���v����
k <- 4   #�����ϐ���

####�����ϐ��̔���####
cont1 <- 3; bin1 <- 3; multi1 <- 4
X.cont <- matrix(rnorm(N*cont1), nrow=N, ncol=cont1)
X.bin <- matrix(0, nrow=N, ncol=bin1)
X.multi <- matrix(0, nrow=N, ncol=multi1)

#��l�����ϐ���ݒ�
for(i in 1:bin1){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(N, 1, p)
}

#���l�����ϐ���ݒ�
p <- runif(multi1)
X.multi <- t(rmultinom(N, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))]   #�璷�ȕϐ��͍폜

#�f�[�^������
X <- cbind(1, X.cont, X.bin, X.multi)


####�����ϐ��̔���####
##�p�����[�^�̐ݒ�
#��A�p�����[�^�̐ݒ�
beta01 <- runif(k, 0.2, 1.2)
beta02 <- matrix(runif(k*cont1, 0, 0.6), nrow=cont1, ncol=k)
beta03 <- matrix(runif(k*bin1, -0.5, 0.6), nrow=bin1, ncol=k)
beta04 <- matrix(runif(k*(multi1-1), -0.5, 0.8), nrow=multi1-1, ncol=k)
beta0 <- rbind(beta01, beta02, beta03, beta04)
rownames(beta0) <- 1:nrow(beta0)

#���U�����U�p�����[�^�̐ݒ�
corM <- corrM(col=k, lower=-0.5, upper=0.8, eigen_lower=0.01, eigen_upper=0.2)   #���֍s����쐬
Sigma <- covmatrix(col=k, corM=corM, lower=0.5, upper=0.75)   #���U�����U�s��
Cov0 <- Sigma$covariance
cov2cor(Cov0)

##�����ϐ��̔���
#���ϗʐ��K���z���lambda�𔭐�
mu <- X %*% beta0
lambda <- t(apply(mu, 1, function(x) mvrnorm(1, x, Cov0)))

#�|�A�\�����z��艞���ϐ��𔭐�
lambda_exp <- exp(lambda)
Y <- apply(lambda_exp, 2, function(x) rpois(N, x))

#���������������̕��z
hist(Y[, 1], col="grey", main="���������������ϐ��̕��z1", xlab="�����ϐ�")
hist(Y[, 2], col="grey", main="���������������ϐ��̕��z1", xlab="�����ϐ�")
hist(Y[, 3], col="grey", main="���������������ϐ��̕��z1", xlab="�����ϐ�")
hist(Y[, 4], col="grey", main="���������������ϐ��̕��z1", xlab="�����ϐ�")


####�}���R�t�A�������e�J�����@�ő��ϗʃ|�A�\����A���f���𐄒�####
##�|�A�\�����z�̑ΐ��ޓx�֐����`
fr <- function(theta, Y, Y_factorial){
  #�p�����[�^�̐ݒ�
  lambda <- exp(theta)
  
  #�ΐ��ޓx���v�Z
  LLi <- rowSums(Y*log(lambda)-lambda - Y_factorial)
  LL <- sum(LLi)
  LL_val <- list(LLi=LLi, LL=LL)
  return(LL_val)
}

##�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- c()   #�ΐ��ޓx�̕ۑ��p
Y_factorial <- lfactorial(Y)   #Y�̊K��̑ΐ����v�Z���Ă���

##���O���z�̐ݒ�
#�K�w���f���̎��O���z
Deltabar <- matrix(0, nrow=ncol(X), ncol=k)
Adelta <- 0.01 * diag(ncol(X))
nu <- k + ncol(X)
V <- nu * diag(k)

##�T���v�����O���ʂ̕ۑ��p�z��
BETA <- matrix(0, nrow=R/keep, ncol=k*ncol(X))
SIGMA <- matrix(0, nrow=R/keep, ncol=k^2)
THETA <- matrix(0, nrow=R/keep, ncol=N)

##�����l�̐ݒ�
#��A�p�����[�^�̏����l��ݒ�
beta1 <- runif(k, -0.6, 1.2)
beta2 <- matrix(runif(k*cont1, 0, 0.6), nrow=cont1, ncol=k)
beta3 <- matrix(runif(k*bin1, -0.6, 0.8), nrow=bin1, ncol=k)
beta4 <- matrix(runif(k*(multi1-1), -0.6, 0.9), nrow=multi1-1, ncol=k)
oldbeta <- rbind(beta01, beta02, beta03, beta04)
rownames(oldbeta) <- 1:nrow(beta0)

#���U�����U�p�����[�^�̏����l��ݒ�
corM <- corrM(col=k, lower=-0.4, upper=0.4, eigen_lower=0.01, eigen_upper=0.2)   #���֍s����쐬
Sigma <- covmatrix(col=k, corM=corM, lower=0.3, upper=0.45)   #���U�����U�s��
oldcov <- Sigma$covariance
cov_inv <- solve(oldcov)
cov2cor(oldcov)

#lambda�̔���
#���ϗʐ��K���z���lambda�𔭐�
theta_mu <- X %*% beta0
oldtheta <- t(apply(theta_mu, 1, function(x) mvrnorm(1, x, Cov0)))


####MCMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){

  ##MH�@�ŃT���v�����Ƃ�theta���T���v�����O
  thetad <- oldtheta
  thetan <- thetad + 0.3 * mvrnorm(1, rep(0, k), diag(k))
  
  #���O���z�̌덷���v�Z
  er_new <- thetan - theta_mu
  er_old <- thetad - theta_mu
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew <- fr(thetan, Y, Y_factorial)$LLi
  logold <- fr(thetad, Y, Y_factorial)$LLi
  logpnew <- apply(er_new, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  logpold <- apply(er_old, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  
  ##MH�T���v�����O
  #�T���v�����O���̑����邩�ǂ���������
  rand <- runif(N)   #��l�������痐���𔭐�
  LLind.diff <- exp(lognew + logpnew - logold - logpold)   #�̑𗦂��v�Z
  alpha <- ifelse(LLind.diff > 1, 1, LLind.diff)
  alpha <- matrix(alpha, nrow=N, ncol=k)
  
  #alpha�Ɋ�Â�beta���̑�
  oldtheta.r <- ifelse(alpha > rand, thetan, thetad)   #alpha��rand�������Ă�����̑�
  adopt <- sum(oldtheta[, 1]!=oldtheta.r[, 1])/N   #�̑�
  oldtheta <- oldtheta.r   #�p�����[�^���X�V
  
  
  ##���ϗʉ�A���f���ɂ��K�w���f���̃T���v�����O
  out <- rmultireg(Y=oldtheta, X=X, Bbar=Deltabar, A=Adelta, nu=nu, V=V)
  oldbeta <- out$B
  oldcov <- out$Sigma
  
  #�K�w���f���̃p�����[�^���X�V
  cov_inv <- solve(oldcov)
  theta_mu <- X %*% oldbeta
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep  
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    
    print(rp)
    print(round(sum(lognew), 2))
    print(round(cbind(oldbeta, beta0), 3))
    print(round(cbind(oldcov, Cov0), 3))
    print(round(cbind(cov2cor(oldcov), cov2cor(Cov0)), 3))
    print(round(adopt, 3))
  }
}

matplot(BETA[, 1:5], type="l")
matplot(BETA[, 6:10], type="l")
matplot(BETA[, 11:15], type="l")
matplot(BETA[, 16:20], type="l")
matplot(BETA[, 2:10], type="l")
matplot(SIGMA[, 1:4], type="l")
matplot(SIGMA[, 5:8], type="l")
matplot(SIGMA[, 9:12], type="l")
matplot(SIGMA[, 13:16], type="l")