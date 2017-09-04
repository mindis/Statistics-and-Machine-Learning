#####�T���v���Z���N�V�������f��#####
library(MASS)
library(sampleSelection)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####�C�ӂ̕��U�����U�s����쐬������֐�####
##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
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


####�����ϐ��̔���####
N <- 5000   #�T���v����

##���p�������邩�ǂ��������肷������ϐ�
cont1 <- 3; bin1 <- 4; multi1 <- 4
Z.cont <- matrix(rnorm(N*cont1), nrow=N, ncol=cont1)
Z.bin <- matrix(0, nrow=N, ncol=bin1)
Z.multi <- matrix(0, nrow=N, ncol=multi1)

#��l�����ϐ���ݒ�
for(i in 1:bin1){
  p <- runif(1, 0.3, 0.7)
  Z.bin[, i] <- rbinom(N, 1, p)
}

#���l�����ϐ���ݒ�
p <- runif(multi1)
Z.multi <- t(rmultinom(N, 1, p))
Z.multi <- Z.multi[, -which.min(colSums(X.multi))] #�璷�ȕϐ��͍폜

#�f�[�^������
Z <- cbind(1, Z.cont, Z.bin, Z.multi)


##���p�����Ȃ�ǂꂭ�炢���p�����̂������肷������ϐ�
#�A���ϐ��̔���
cont2 <- 3
X.cont <- matrix(rnorm(N*cont2, 0, 1), nrow=N, ncol=cont2)

#��l�ϐ��̔���
bin2 <- 5
X.bin <- matrix(0, nrow=N, ncol=bin2)
for(i in 1:bin2){
  X.bin[, i] <- rbinom(N, 1, runif(1, 0.35, 0.65))
}

#�f�[�^�̌���
X <- cbind(1, X.cont, X.bin)


####��ϗʐ��K���z���牞���ϐ��𔭐�####
#���U�����U�s���ݒ�
sigma0 <- 1.5
rho0 <- 0.6
Cov0 <- matrix(c(1, rho0*sqrt(sigma0), rho0*sqrt(sigma0), sigma0), nrow=2, ncol=2)


for(i in 1:1000){
  print(i)
  #�p�����[�^��ݒ�
  alpha0 <- c(runif(1, -0.8, -0.5), runif(cont1, 0, 0.8), runif(bin1+multi1-1, -0.5, 1.0))
  beta0 <- c(runif(1, 0.3, 0.5), runif(cont2, 0, 0.5), runif(bin2, -0.3, 0.5))
  
  #�����ϐ��̕��ύ\��
  y1 <- Z %*% alpha0
  y2 <- X %*% beta0
  
  #���ϗʐ��K���z���牞���ϐ��𔭐�
  Y <- t(apply(cbind(y1, y2), 1, function(x) mvrnorm(1, x, Cov0)))
  y2 <- exp(Y[, 2])
  y1 <- ifelse(Y[, 1] > 0, 1, 0)
  
  if(max(y2) < 150) break
}
summary(cbind(y1, y2))
cor(Y)

####�Ŗޖ@�ŃT���v���Z���N�V�������f���𐄒�####
##�T���v���Z���N�V�������f���̑ΐ��ޓx���`
loglike <- function(theta, Z, X, y1, y2, k1, k2){

  #�p�����[�^�̐ݒ�  
  alpha <- theta[1:k1]
  beta <- theta[(k1+1):(k2)]
  sigma <- theta[length(theta)-1]
  rho <- theta[length(theta)]
  
  #���ύ\��
  alphaZ <- Z %*% alpha
  betaX <- X %*% beta
  
  #�ΐ��ޓx�̌v�Z
  LL1 <- sum((1-y1) * log(pnorm(-alphaZ, 0, 1))) 
  LL2 <- sum(y1 * log(1/sigma * dnorm((y2-betaX)/sigma, 0, 1)))
  LL3 <- sum(y1 * log(pnorm((alphaZ + rho*(y2-betaX)/sigma)/sqrt(1-rho^2))))
  LL <- LL1 + LL2 + LL3
  return(LL)
}

pnorm(-0.2)

##���j���[�g���@�őΐ��ޓx���ő剻
k1 <- ncol(Z)
k2 <- k1 + ncol(X)

for(i in 1:1000){
  print(i)
  #�����p�����[�^�̐ݒ�
  alpha1 <- c(runif(1, -0.8, -0.5), runif(cont1, 0, 0.8), runif(bin1+multi1-1, -0.5, 1.0))
  beta1 <- c(runif(1, 0.3, 0.5), runif(cont2, 0, 0.5), runif(bin2, -0.3, 0.5))
  sigma1 <- runif(1, 1.5, 2.5)
  rho1 <- runif(1, 0.3, 0.7)
  theta0 <- c(alpha1, beta1, sigma1, rho1)

  #���j���[�g���@�őΐ��ޓx���ő剻
  res <- try(optim(theta0, loglike, Z=Z, X=X, y1=y1, y2=log(y2), k1=k1, k2=k2, method="BFGS", 
                   hessian=TRUE, control=list(fnscale=-1, trace=TRUE)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #�G���[����
}

round(rbind(res$par, c(alpha0, beta0, sigma0, rho)), 2)

