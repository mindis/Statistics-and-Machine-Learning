#####���t���背�R�����f�[�V����#####
library(MASS)
library(nlme)
library(ncvreg)
library(glmnet)
library(lars)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)


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
##�f�[�^�̐ݒ�
hh <- 10000
item <- 500

##�A�C�e�����N�s��𔭐�������
k <- 4
Pr <- matrix(0, nrow=hh, ncol=item)
Z <- matrix(0, nrow=hh, ncol=item)
X <- matrix(rnorm(hh*k), nrow=hh, ncol=k)

for(i in 1:item){
  print(i)
  logit <- runif(1, -5.5, -1.0) + X %*% rnorm(k, 0, 0.3) + rnorm(hh, 0, 0.5)
  Pr[, i] <- exp(logit) / (1+exp(logit))
  Z[, i] <- rbinom(hh, 1, Pr[, i])
}

#�������������N�s��̏W�v�Ɖ���
colMeans(Z)
hist(colMeans(Z), col="grey", main="�A�C�e���̏o����", xlab="")

##�����ϐ��𔭐�������
#�����m�����Ó��Ȑ����ɂȂ�܂ŌJ��Ԃ�
#�p�����[�^�̐ݒ�
for(i in 1:1000){
  print(i)
  beta00 <- runif(1, -2.5, -1.5)
  beta01 <- rnorm(item, 0.2, 0.4)
  beta11 <- ifelse(beta01 > 0, beta01-0.2, beta01)
  
  #���W�b�g�Ɗm���̌v�Z
  logit <- beta00 + Z %*% beta11
  P <- exp(logit)/(1+exp(logit))
  y <- rbinom(hh, 1, P)
  if(mean(y) < 0.45 & mean(y) > 0.4) break
}

####���t���背�R�����f�[�V�������f���𐄒�####
##L1���������W�X�e�B�b�N��A���f����CV������A
#�N���X�o���f�[�V�����ōœK�ȃp�����[�^��ݒ�
k <- 5   #CV������
p1 <- 5
p2 <- 10
alpha_par <- seq(0.4, 1, length=p1)
lambda_par <- seq(0.001, 0.02, length=p2)
pars <- matrix(0, nrow=p1*p2, ncol=2)
cv_right <- c()

#�T���v���𕪊������Ă���
split0 <- split(sample(1:length(y), length(y)), 1:k)
index <- c()
split <- c()

for(i in 1:k){
  split <- c(split, split0[[i]])
  index <- c(index, ifelse(split0[[i]] > 0, i, 0))
}

##�N���X�o���f�[�V�����Ńn�C�p�[�p�����[�^�𐄒�
for(i in 1:length(alpha_par)){
  print(i)
  alpha <- alpha_par[i]   #L1��L2�̏d�݃p�����[�^��I��
  for(j in 1:length(lambda_par)){
    r <- (i-1)*p1+j
    lambda <- lambda_par[j]   #�������p�����[�^��I��
    pars[r, ] <- c(alpha, lambda)   #�n�C�p�[�p�����[�^���i�[
    
    ##5�����N���X�o���f�[�V������CV�G���[�����Z�o
    cv_vec <- c() 
    
    for(v in 1:k){
      #�T���v���𕪊�
      z1 <- Z[split[index!=v], ]
      z2 <- Z[split[index==v], ]
      y1 <- y[split[index!=v]]
      y2 <- y[split[index==v]]
      
      #glmnet��L1���������W�X�e�B�b�N��A���f���𐄒�
      res <- glmnet(x=z1, y=y1, lambda=lambda, family="binomial", alpha=alpha)
      alpha <- res$a0
      beta <- res$beta
      
      #CV�G���[���Z�o
      logit <- alpha + z2 %*% beta
      Pr <- exp(logit) / (1+exp(logit))
      cv_vec <- c(cv_vec, y2==ifelse(Pr > 0.5, 1, 0))
    }
    cv_right <- c(cv_right, mean(cv_vec))
    print(cv_right[r])
  }
}
