#####�����������ϗʐ��K���z#####
library(MASS)
library(matrixStats)
library(FAdist)
library(mclust)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####�C�ӂ̕��U�����U�s����쐬������֐�####
##���ϗʐ��K���z����̗����𔭐�������
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
##�f�[�^�̐ݒ�
k <- 12   #������
N <- 100000   #�T���v����
v <- 7   #������
vec <- rep(1, k)

##�p�����[�^�̐ݒ�
#���σx�N�g���𐶐�
mu <- matrix(0, nrow=k, ncol=v)
for(j in 1:k){
  lower <- runif(v, 30, 70)   #�ϐ��̉����l
  upper <- lower + runif(v, 30, 60)   #�ϐ��̏���l
  mu[j, ] <- runif(v, lower, upper)   #�ϐ��̕��ϒl�𔭐�
}
mut <- mu

#���U�����U�s��𐶐�
Cov <- array(0, dim=c(v, v, k))
lower <- 4; upper <- 25
for(j in 1:k){
  Cor <- corrM(v, -0.6, 0.9, 0.1, 0.5)
  Cov[, , j] <- covmatrix(v, Cor, lower, upper)$covariance
}
Covt <- Cov

##�����ϐ��𐶐�
#�Z�O�����g�����𐶐�
prob <- extraDistr::rdirichlet(N, rep(25.0, k))
Z <- rmnom(N, 1, prob)
z_vec <- as.numeric(Z %*% 1:k)
r <- colMeans(Z)

#���ϗʐ��K���z����f�[�^�𐶐�
y <- matrix(0, nrow=N, ncol=v)
for(i in 1:N){
  y[i, ] <- mvrnorm(1, mu[z_vec[i], ], Cov[, , z_vec[i]])
}


####EM�A���S���Y���ō������ϗʐ��K���z�𐄒�####
##���ϗʐ��K���z�̖��x�֐�
mvdnorm <- function(u, mu, Cov, N, s){
  er <- y - matrix(mu, nrow=N, ncol=v, byrow=T)
  Lho <- 1 / (sqrt(2*pi)^s*sqrt(abs(det(Cov)))) * exp(-1/2 * as.numeric((er %*% ginv(Cov) * er) %*% rep(1, s)))
  return(Lho)
}

##�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̒�`
LLobz <- function(y, mu, Cov, r, v, k, vec){
  #�������ϗʐ��K���z�̃Z�O�����g���Ƃ̖ޓx
  LLind <- matrix(0, nrow=N, ncol=k)
  for(j in 1:k){
    LLind[, j] <- mvdnorm(y, mu[j, ], Cov[, , j], N, v)   #���ϗʐ��K���z�̖ޓx
  }

  #�ΐ��ޓx�Ɛ��ݕϐ�z���`
  LLho <- matrix(r, nrow=N, ncol=k, byrow=T) * LLind; LLho_vec <- as.numeric(LLho %*% vec)
  z <- LLho / LLho_vec   #���ݕϐ�z�̊����m��
  LL <- sum(log(LLho_vec))   #�ϑ��f�[�^�̑ΐ��ޓx�̘a
  value <- list(LL=LL, z=z)
  return(value)
}

##EM�A���S���Y���̐ݒ�
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏�����
tol <- 0.001 
iter <- 1

##�p�����[�^�̏����l
Zi <- rmnom(N, 1, rep(20, k))
mu <- matrix(0, nrow=k, ncol=v)
Cov <- array(0, dim=c(v, v, k))
r <- rep(1/k, k)
for(j in 1:k){
  data <- y[Zi[, j]==1, ]
  mu[j, ] <- colMeans(data) + runif(v, -15, 15)
  Cov[, , j] <- var(data)
}

#�ΐ��ޓx�̏�����
Lobz <- LLobz(y, mu, Cov, r, v, k, vec)
LL1 <- Lobz$LL


####EM�A���S���Y���ɂ��p�����[�^�̍X�V####
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  ##M�X�e�b�v�Ńp�����[�^���Ŗސ���
  z <- Lobz$z   #���ݕϐ�z�̏o��
  
  for(j in 1:k){
    #���σx�N�g�����Ŗސ���
    weighted_z <- sum(z[, j])
    weighted_data <- z[, j] * y   #�d�ݕt���̉����ϐ�
    mu[j, ] <- colSums(weighted_data) / weighted_z
    
    #���U�����U�s����Ŗސ���
    er <- weighted_data - (z[, j] * matrix(mu[j, ], nrow=N, ncol=v, byrow=T))   #�d�ݕt���덷
    Cov[, , j] <- t(er) %*% er / weighted_z
  }
  #���������X�V
  r <- colSums(z) / N
  
  ##E�X�e�b�v�Ŋϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̒�`
  Lobz <- LLobz(y, mu, Cov, r, v, k, vec)   #�ϑ��f�[�^�̑ΐ��ޓx
  LL <- Lobz$LL
  
  ##�A���S���Y���̍X�V
  ite <- iter + 1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}
res <- Mclust(y, k)


LLobz(y, mut, Covt, r, v, k, vec)$LL
t(res$parameters$mean)
mut
mu