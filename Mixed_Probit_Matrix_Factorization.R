#####Mixed Probit Matrix Factorization#####
library(MASS)
library(matrixStats)
library(FAdist)
library(mnormt)
library(NMF)
library(bayesm)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

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
  D <- ifelse(val < 0, val + abs(val) + 0.1, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}


####�f�[�^�̔���####
hh <- 3000
item <- 300
k <- 10   #���ݕϐ���
n1 <- rep(item, hh)
n2 <- rep(hh, item)

##���f���̉���ɏ]���f�[�^�𐶐�
#���֍s���ݒ�
Cor <- Cort <- corrM(k, -0.7, 1.0)

#���[�U�[�A�C�e�������s��̃p�����[�^�𐶐�
W <- WT <- mvrnorm(hh, rep(0, k), Cor)
H <- HT <- matrix(0, nrow=k, ncol=item)
for(j in 1:k){
  par <- rbeta(1, 20.0, 20.0)
  H[j, ] <- HT[j, ] <- rbinom(item, 1, par)
}

#���[�U�[�ƃA�C�e���̕ϗʌ��ʃp�����[�^�𐶐�
theta1 <- thetat1 <- -1.25
theta2 <- thetat2 <- -1.25
sigma1 <- sigmat1 <- 1
sigma2 <- sigmat2 <- 1
alpha <- alphat <- rnorm(hh, theta1, sigma1)
beta <- betat <- rnorm(item, theta2, sigma2)


#�v���r�b�g���f������A�C�e���w���s��𐶐�
alpha_matrix <- matrix(alpha, nrow=hh, ncol=item)
beta_matrix <- matrix(beta, nrow=hh, ncol=item, byrow=T)
Util <- alpha_matrix + beta_matrix + W %*% H   #���p�֐�
y0 <- rnorm(hh*item, as.numeric(Util), 1)   #���K���z����w���L���𐶐�
y <- ifelse(y0 > 0, 1, 0)
Data <- matrix(y, nrow=hh, ncol=item)
storage.mode(Data) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")   #�X�p�[�X�s��


####�}���R�t�A�������e�J�����@��Mixed Probit Matrix Factorization�𐄒�####
##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b, hh, item){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(matrix(runif(length(mu)), nrow=hh, ncol=item)*(FB-FA)+FA, mu, sigma))
}


##���ϗʐ��K���z�̏����t�����Ғl�Ə����t�����U���v�Z����֐�
cdMVN <- function(mean, Cov, dependent, U){
  
  #���U�����U�s��̃u���b�N�s����`
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent, drop=FALSE]
  Cov21 <- Cov[-dependent, dependent, drop=FALSE]
  Cov22 <- Cov[-dependent, -dependent]
  
  #�����t�����U�Ə����t�����ς��v�Z
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mean[, dependent] + t(CDinv %*% t(U[, -dependent] - mean[, -dependent]))   #�����t�����ς��v�Z
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #�����t�����U���v�Z
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}


##�A���S���Y���̐ݒ�
R <- 10000
keep <- 2
disp <- 10
iter <- 0

##���O���z�̐ݒ�
#�ϗʌ��ʂ̕��U�̎��O���z
alpha01 <- 100
beta01 <- 100
alpha02 <- 0.01
beta02 <- 0.01

#�s�񕪉��̃p�����[�^�̎��O���z
nu <- k   #�t�E�B�V���[�ƕ��z�̎��R�x
V <- nu * diag(rep(1, k))   #�t�E�B�V���[�g���z�̃p�����[�^


##�p�����[�^�̐^�l
W <- WT
H <- HT
r <- rowMeans(H)
alpha <- alphat
alpha_matrix <- matrix(alpha, nrow=hh, ncol=item)
beta <- betat
beta_matrix <- matrix(beta, nrow=hh, ncol=item, byrow=T)
util_mu <- alpha_matrix + beta_matrix + W %*% H   #���p�֐�
Cor <- Cort
cov_inv <- solve(Cor)
theta1 <- -1.25
theta2 <- -1.25
sigma1 <- 1
sigma2 <- 1



##�����l�̐ݒ�
#�ϗʌ��ʂ̏����l
user_rnorm <- rnorm(hh, 0, 1)
item_rnorm <- rnorm(item, 0, 1)
user_sort <- order(user_rnorm)
item_sort <- order(item_rnorm)
alpha <- user_rnorm[user_sort][ceiling(rank(rowSums(Data)))]
alpha_matrix <- matrix(alpha, nrow=hh, ncol=item)
beta <- item_rnorm[item_sort][ceiling(rank(colSums(Data)))]
beta_matrix <- matrix(beta, nrow=hh, ncol=item, byrow=T)
Cor <- corrM(k, -0.3, 0.3)
cov_inv <- solve(Cor)
theta1 <- -1
theta2 <- -1
sigma1 <- 1
sigma2 <- 1

#�s�񕪉��̃p�����[�^�̏����l
res <- nmf(Data+0.1, k)
W <- scale(basis(res))
H <- round(coef(res)*k)
H[H > 1] <- 1


##�T���v�����O���ʂ̕ۑ��p�z��
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, item, R/keep))
COR <- array(0, dim=c(k, k, R/keep))
ALPHA <- matrix(0, nrow=R/keep, ncol=hh)
BETA <- matrix(0, nrow=R/keep, ncol=item)
THETA <- matrix(0, nrow=R/keep, ncol=2)
SIGMA <- matrix(0, nrow=R/keep, ncol=2)

##�ؒf�̈���`
a <- ifelse(Data==0, -100, 0)
b <- ifelse(Data==1, 100, 0)
a_vec <- as.numeric(a)
b_vec <- as.numeric(b)
const <- -1/2*log(2*pi)   #�W�����K���z�̑ΐ��ޓx�̒萔



####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�ؒf���K���z���A�C�e���w���̌��p�l���T���v�����O
  util_mu <- alpha_matrix + beta_matrix + W %*% H   #���p�֐�
  util <- rtnorm(util_mu, 1, a, b, hh, item)
  util[is.infinite(util)] <- 0
  
  ##���[�U�[�����s����T���v�����O
  #�A�C�e�������s��̌덷��ݒ�
  W0 <- matrix(0, nrow=hh, ncol=k)
  util_fearture <- util - alpha_matrix - beta_matrix
  
  #�A�C�e�������s��̎��㕪�z�̃p�����[�^
  XX <- H %*% t(H)
  XXV <- solve(XX + cov_inv)
  
  for(i in 1:hh){
    XXb <- H %*% util_fearture[i, ]
    beta_mean <- beta_mean <- XXV %*% XXb
    
    #���ϗʐ��K���z����beta���T���v�����O
    W0[i, ] <- mnormt::rmnorm(1, beta_mean, XXV)
  }
  W <- scale(W0)
  
  ##�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^
  V_par <- V + t(W) %*% W
  Sn <- nu + hh
  
  #�t�E�B�V���[�g���z���番�U�����U�s��𔭐�
  Cor <- cov2cor(rwishart(Sn, solve(V_par))$IW)
  inv_cov <- solve(Cor)
  
  
  ##�A�C�e�������s����T���v�����O
  for(j in 1:k){
    
    #�ϐ��p�^�[����ݒ�
    H1 <- H0 <- H
    H1[j, ] <- 1
    H0[j, ] <- 0
    
    #�p�^�[�����Ƃ̑ΐ��ޓx���v�Z
    WH0 <- W %*% H0
    WH1 <- W %*% H1
    LL_item0 <- colSums(const -1/2*(util_fearture - WH0)^2)   
    LL_item1 <- colSums(const -1/2*(util_fearture - WH1)^2)
    LLi0 <- cbind(LL_item0, LL_item1)
    LLi <- exp(LLi0 - rowMaxs(LLi0)) * matrix(c(1-r[j], r[j]), nrow=item, ncol=2, byrow=T)   #�ޓx�ɕϊ�
    
    #���ݕϐ��̊����m������H���T���v�����O
    z_rate <- LLi[, 2] / rowSums(LLi)  #���ݕϐ��̊����m�� 
    H[j, ] <- rbinom(item, 1, z_rate)
    
    #���������X�V
    r[j] <- mean(H[j, ])
  }
  
  
  ##���[�U�[�ԕϗʌ��ʂ��T���v�����O
  util_user <-  util - beta_matrix - W %*% H
  
  #���[�U�[���Ƃ̕��ς𐄒�
  mu <- rowMeans(util_user)
  
  #���K���z����ϗʌ��ʂ��T���v�����O
  mu_par <- (1/sigma1)/(1/sigma1 + n1)*theta1 + n1/(1/sigma1 + n1)*mu
  sigma_par <-  1 / (1/sigma1 + n1[1])
  alpha <- rnorm(hh, mu_par, sigma_par)
  alpha_matrix <- matrix(alpha, nrow=hh, ncol=item)
  
  ##�K�w���f���̕��ςƕ��U���T���v�����O
  #���ς��T���v�����O
  mu <- mean(alpha)
  mu_par <- (1/alpha01)/(1/alpha01 + hh)*beta01 + hh/(1/alpha01 + hh)*mu
  sigma_par <- sigma1 / (1/alpha01+hh)
  theta1 <- rnorm(1, mu_par, sigma_par)
  
  #���U���T���v�����O
  r0 <- alpha02 + hh
  s0 <- beta02 + (hh-1)*var(alpha)
  sigma1 <- 1/rgamma(1, r0/2, s0/2)
  
  
  ##�A�C�e���ԕϗʌ��ʂ��T���v�����O
  util_item <-  util - alpha_matrix - W %*% H
  
  #���[�U�[���Ƃ̕��ς𐄒�
  mu <- colMeans(util_item)
  
  #���K���z����ϗʌ��ʂ��T���v�����O
  mu_par <- (1/sigma2)/(1/sigma2 + n2)*theta2 + n2/(1/sigma2 + n2)*mu
  sigma_par <-  1 / (1/sigma2 + n2[1])
  beta <- rnorm(item, mu_par, sigma_par)
  beta_matrix <- matrix(beta, nrow=hh, ncol=item, byrow=T)
  
  ##�K�w���f���̕��ςƕ��U���T���v�����O
  #���ς��T���v�����O
  mu <- mean(beta)
  mu_par <- (1/alpha01)/(1/alpha01 + item)*beta01 + item/(1/alpha01 + item)*mu
  sigma_par <- sigma1 / (1/alpha01+item)
  theta2 <- rnorm(1, mu_par, sigma_par)
  
  #���U���T���v�����O
  r0 <- alpha02 + item
  s0 <- beta02 + (item-1)*var(beta)
  sigma2 <- 1/rgamma(1, r0/2, s0/2)
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    W_array[, , mkeep] <- W
    H_array[, , mkeep] <- H
    COR[, , mkeep] <- Cor
    ALPHA[mkeep, ] <- alpha
    BETA[mkeep, ] <- beta
    THETA[mkeep, ] <- c(theta1, theta2)
    SIGMA[mkeep, ] <- c(sigma1, sigma2)
  
    
    #�T���v�����O���ʂ��m�F
    if(rp%%disp==0){
      print(rp)
      print(c(sum(const -1/2*(util - util_mu)^2), sum(const -1/2*(util - mean(util))^2)))
      #print(round(cbind(W[1:5, ], WT[1:5, ]), 2))
      print(cbind(H[, 1:10], HT[, 1:10]))
      print(round(cbind(Cor, Cort), 2))
      print(round(c(theta1, theta2, sigma1, sigma2, thetat1, thetat2, sigmat1, sigmat2), 3))
    }
  }
}
