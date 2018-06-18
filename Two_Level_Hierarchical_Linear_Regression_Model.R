#####Two Level Hierarchical Linear Regression Model#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(actuar)
library(extraDistr)
library(caret)
library(dplyr)
library(foreach)
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
item <- 2500   #�A�C�e����
year <- 10   #�ϑ��N��

#�A�C�e�����Ƃ̊ϑ�����ݒ�
item_id0 <- rep(1:item, rep(year, item))
year_id0 <- rep(1:year, item)

#�|�A�\�����z����ϑ����𐶐�
repeat {
  alpha <- rnorm(item, 1.7, 1.0)
  beta <- sort(rnorm(year, 0, 1.0))
  n0 <- rep(0, item*year)
  
  for(j in 1:year){
    lambda <- exp(alpha + beta[j])
    n0[year_id0==j] <- rpois(item, lambda)
  }
  n0[n0 < 0] <- 0
  if(max(n0) <= 500 & min(plyr::count(rep(item_id0, n0))$freq) >= 2 & sum(n0==0) >= 2500) break
}

##��W�v�f�[�^��ID��ݒ�
item_id <- rep(item_id0, n0)
pt_id <- rep(year_id0, n0)
id0 <- paste(item_id, pt_id, sep="-")
n_id <- left_join(data.frame(id=id0, stringsAsFactors=FALSE),
                  data.frame(id=unique(id0), no=1:length(unique(id0)), stringsAsFactors=FALSE), by="id")$no
N <- length(n_id)
context <- length(unique(n_id))


##��W�v���x���̐����ϐ��𐶐�
#���[�U�[�̐����ϐ�
k1 <- 4; k2 <- 5; k3 <- 5
u1 <- matrix(runif(N*k1, 0, 1), nrow=N, ncol=k1)
u2 <- matrix(0, nrow=N, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(N, 1, pr)
}
u3 <- rmnom(N, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #�f�[�^������

##�����ϐ����Ó��Ȓl�ɂȂ�܂Ŕ���������
repeat {

  ##�p�����[�^��ݒ�
  #��A�׃N�g���𐶐�
  beta0 <- 5.5   #�ؕ�
  beta1 <- rnorm(ncol(u)-1, 0, 0.7)
  beta <- betat <- c(beta0, beta1)
  
  #�}���`���x���ϗʌ��ʂ𐶐�
  tau1 <- taut1 <- 1.5
  tau2 <- taut2 <- 0.4
  theta1 <- thetat1 <- rnorm(item, 0, tau1)
  theta2 <- thetat2 <- rnorm(unique(n_id), 0, tau2)
  theta <- thetat <- theta1[item_id] + theta2[n_id]
  
  ##�����ϐ��𐶐�
  sigma <- sigmat <- 0.55
  y0 <- u %*% beta + theta + rnorm(N, 0, sigma)
  
  #���[�v�I������
  if(min(y0) > -2.0 & max(y0) < 13.0) break
}

##�f�[�^��1~10�͈̔͂Ɏ��߂�
y0[y0 > 10] <- 10; y0[y0 < 1] <- 1
y <- as.numeric(round(y0))
hist(y0, col="grey", main="��W�v���x���̕]���X�R�A�̎����l���z")
hist(y, col="grey", main="��W�v���x���̕]���X�R�A�̐����l���z")


####�}���R�t�A�������e�J�����@��multilevel model�𐄒�####
##�A���S���Y���̐ݒ�
LL1 <- -100000000   #�ΐ��ޓx�̏����l
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##�f�[�^�ƃC���f�b�N�X��ݒ�
#�C���f�b�N�X�̐ݒ�
index_n <- list()
index_item <- list()
n_item <- rep(0, item)
n_context <- rep(0, context)
id_vec <- rep(0, context)

for(i in 1:context){
  index_n[[i]] <- which(n_id==i)
  n_context[i] <- length(index_n[[i]])
  id_vec[i] <- item_id[index_n[[i]]][1]
}
for(i in 1:item){
  index_item[[i]] <- which(item_id==i)
  n_item[i] <- length(index_item)
}

#�f�[�^�̐ݒ�
uu <- t(u) %*% u
inv_uu <- solve(uu)

##���O���z�̐ݒ�
alpha0 <- rep(0, ncol(u))
tau0 <- 100 * diag(ncol(u))
inv_tau0 <- solve(tau0)
s0 <- 0.01
v0 <- 0.01

##�����l�̐ݒ�
#�̓����f���̏����l
sigma <- 0.5
beta <- rep(0, ncol(u))

#�ϗʌ��ʂ̏����l
tau1 <- 0.5
tau2 <- 0.25
theta1 <- rnorm(item, 0, tau1)
theta2 <- rnorm(unique(n_id), 0, tau2)
theta <- thetat <- theta1[item_id] + theta2[n_id]


##�p�����[�^�̊i�[�p�z��
BETA <- matrix(0, nrow=R/keep, ncol=ncol(u))
THETA1 <- matrix(0, nrow=R/keep, ncol=item)
THETA2 <- matrix(0, nrow=R/keep, ncol=context)
COV <- matrix(0, nrow=R/keep, ncol=3)

##�ΐ��ޓx�̊�l
#���ύ\�����f���̑ΐ��ޓx
LLst1 <- sum(dnorm(y, mean(y), sd(y), log=TRUE))

#���`��A���f���̑ΐ��ޓx
beta0 <- solve(t(u) %*% u) %*% t(u) %*% y
mu <- u %*% beta0
LLst2 <- sum(dnorm(y, as.numeric(mu), sqrt(sum((y-mu)^2)/N), log=TRUE))



####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�̓���A�x�N�g�����T���v�����O
  #�����ϐ���ݒ�
  y_er <- y - theta1[item_id] - theta2[n_id]
  
  #��A�x�N�g���̎��㕪�z�̃p�����[�^
  Xy <- t(u) %*% y_er
  XXV <- uu + inv_tau0
  inv_XXV <- solve(XXV)
  sigma_par <- sigma^2 * inv_XXV
  mu_par <- inv_XXV %*% (Xy + inv_tau0 %*% alpha0)
  
  #���ϗʐ��K���z����A�x�N�g�����T���v�����O
  beta <- mvrnorm(1, mu_par, sigma_par)
  u_mu <- as.numeric(u %*% beta)
  
  ##�̓��W���΍����T���v�����O
  #�t�K���}���z�̃p�����[�^
  s <- s0 + sum((y_er - u_mu)^2)
  v <- v0 + N
  
  #�t�K���}���z���W���΍����T���v�����O
  sigma <- sqrt(1/rgamma(1, v/2, s/2))
  
  
  ##���Ԉˑ��̃A�C�e���ϗʌ��ʂ��T���v�����O
  #�����ϐ���ݒ�
  y_er <- y - u_mu - theta1[item_id]
  
  #�������Ƃ̎��㕪�z�̃p�����[�^
  mu_context <- rep(0, context)
  for(i in 1:context){
    mu_context[i] <- mean(y_er[index_n[[i]]])
  }
  weights <- tau2^2 / (sigma^2/n_context + tau2^2)   #�d�݌W��
  mu_par <- weights*mu_context   #���㕪�z�̕���
  sigma_par <- sqrt(1 / (n_context/sigma^2 + 1/tau2^2))
  
  #���K���z��莖�㕪�z���T���v�����O
  theta2 <- rnorm(context, mu_par, sigma_par)
  
  ##���Ԉˑ��̃A�C�e���ϗʌ��ʂ̕W���΍����T���v�����O
  #�t�K���}���z�̃p�����[�^
  s <- s0 + sum((theta2 - mean(theta2))^2)
  v <- v0 + context
  
  #�t�K���}���z���W���΍����T���v�����O
  tau2 <- sqrt(1/rgamma(1, v/2, s/2))
  
  
  ##�A�C�e���ϗʌ��ʂ��T���v�����O
  #�����ϐ���ݒ�
  y_er <- y - u_mu - theta2[n_id]
  
  #�A�C�e�����Ƃ̎��㕪�z�̃p�����[�^
  mu_item <- rep(0, item)
  for(i in 1:item){
    mu_item[i] <- mean(y_er[index_item[[i]]])
  }
  weights <- tau1^2 / (sigma^2/n_item + tau1^2)   #�d�݌W��
  mu_par <- weights*mu_item   #���㕪�z�̕���
  sigma_par <- sqrt(1 / (n_item/sigma^2 + 1/tau1^2))
  
  #���K���z��莖�㕪�z���T���v�����O
  theta1 <- rnorm(item, mu_par, sigma_par)
  
  ##�A�C�e���ϗʌ��ʂ̕W���΍����T���v�����O
  #�t�K���}���z�̃p�����[�^
  s <- s0 + sum((theta1 - mean(theta1))^2)
  v <- v0 + item
  
  #�t�K���}���z���W���΍����T���v�����O
  tau1 <- sqrt(1/rgamma(1, v/2, s/2))
  
  
  ##�T���v�����O���ʂ̊i�[�ƕ\��
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
    THETA1[mkeep, ] <- theta1
    THETA2[mkeep, ] <- theta2
    COV[mkeep, ] <- c(sigma, tau1, tau2)
  }
  
  #�T���v�����O���ʂ̕\��
  if(rp%%disp==0){
    #�ΐ��ޓx�𐄒�
    mu <- u_mu + theta1[item_id] + theta2[n_id]
    LL <- sum(dnorm(y, mu, sigma, log=TRUE))
    
    #�T���v�����O���ʂ̕\��
    print(rp)
    print(c(LL, LLst1, LLst2))
    print(round(rbind(beta, betat), 3))
    print(round(c(sigma, tau1, tau2), 3))
  }
}

####�T���v�����O���ʂ̉����Ɨv��####
##�T���v�����O���ʂ̉���
matplot(BETA, type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="��A�x�N�g���̃T���v�����O����")
matplot(COV, type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�W���΍��̃T���v�����O����")
matplot(THETA1[, 1:20], type="l", xlab="�T���v�����O��", ylab="�p�����[�^", 
        main="�A�C�e���ϗʌ��ʂ̃T���v�����O����")
matplot(THETA1[, 101:120], type="l", xlab="�T���v�����O��", ylab="�p�����[�^", 
        main="�A�C�e���ϗʌ��ʂ̃T���v�����O����")
matplot(THETA1[, 501:520], type="l", xlab="�T���v�����O��", ylab="�p�����[�^", 
        main="�A�C�e���ϗʌ��ʂ̃T���v�����O����")
matplot(THETA1[, 1001:1020], type="l", xlab="�T���v�����O��", ylab="�p�����[�^", 
        main="�A�C�e���ϗʌ��ʂ̃T���v�����O����")
matplot(THETA2[, 1:20], type="l", xlab="�T���v�����O��", ylab="�p�����[�^", 
        main="���Ԉˑ��A�C�e���ϗʌ��ʂ̃T���v�����O����")
matplot(THETA2[, 1001:1020], type="l", xlab="�T���v�����O��", ylab="�p�����[�^", 
        main="���Ԉˑ��A�C�e���ϗʌ��ʂ̃T���v�����O����")
matplot(THETA2[, 5001:5020], type="l", xlab="�T���v�����O��", ylab="�p�����[�^",
        main="���Ԉˑ��A�C�e���ϗʌ��ʂ̃T���v�����O����")
matplot(THETA2[, 10000:10021], type="l", xlab="�T���v�����O��", ylab="�p�����[�^",
        main="���Ԉˑ��A�C�e���ϗʌ��ʂ̃T���v�����O����")

##���㕽�ς��v�Z
burnin <- 1000/keep
RS <- R/keep

#�p�����[�^���Ƃ̎��㕽��
beta <- colMeans(BETA[burnin:RS, ])   #��A�x�N�g���̎��㕽��
theta1 <- colMeans(THETA1[burnin:RS, ])   #�A�C�e���ϗʌ��ʂ̎��㕽��
theta2 <- colMeans(THETA2[burnin:RS, ])   #���Ԉˑ��̃A�C�e���ϗʌ��ʂ̎��㕽��
cov <- colMeans(COV[burnin:RS, ])   #�W���΍��̎��㕽��

#���f���̓K���x
mu <- u %*% beta + theta1[item_id] + theta2[n_id]   #���ύ\��
sum(dnorm(y, mu, cov[1], log=TRUE))   #�ΐ��ޓx
sum((y - mu)^2)   #���덷

#���茋�ʂƐ^�l�̔�r
round(cbind(y, mu, abs(y-mu)), 3)
round(cbind(beta, betat), 3)
round(cbind(cov, c(sigmat, taut1, taut2)), 3)
round(cbind(theta1, thetat1), 3)
round(cbind(theta2, thetat2), 3)


