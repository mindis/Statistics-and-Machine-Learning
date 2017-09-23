#####�����f�[�^�̂���ϗʌ��ʑ��ϗʉ�A���f��#####
library(MASS)
library(nlme)
library(glmm)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(5783)

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
n <- 1000   #�]���Ώې�
g <- rpois(n, rgamma(n, 13.5, 1.1))
g <- ifelse(g==0, 1, g)
hh <- sum(g)   #�]���Ґ�
k <- 8   #�����ϐ���

##id�̐ݒ�
c.id <- rep(1:length(g), g)   #�]���Ώ�ID

u.id <- c()   #���[�U�[ID
for(i in 1:length(g)){ 
  u.id <- c(u.id, 1:g[i])
}

ID <- data.frame(no=1:sum(g), c.id=c.id, u.id=u.id)


####�����ϐ��̔���####
##�K�w���f���̐����ϐ�
cont1 <- 3; bin1 <- 4; multi1 <- 4
X.cont <- matrix(rnorm(hh*cont1), nrow=hh, ncol=cont1)
X.bin <- matrix(0, nrow=hh, ncol=bin1)
X.multi <- matrix(0, nrow=hh, ncol=multi1)

#��l�����ϐ���ݒ�
for(i in 1:bin1){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(hh, 1, p)
}

#���l�����ϐ���ݒ�
p <- runif(multi1)
X.multi <- t(rmultinom(hh, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))] #�璷�ȕϐ��͍폜

#�f�[�^������
X <- cbind(1, X.cont, X.bin, X.multi)


##�ϗʌ��ʂ̃f�U�C���s���ݒ�
Z <- matrix(0, nrow=sum(g), ncol=n)

for(j in 1:n){
  index <- subset(1:nrow(ID), ID$c.id==j)
  Z[index, j] <- 1 
}

####�����ϐ��̔���####
##�p�����[�^�̐ݒ�
#���U�����U�s��̐ݒ�
Cor0 <- corrM(k, -0.6, 0.9, 0.01, 0.2)   #�̓����f���̑��֍s��
Cov0 <- covmatrix(k, Cor0, 0.6, 0.6)$covariance   #���U�����U�s��ɕϊ�
CorH <- diag(runif(k, 0.5, 0.75), k)   #�K�w���f���̕��U�����U�s��

#�ϗʌ��ʂ̐ݒ�
b.random <- matrix(0, nrow=n, ncol=k)
B.Random <- matrix(0, nrow=sum(g), ncol=k)

for(i in 1:n){
  b.random[i, ] <- mvrnorm(1, rep(0, k), CorH)
  B.Random[ID$c.id==i, ] <- matrix(b.random[i, ], nrow=g[i], ncol=k, byrow=T)
}

#�K�w���f���̃p�����[�^��ݒ�
mu_score <- rnorm(k, 3.2, 0.25)   #�X�R�A�̕��ύ\��
b1 <- matrix(runif(k*cont1, 0, 0.6), nrow=cont1, ncol=k)
b2 <- matrix(runif(k*(bin1+multi1-1), -0.7, 0.7), nrow=bin1+multi1-1, ncol=k)

BETA <- rbind(mu_score, b1, b2)   #�p�����[�^������
rownames(BETA) <- c()
BETAT <- BETA


##�����ϐ��̔���
Mu <- X %*% BETA + Z %*% b.random   #���ύ\��
Y.comp <- Mu + mvrnorm(hh, rep(0, k), Cov0)   #���ύ\��+�덷

##�����ϐ�������������
#�����p�����[�^
#�ϗʌ��ʂ̃p�����[�^
CorA <- diag(runif(k, 0.55, 0.8))
a.random <- matrix(0, nrow=n, ncol=k)
for(i in 1:n){a.random[i, ] <- mvrnorm(1, rep(0, k), CorA)}

#�Œ���ʂ̃p�����[�^
a0 <- runif(k, -1.8, -0.9)   #�X�R�A�̕��ύ\��
a1 <- matrix(runif(k*cont1, 0, 0.6), nrow=cont1, ncol=k)
a2 <- matrix(runif(k*(bin1+multi1-1), -0.9, 0.5), nrow=bin1+multi1-1, ncol=k)
alpha0 <- rbind(a0, a1, a2)   #�p�����[�^������
rownames(alpha0) <- c()

#���W�b�g�ƌ����m�����v�Z
logit <- X %*% alpha0 + Z %*% a.random
Pr <- exp(logit)/(1+exp(logit))

#�����ϐ��𔭐�
Z.na <- 1-apply(Pr, 2, function(x) rbinom(length(x), 1, x))
Y <- Y.comp * Z.na
Y <- ifelse(Y==0, NA, Y)

#���������W�v
colSums(is.na(Y)); round(colMeans(is.na(Y)), 2)


####�}���R�t�A�������e�J�����@�ŕϗʌ��ʑ��ϗʉ�A���f���𐄒�####
#�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4

##���O���z�̐ݒ�
#�Œ����(���ϗʉ�A)�̎��O���z
Deltabar <- matrix(0, nrow=ncol(X), ncol=k)   #��A�p�����[�^�̕��ς̎��O���z
Adelta <- 0.01 * diag(1, ncol(X))   #��A�p�����[�^�̕��U�̎��O���z
nu <- (ncol(X)+1)+k   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(k)   #�t�E�B�V���[�g���z�̃p�����[�^


#�ϗʌ��ʂ̎��O���z
Bbar <- rep(0, k)
A <- 0.01 * diag(1, k)
nu.random <- k
V.random <- nu.random * diag(k)

##�T���v�����O���ʂ̊i�[�p�z��
BETA <- matrix(0, nrow=R/keep, ncol=ncol(X)*k)
SIGMA <- matrix(0, nrow=R/keep, ncol=k*k)
Random <- array(0, dim=c(n, k, R/keep))
Cov.Random <- matrix(0, nrow=R/keep, ncol=k)
Mu.random <- matrix(0, nrow=R/keep, ncol=k)

##MCMC����̂��߂̒萔�̌v�Z
mu_random <- matrix(0, nrow=n, ncol=k)
sigma_random <- array(0, dim=c(k, k, n))


##�����l�̐ݒ�
oldbeta <- solve(t(X) %*% X) %*% t(X) %*% Y.comp
oldsigma <- t(Y.comp - X %*% oldbeta) %*% (Y.comp - X %*% oldbeta)/nrow(X)
beta_random <- solve(t(Z) %*% Z) %*% t(Z) %*% (Y.comp - X %*% oldbeta)
cov_random <- t((Y.comp - X %*% oldbeta) - Z %*% beta_random) %*% ((Y.comp - X %*% oldbeta) - Z %*% beta_random)/nrow(Z)


####MCMC�ō������ϗʉ�A���f���𐄒�####
for(rp in 1:R){
  
  ##�M�u�X�T���v�����O�ŌŒ����beta��sigma���T���v�����O
  y.er <- Y.comp - Z %*% beta_random   #�����ϐ��ƕϗʌ��ʂ̌덷���v�Z
  
  #�x�C�W�A�����ϗʉ�A���f���𐄒�
  out <- rmultireg(y.er, X, Deltabar, Adelta, nu, V)   
  oldbeta <- out$B
  oldsigma <- out$Sigma
  
  ##�M�u�X�T���v�����O�ŕϗʌ��ʂ��T���v�����O
  z.er <- Y.comp - X %*% oldbeta
  
  #ID���Ƃɕ��ς��v�Z
  mu <- as.matrix(data.frame(id=ID$c.id, z=z.er) %>%
                    dplyr::group_by(id) %>%
                    dplyr::summarize_each(funs(mean), everything()))[, -1]
  
  #�x�C�Y����̂��߂̌v�Z
  inv_random <- solve(cov_random) 
  inv_sigma <- solve(oldsigma)
  
  for(i in 1:n){
    n.inv_sigma <- g[i]*inv_sigma
    inv_mixed <- solve(inv_random + n.inv_sigma)
    mu_random[i, ] <- inv_mixed %*% (n.inv_sigma %*% mu[i, ])
    sigma_random[, , i] <- inv_mixed
    beta_random[i, ] <- mvrnorm(1, mu_random[i, ], sigma_random[, , i])
  }
  
  #�K�w���f���̕��U���T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^���v�Z
  R <- solve(V.random) + matrix(rowSums(apply(beta_random, 1, function(x) x %*% t(x))), nrow=k, ncol=k, byrow=T)
  Sn <- nu.random + n
  
  #�t�E�B�V���[�g���z����K�w���f���̕��U���T���v�����O
  cov_random <- rwishart(Sn, solve(R))$IW
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- as.numeric(oldbeta)
    SIGMA[mkeep, ] <- as.numeric(oldsigma)
    Random[, , mkeep] <- beta_random
    Cov.Random[mkeep, ] <- diag(cov_random)
    
    print(rp)
    print(round(cbind(oldbeta, BETAT), 2))
    print(round(cbind(cov2cor(oldsigma), Cor0), 2))
    print(round(rbind(diag(cov_random), diag(CorH)), 2))
  }
}

matplot(Cov.Random[, 1:4], type="l")
matplot(SIGMA[, 1:20], type="l")
matplot(BETA[, 1:8], type="l")

round(cbind(mu, mu_random), 2)