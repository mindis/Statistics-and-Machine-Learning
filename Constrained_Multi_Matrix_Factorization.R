#####Constrained Multi Matrix Factorization#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
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
k <- 10   #��ꐔ
hh <- 5000   #���[�U�[��
item <- 2000   #�A�C�e����
context <- 200   #�R���e�L�X�g��

#�R���e�L�X�g�𐶐�
context_list <- list()
prob <- as.numeric(extraDistr::rdirichlet(1, rep(20.0, context)))
for(j in 1:item){
  if(j%%100==0){
    print(j)
  }
  context_list[[j]] <- as.numeric(rmnom(hh, 1, prob) %*% 1:context)
}

##ID�̐ݒ�
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)
context_id0 <- unlist(context_list)


##�����ϐ����Ó��ɂȂ�܂Ńp�����[�^�̐������J��Ԃ�
for(rp in 1:1000){
  print(rp)
  
  ##CMMF���f���̉����ϐ��𐶐�
  sigmat <- sigma <- 0.5   #���f���̊ϑ��덷 
  
  #����t���̓����s��
  W <- WT <- mvrnorm(hh, rep(0.675, k), diag(0.25, k))
  H <- HT <- mvrnorm(item, rep(0.5, k), diag(0.2, k))
  V <- VT <- mvrnorm(context, rep(0.25, k), diag(0.15, k))
  
  #���f���̕��ύ\�����牞���ϐ��𐶐�
  mu <- rowSums(W[user_id0, ] * H[item_id0, ]) + rowSums(W[user_id0, ] * V[context_id0, ])
  y0 <- rtnorm(hh*item, mu, sigma, a=1, b=10)
  
  #�����ϐ���break����
  if(mean(y0) < 6.5 & mean(y0) > 4.5) break
}

#���������X�R�A��]���f�[�^�ɕϊ�
y0_censor <- ifelse(y0 < 1, 1, ifelse(y0 > 10, 10, y0)) 
y_vec <- round(y0_censor, 0)   #�X�R�A���ۂ߂�

##�����x�N�g���𐶐�
#�����L���̃x�[�^���z�̃p�����[�^��ݒ�
beta1 <- rbeta(hh, 8.5, 10.0)   #���[�U-�w���m��
beta2 <- rbeta(item, 6.5, 8.0)   #�A�C�e���w���m��

#����������w���f�[�^�𐶐�
Z <- matrix(0, nrow=hh, ncol=item)
for(j in 1:item){
  deficit <- rbinom(hh, 1, beta1 * beta2[j])
  Z[, j] <- deficit   #��������
}

#�����C���f�b�N�X
z_vec <- as.numeric(t(Z))
index_z1 <- which(z_vec==1)
index_z0 <- which(z_vec==0)
N <- length(index_z1)

#�����x�N�g���ɉ����ăf�[�^�𒊏o
user_id <- user_id0[index_z1]
item_id <- item_id0[index_z1]
context_id <- context_id0[index_z1]
y <- y_vec[index_z1]
n1 <- plyr::count(user_id)$freq
n2 <- plyr::count(item_id)$freq
n3 <- plyr::count(context_id)$freq

#�������������ϐ��̃q�X�g�O����
hist(y0, col="grey", xlab="�X�R�A", main="���[�U�[�~�A�C�e���̃X�R�A���z")   #���f�[�^
hist(y_vec, col="grey", xlab="�X�R�A", main="���[�U�[�~�A�C�e���̃X�R�A���z")   #���S�f�[�^�̃X�R�A���z
hist(y, col="grey", xlab="�X�R�A", main="���[�U�[�~�A�C�e���̃X�R�A���z")   #�w���f�[�^�̃X�R�A���z


####�}���R�t�A�������e�J�����@��CMMF�𐄒�####
##�A���S���Y���̐ݒ�
R <- 2000
keep <- 2
disp <- 10
iter <- 0

##���O���z�̐ݒ�
theta <- rep(0, k)
tau <- 100 * diag(k)
inv_tau <- solve(tau)
s0 <- 1.0
v0 <- 1.0

##�^�l�̐ݒ�
W <- WT
H <- HT
V <- VT
sigma <- sigmat

##�����l�̐ݒ�
sigma <- sd(y)
W <- mvrnorm(hh, rep(0, k), 0.2 * diag(k))
H <- mvrnorm(item, rep(0, k), 0.2 * diag(k))
V <- mvrnorm(context, rep(0, k), 0.2 * diag(k))

##�p�����[�^�̊i�[�p�z��
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(item, k, R/keep))
V_array <- array(0, dim=c(context, k, R/keep))
SIGMA <- rep(0, R/keep)


##�C���f�b�N�X��ݒ�
user_index <- item_index <- context_index <- list()
ui_id <- uc_id <- iu_id <- ic_id <- cu_id <- ci_id <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
  ui_id[[i]] <- item_id[user_index[[i]]]
  uc_id[[i]] <- context_id[user_index[[i]]]
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
  iu_id[[j]] <- user_id[item_index[[j]]]
  ic_id[[j]] <- context_id[item_index[[j]]]
}
for(j in 1:context){
  context_index[[j]] <- which(context_id==j)
  cu_id[[j]] <- user_id[context_index[[j]]]
  ci_id[[j]] <- item_id[context_index[[j]]]
}
vec <- rep(1, k)
const <- hh / 100


##�ΐ��ޓx�̊�l
LLst <- sum(dnorm(y, mean(y), sd(y), log=TRUE))


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##���[�U�[�̓����s����T���v�����O
  for(i in 1:hh){
    #�����x�N�g���̃p�����[�^
    X <- H[ui_id[[i]], ] + V[uc_id[[i]], ]
    inv_XXV <- solve(t(X) %*% X + inv_tau)
    beta_mu <- inv_XXV %*% t(X) %*% y[user_index[[i]]]   #���ϗʐ��K���z�̕��σx�N�g��
    cov <- sigma^2 * inv_XXV
    
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    W[i, ] <- mvrnorm(1, beta_mu, cov)  
  }

  
  ##�A�C�e���̓����s����T���v�����O
  W_vec <- W[user_id, ]
  er <- as.numeric(y - (W_vec * V[context_id, ]) %*% vec)   #�덷��ݒ�
  
  for(j in 1:item){
    #�����x�N�g���̃p�����[�^
    X <- W[iu_id[[j]], ]
    inv_XXV <- solve(t(X) %*% X + inv_tau)
    beta_mu <- inv_XXV %*% t(X) %*% er[item_index[[j]]]   #���ϗʐ��K���z�̕��σx�N�g��
    cov <- sigma^2 * inv_XXV
    
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    H[j, ] <- mvrnorm(1, beta_mu, cov)
  }
  
  ##�R���e�L�X�g�̓����s����T���v�����O
  er <- as.numeric(y - (W_vec * H[item_id, ]) %*% vec)   #�덷��ݒ�
  for(j in 1:context){
    #�����x�N�g���̃p�����[�^
    X <- W[cu_id[[j]], ]
    inv_XXV <- solve(t(X) %*% X + inv_tau)
    beta_mu <- inv_XXV %*% t(X) %*% er[context_index[[j]]]   #���ϗʐ��K���z�̕��σx�N�g��
    cov <- sigma^2 * inv_XXV
    
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    V[j, ] <- mvrnorm(1, beta_mu, cov)
  }
  
  ##���f���p�����[�^���T���v�����O
  #���f���̕W���΍����T���v�����O
  WHV <- as.numeric((W_vec * H[item_id, ]) %*% vec + (W_vec * V[context_id, ]) %*% vec)
  er <- y - WHV
  s <- s0 + t(er) %*% er
  v <- v0 + N
  sigma <- sqrt(1/(rgamma(1, v/2, s/2)))   #�t�K���}���z����sigma���T���v�����O
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    #�p�����[�^�̊i�[
    mkeep <- rp/keep
    W_array[, , mkeep] <- W
    H_array[, , mkeep] <- H
    V_array[, , mkeep] <- V
    SIGMA[mkeep] <- sigma
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    LL <- sum(dnorm(y, WHV, sigma, log=TRUE))
    print(rp)
    print(c(LL, LLst))
    print(round(c(sigma, sigmat), 3))
  }
}

matplot(t(W_array[1, , ]), type="l")
