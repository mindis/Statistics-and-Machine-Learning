#####Collaborative Latent Dirichlet Allocation#####
options(warn=0)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(2506787)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
k <- 20   #�g�s�b�N��
hh <- 5000   #���[�U�[��
item <- 2000   #�A�C�e����
g <- 1000   #��b��
pt <- rtpois(hh, rgamma(hh, 20, 0.2), a=0, b=Inf)   #�]������
hhpt <- sum(pt)   #���]������
w <- extraDistr::rtpois(item, rgamma(item, 30, 0.2), a=0, b=Inf)   #�A�C�e���̒P�ꐔ
f <- sum(w)   #���P�ꐔ


##ID�ƃC���f�b�N�X�̐ݒ�
#ID�̐ݒ�
user_id <- rep(1:hh, pt)
d_id <- rep(1:item, w)

#�C���f�b�N�X�̐ݒ�
user_index <- d_index <- list()
vec <- c(0, cumsum(pt))
for(i in 1:hh){
  user_index[[i]] <- (1:hhpt)[(vec[i]+1):vec[i+1]]
}
vec <- c(0, cumsum(w))
for(j in 1:item){
  d_index[[j]] <- (1:f)[(vec[j]+1):vec[j+1]]
}

##�f���x�N�g���𐶐�
k1 <- 3; k2 <- 5; k3 <- 5
x1 <- matrix(runif(hhpt*k1, 0, 1), nrow=hhpt, ncol=k1)
x2 <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #�f�[�^������

##�K�w���f���̐����ϐ��𐶐�
#���[�U�[�̐����ϐ�
k1 <- 2; k2 <- 4; k3 <- 5
u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #�f�[�^������

#�A�C�e���̐����ϐ�
k1 <- 2; k2 <- 4; k3 <- 4
v1 <- matrix(runif(item*k1, 0, 1), nrow=item, ncol=k1)
v2 <- matrix(0, nrow=item, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  v2[, j] <- rbinom(item, 1, pr)
}
v3 <- rmnom(item, 1, runif(k3, 0.2, 1.25)); v3 <- v3[, -which.min(colSums(v3))]
v <- cbind(1, v1, v2, v3)   #�f�[�^������

#�p�����[�^��
k1 <- ncol(x); k2 <- ncol(u); k3 <- ncol(v)


##�A�C�e���̊����𐶐�
#�Z�O�����g�����𐶐�
topic <- 25
gamma <- extraDistr::rdirichlet(topic, rep(0.5, item))
z <- as.numeric(rmnom(hh, 1, extraDistr::rdirichlet(hh, rep(2.5, topic))) %*% 1:topic)

#�������z����A�C�e���𐶐�
item_id_list <- list()
for(i in 1:hh){
  if(i%%100==0){
    print(i)
  }
  item_id_list[[i]] <- as.numeric(rmnom(pt[i], 1, gamma[z[user_id[user_index[[i]]]], ]) %*% 1:item)
}
item_id <- unlist(item_id_list)
item_index <- list()
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
}


##�����ϐ����Ó��ɂȂ�܂Ńp�����[�^�̐������J��Ԃ�
rp <- 0
repeat { 
  rp <- rp + 1
  
  ##LDA�̃p�����[�^�𐶐�
  #LDA�̃f�B���N�����z�̃p�����[�^��ݒ�
  alpha01 <- rep(0.1, k)
  alpha02 <- matrix(0.015, nrow=k, ncol=g)
  for(j in 1:k){
    alpha02[j, matrix(1:item, nrow=k, ncol=g/k, byrow=T)[j, ]] <- 0.125
  }
  
  #�f�B���N�����z����p�����[�^�𐶐�
  theta <- thetat <- extraDistr::rdirichlet(item, alpha01)
  phi <- phit <- extraDistr::rdirichlet(k, alpha02)
  
  #�A�C�e���o���m�����Ⴂ�g�s�b�N�����ւ���
  index <- which(colMaxs(phi) < (k*5)/f)
  for(j in 1:length(index)){
    phi[as.numeric(rmnom(1, 1, extraDistr::rdirichlet(1, alpha01)) %*% 1:k), index[j]] <- (k*5)/f
  }
  
  ##�s�񕪉��̃p�����[�^�𐶐�
  #�f���x�N�g���̃p�����[�^
  sigma <- sigmat <- 0.3
  beta <- betat <- c(5.5, rnorm(k1-1, 0, 0.5))
  
  #�K�w���f���̕��U�p�����[�^
  Cov_u <- Cov_ut <- diag(runif(k, 0.005, 0.1), k)   #���[�U�[-�A�C�e���̊K�w���f���̕��U
  Cov_v <- Cov_vt <- diag(runif(k, 0.005, 0.1), k)   #�A�C�e���̊K�w���f���̕��U
  
  #�K�w���f���̉�A�W����ݒ�
  alpha_u <- alpha_ut <- matrix(rnorm(k*k2, 0, 0.3), nrow=k2, ncol=k)
  alpha_v <- alpha_vt <- matrix(rnorm(k*k3, 0, 0.3), nrow=k3, ncol=k)
  
  #�s�񕪉��̃p�����[�^�𐶐�
  theta_u <- theta_ut <- u %*% alpha_u + mvrnorm(hh, rep(0, k), Cov_u)
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, k), Cov_v)
  
  
  ##�g�s�b�N�ƒP��𐶐�
  #�������z����g�s�b�N�𐶐�
  Z <- rmnom(f, 1, theta[d_id, ])
  z_vec <- as.numeric(Z %*% 1:k)
  
  #�g�s�b�N�Ɋ�Â��P��𐶐�
  word <- rmnom(f, 1, phi[z_vec, ])
  wd <- as.numeric(word %*% 1:g)
  
  
  ##���K���z����]���x�N�g���𐶐�
  #�]���X�R�A�̊��Ғl
  vec_topic <- rep(1, k)
  lambda <- as.numeric(x %*% beta)   #�f���x�N�g���̊��Ғl
  uv <- as.numeric((theta_u[user_id, ] * (theta[item_id, ] + theta_v[item_id, ])) %*% vec_topic)   #�s�񕪉��̃p�����[�^
  mu <- lambda + uv   #���Ғl
  
  #�]���x�N�g���𐶐�
  y0 <- rnorm(hhpt, mu, sigma)
  
  #break����
  print(sum(colSums(word)==0))
  print(c(max(y0), min(y0)))
  if(max(y0) < 15.0 & min(y0) > -4.0 & max(y0) > 12.0 & min(y0) < -2.0 & sum(colSums(word)==0)==0){
    break
  }
}

#���������X�R�A��]���f�[�^�ɕϊ�
y0_censor <- ifelse(y0 < 1, 1, ifelse(y0 > 10, 10, y0)) 
y <- round(y0_censor, 0)   #�X�R�A���ۂ߂�

#�X�R�A���z�ƒP�ꕪ�z
hist(colSums(word), col="grey", breaks=25, xlab="�P��p�x", main="�P��p�x���z")
hist(y0, col="grey", breaks=25, xlab="�X�R�A", main="���S�f�[�^�̃X�R�A���z")
hist(y, col="grey", breaks=25, xlab="�X�R�A", main="�ϑ����ꂽ�X�R�A���z")

#�X�p�[�X�s��ɕϊ�
item_data <- sparseMatrix(1:hhpt, item_id, x=rep(1, hhpt), dims=c(hhpt, item))
item_data_T <- t(item_data)
word_data <- sparseMatrix(1:f, wd, x=rep(1, f), dims=c(f, g))
word_data_T <- t(word_data)


####�M�u�X�T���v�����O��Collaborative Latent Dirichlet Allocation�𐄒�####
##�A�C�e�����Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, w, k){
  #���S�W�����v�Z
  Bur <- theta[w, ] * t(phi)[wd, ]   #�ޓx
  Br <- Bur / rowSums(Bur)   #���S��
  r <- colSums(Br) / sum(Br)   #������
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##�A���S���Y���̐ݒ�
R <- 3000
keep <- 2  
iter <- 0
burnin <- 500
disp <- 10

##���O���z�̐ݒ�
#LDA�̃f�B���N�����z�̎��O���z
alpha01 <- 0.1
alpha02 <- 0.1

#�f���x�N�g���̎��O���z
delata <- rep(0, k1)
tau <- 100 * diag(k1)
inv_tau <- solve(tau)

#�t�K���}���z�̎��O���z
s0 <- 1.0
v0 <- 1.0

#���[�U�[�̊K�w���f���̎��O���z
Deltabar1 <- matrix(rep(0, k2*k), nrow=k2, ncol=k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta1 <- 0.01 * diag(rep(1, k2))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu1 <- k2   #�t�E�B�V���[�g���z�̎��R�x
V1 <- nu1 * diag(rep(1, k)) #�t�E�B�V���[�g���z�̃p�����[�^

#�A�C�e���̊K�w���f���̎��O���z
Deltabar2 <- matrix(rep(0, k3*k), nrow=k3, ncol=k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta2 <- 0.01 * diag(rep(1, k3))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu2 <- k3   #�t�E�B�V���[�g���z�̎��R�x
V2 <- nu2 * diag(rep(1, k)) #�t�E�B�V���[�g���z�̃p�����[�^

##�p�����[�^�̐^�l
#LDA�̃p�����[�^
theta <- thetat
phi <- phit

#�f���x�N�g���̃p�����[�^
beta <- betat
sigma <- sigmat
lambda <- as.numeric(x %*% beta)

#�K�w���f���̃p�����[�^
alpha_u <- alpha_ut; Cov_u <- Cov_ut
mu_u <- u %*% alpha_u; inv_Cov_u <- solve(Cov_u)
alpha_v <- alpha_vt; Cov_v <- Cov_vt
mu_v <- v %*% alpha_v; inv_Cov_v <- solve(Cov_v)

#�s�񕪉��̃p�����[�^
theta_u <- theta_ut
theta_v <- theta_vt
uv <- as.numeric((theta_u[user_id, ] * (theta[item_id, ] + theta_v[item_id, ])) %*% vec_topic)


##�p�����[�^�̏����l��ݒ�
#LDA�̃p�����[�^
theta <- extraDistr::rdirichlet(item, rep(2.0, k))
phi <- extraDistr::rdirichlet(k, rep(2.0, g))

#�f���x�N�g���̃p�����[�^
beta <- as.numeric(solve(t(x) %*% x) %*% t(x) %*% y)
sigma <- 0.5
lambda <- as.numeric(x %*% beta)

#�K�w���f���̃p�����[�^
alpha_u <- matrix(0, nrow=k2, ncol=k); Cov_u <- 0.01 * diag(k)
mu_u <- u %*% alpha_u; inv_Cov_u <- solve(Cov_u)
alpha_v <- matrix(0, nrow=k3, ncol=k); Cov_v <- 0.01 * diag(k)
mu_v <- v %*% alpha_v; inv_Cov_v <- solve(Cov_v)

#�s�񕪉��̃p�����[�^
theta_u <- mu_u + mvrnorm(hh, rep(0, k), Cov_u)
theta_v <- mu_v + mvrnorm(item, rep(0, k), Cov_v)
uv <- as.numeric((theta_u[user_id, ] * (theta[item_id, ] + theta_v[item_id, ])) %*% vec_topic)

##�p�����[�^�̊i�[�p�z��
#LDA�̊i�[�p�z��
THETA <- array(0, dim=c(item, k, R/keep))
PHI <- array(0, dim=c(k, g, R/keep))
SEG <- matrix(0, nrow=f, ncol=k)
storage.mode(SEG) <- "integer"

#�s�񕪉��̊i�[�p�z��
BETA <- matrix(0, nrow=R/keep, ncol=k1)
SIGMA <- rep(0, R/keep)
ALPHA_U <- array(0, dim=c(k2, k, R/keep))
ALPHA_V <- array(0, dim=c(k3, k, R/keep))
COV_U <- array(0, dim=c(k, k, R/keep))
COV_V <- array(0, dim=c(k, k, R/keep))
THETA_U <- array(0, dim=c(hh, k, R/keep))
THETA_V <- array(0, dim=c(item, k, R/keep))


##�f�[�^�ƃC���f�b�N�X�̐ݒ�
#�f�[�^�̐ݒ�
xx <- t(x) %*% x + inv_tau; inv_xx <- solve(xx)

#�C���f�b�N�X�̐ݒ�
ui_id <- iu_id <- list()
for(i in 1:hh){
  ui_id[[i]] <- item_id[user_index[[i]]]
}
for(j in 1:item){
  iu_id[[j]] <- user_id[item_index[[j]]]
}

user_data <- sparseMatrix(1:hhpt, user_id, x=rep(1, hhpt), dims=c(hhpt, hh))
user_data_T <- t(user_data)
d_data <- sparseMatrix(1:f, d_id, x=rep(1, f), dims=c(f, item))
d_data_T <- t(d_data)


##�ΐ��ޓx�̊�l��ݒ�
#LDA�̑ΐ��ޓx�̊�l
LLst1 <- sum(d_data %*% log(colSums(d_data)/f))   #1�p�����[�^�̑ΐ��ޓx
LLbest1 <- sum(log((thetat[d_id, ] * t(phit)[wd, ]) %*% vec_topic))   #�x�X�g�ȑΐ��ޓx

##�s�񕪉����f���̑ΐ��ޓx�̊�l
LLst2 <- sum(dnorm(y, mean(y), sd(y), log=TRUE))   #1�p�����[�^���f���̑ΐ��ޓx
mu <- as.numeric(x %*% betat) + as.numeric((theta_ut[user_id, ] * (thetat[item_id, ] + theta_vt[item_id, ])) %*% vec_topic)
LLbest2 <- sum(dnorm(y, mu, sigmat, log=TRUE))   #�x�X�g�ȑΐ��ޓx



####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##LDA�̃p�����[�^���T���v�����O
  #�������z����g�s�b�N���T���v�����O
  Lho <- theta[d_id, ] * t(phi)[wd, ]   #�g�s�b�N�̖ޓx
  prob <- Lho / as.numeric(Lho %*% vec_topic)   #�g�s�b�N�̊����m��
  Zi <- rmnom(f, 1, prob)   #�g�s�b�N���T���v�����O
  z_vec <- as.numeric(Zi %*% 1:k)
  
  #�f�B���N�����z����g�s�b�N���z���T���v�����O
  dsums <- d_data_T %*% Zi + alpha01   #�f�B���N�����z�̃p�����[�^
  theta <- extraDistr::rdirichlet(item, dsums)   #�p�����[�^���T���v�����O
  
  #�f�B���N�����z����P�ꕪ�z���T���v�����O
  wsums <- t(word_data_T %*% Zi) + alpha02   #�f�B���N�����z�̃p�����[�^
  phi <- extraDistr::rdirichlet(k, wsums)   #�p�����[�^���T���v�����O
  
  
  ##�f���x�N�g���̃p�����[�^���T���v�����O
  #�����ϐ��̐���
  y_er <- y - uv   #���f���덷
  
  #���ϗʐ��K���z�̃p�����[�^��ݒ�
  xy <- t(x) %*% y_er
  mu_vec <- as.numeric(inv_xx %*% xy)   #���ϗʐ��K���z�̕��σx�N�g��
  
  #���ϗʐ��K���z����f���x�N�g�����T���v�����O
  beta <- mvrnorm(1, mu_vec, sigma^2*inv_xx)
  lambda <- as.numeric(x %*% beta)
  
  ##���f���̕W���΍����T���v�����O
  #�t�K���}���z�̃p�����[�^
  er <- y - lambda - uv   #���f���̌덷
  s1 <- as.numeric(t(er) %*% er) + s0
  v1 <- hhpt + v0
  
  #�t�K���}���z����W���΍����T���v�����O
  sigma <- sqrt(1/rgamma(1, v1/2, s1/2))
  
  
  ##���[�U�[�̓����s����T���v�����O
  #���f���̉����ϐ�
  y_er <- y - lambda
  
  for(i in 1:hh){
    #�����x�N�g���̎��㕪�z�̃p�����[�^
    X <- theta_v[ui_id[[i]], ] + theta[ui_id[[i]], ]
    Xy <- t(X) %*% y_er[user_index[[i]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_u)
    theta_vec <- inv_XXV %*% (Xy + inv_Cov_u %*% mu_u[i, ])
    
    #���ϗʐ��K���z���烆�[�U�[�����x�N�g�����T���v�����O
    theta_u[i, ] <- mvrnorm(1, theta_vec, sigma^2*inv_XXV)
  }
  
  ##�A�C�e���̓����s����T���v�����O
  #���f���̉����ϐ�
  y_er <- y - lambda - as.numeric((theta_u[user_id, ]*theta[item_id, ]) %*% vec_topic)
  
  for(j in 1:item){
    #�����x�N�g���̎��㕪�z�̃p�����[�^
    X <- theta_u[iu_id[[j]], ]
    Xy <- t(X) %*% y_er[item_index[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_v)
    theta_vec <- inv_XXV %*% (Xy + inv_Cov_v %*% mu_v[j, ])
    
    #���ϗʐ��K���z���烆�[�U�[�����x�N�g�����T���v�����O
    theta_v[j, ] <- mvrnorm(1, theta_vec, sigma^2*inv_XXV)
  } 
  uv <- as.numeric((theta_u[user_id, ] * (theta[item_id, ] + theta_v[item_id, ])) %*% vec_topic)   #�s�񕪉��̃p�����[�^���X�V
  
  
  ##�K�w���f���̃p�����[�^���T���v�����O
  #���[�U�[�̍s�񕪉��̃p�����[�^���T���v�����O
  out_u <- rmultireg(theta_u, u, Deltabar1, ADelta1, nu1, V1)
  alpha_u <- out_u$B
  Cov_u <- diag(diag(out_u$Sigma))
  mu_u <- u %*% alpha_u
  inv_Cov_u <- solve(Cov_u)
  
  #�A�C�e���̍s�񕪉��̃p�����[�^���T���v�����O
  out_v <- rmultireg(theta_v, v, Deltabar2, ADelta2, nu2, V2)
  alpha_v <- out_v$B
  Cov_v <- diag(diag(out_v$Sigma))
  mu_v <- v %*% alpha_v
  inv_Cov_v <- solve(Cov_v)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    mkeep <- rp/keep
    #LDA�̃T���v�����O���ʂ̊i�[
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    
    #�s�񕪉��̃T���v�����O���ʂ̊i�[
    BETA[mkeep, ] <- beta
    THETA_U[, , mkeep] <- theta_u
    THETA_V[, , mkeep] <- theta_v
    ALPHA_U[, , mkeep] <- alpha_u
    ALPHA_V[, , mkeep] <- alpha_v 
    COV_U[, , mkeep] <- Cov_u
    COV_V[, , mkeep] <- Cov_v
  }
  #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
  if(rp%%keep==0 & rp >= burnin){
    SEG <- SEG + Zi
  }
  
  if(rp%%disp==0){
    #LDA�̑ΐ��ޓx
    LL1 <- sum(log((theta[d_id, ] * t(phi)[wd, ]) %*% vec_topic)) 
    
    #�s�񕪉����f���̑ΐ��ޓx
    mu <- as.numeric(x %*% beta) + as.numeric((theta_u[user_id, ] * (theta[item_id, ] + theta_v[item_id, ])) %*% vec_topic)
    LL2 <- sum(dnorm(y, mu, sigma, log=TRUE))
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(c(LL1, LLbest1, LLst1))
    print(c(LL2, LLbest2, LLst2))
  }
}
