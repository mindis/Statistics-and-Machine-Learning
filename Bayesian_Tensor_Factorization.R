#####Bayesian Tensor Factorization#####
library(MASS)
library(matrixStats)
library(FAdist)
library(NMF)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 5000   #���[�U�[��
item <- 1000   #�A�C�e����
time <- 12   #�ϑ��@�֐�
N0 <- hh*item*time
k <- 10   #��ꐔ
vec <- rep(1, k)

#ID��ݒ�
user_id0 <- rep(1:hh, item*time)
item_id0 <- rep(rep(1:item, rep(hh, item)), time)
time_id0 <- rep(1:time, rep(hh*item, time))
id <- cbind(user_id0, item_id0, time_id0)

##�e���\�������̒�`�Ɋ�Â��f�[�^�𐶐�
for(rp in 1:1000){
  print(rp)
  #���ϗʐ��K���z�̃p�����[�^��ݒ�
  mu01 <- rep(1.2, k); mu02 <- rep(0.9, k); mu03 <- rep(0.5, k)
  tau01 <- tau02 <- tau03 <- diag(0.2, k)
  sigma <- sigmat <- 0.5

  #���ϗʐ��K���z��������s��𐶐�
  W0 <- WT0 <- mvrnorm(hh, mu01, tau01)
  H0 <- HT0 <- mvrnorm(item, mu02, tau02)
  C0 <- CT0 <- mvrnorm(time, mu03, tau03)
  
  #���ϗʐ��K���z����]�_�X�R�A�𐶐�
  WHC0 <- array(0, dim=c(hh, item, time))
  for(j in 1:k){
    WHC0 <- WHC0 + W0[, j] %o% t(H0)[j, ] %o% t(C0)[j, ]
  }
  #whc0 <- as.numeric((W0[user_id0, ] * t(H0)[item_id0, ] * t(C0)[time_id0, ]) %*% vec)   #������ł�ok
  whc0 <- as.numeric(WHC0)   #�e���\�����x�N�g���ɕϊ�
  y_vec0 <- whc0 + rnorm(hh*item*time, 0, sigma)   #�덷�𐶐�
  
  #��������
  if(mean(y_vec0) < 5.5 & mean(y_vec0) > 4.5 & sd(y_vec0) > 1.5 & sd(y_vec0) < 2.0) break
}

#�����ϐ���1�`10�ɕϊ�����
y0 <- round(y_vec0)
y0[y0 > 10] <- 10; y0[y0 < 1] <- 1


##�����x�N�g���𐶐�
#�����m���𐶐�
user_prob <- rbeta(hh, 10, 50)
item_prob <- rbeta(item, 15, 55)
time_prob <- rbeta(time, 60, 140)
prob <- user_prob[user_id0]*item_id0[item_id0]*time_prob[time_id0]

#�x���k�[�C���z���猇���x�N�g���𐶐�
z_vec <- rbinom(N0, 1, prob)
N <- sum(z_vec)
y <- y0[z_vec==1]; y_vec <- y_vec0[z_vec==1]; whc <- whc0[z_vec==1]
hist(y_vec, breaks=25, col="grey", main="���ݓI�ȃX�R�A���z", xlab="�X�R�A")
hist(y, breaks=25, col="grey", main="�����ł̃X�R�A���z", xlab="�X�R�A")

#�����x�N�g������id���č\��
user_id <- user_id0[z_vec==1]
item_id <- item_id0[z_vec==1]
time_id <- time_id0[z_vec==1]


####�}���R�t�A�������e�J�����@�Ńe���\�������𐄒�####
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
W <- WT0
H <- HT0
C <- CT0
sigma <- sigmat

##�����l�̐ݒ�
sigma <- sd(y)
W <- mvrnorm(hh, rep(0.8, k), 0.2 * diag(k))
H <- mvrnorm(item, rep(0.8, k), 0.2 * diag(k))
C <- mvrnorm(time, rep(0.8, k), 0.2 * diag(k))

##�p�����[�^�̊i�[�p�z��
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(item, k, R/keep))
C_array <- array(0, dim=c(time, k, R/keep))
SIGMA <- rep(0, R/keep)


##�C���f�b�N�X��ݒ�
user_index <- item_index <- time_index <- list()
ui_id <- ut_id <- iu_id <- it_id <- tu_id <- ti_id <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
  ui_id[[i]] <- item_id[user_index[[i]]]
  ut_id[[i]] <- time_id[user_index[[i]]]
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
  iu_id[[j]] <- user_id[item_index[[j]]]
  it_id[[j]] <- time_id[item_index[[j]]]
}
for(j in 1:time){
  time_index[[j]] <- which(time_id==j)
  tu_id[[j]] <- user_id[time_index[[j]]]
  ti_id[[j]] <- item_id[time_index[[j]]]
}
const1 <- hh / 1.5  #���K���萔
const2 <- item / 1.5

##�ΐ��ޓx�̊�l
LLst <- sum(dnorm(y, mean(y), sd(y), log=TRUE))


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##���[�U�[�����s����T���v�����O
  for(i in 1:hh){
    #�����x�N�g���̃p�����[�^
    X <- H[ui_id[[i]], ] * C[ut_id[[i]], ]
    inv_XXV <- solve(t(X) %*% X +inv_tau)
    beta_mu <- inv_XXV %*% t(X) %*% y[user_index[[i]]]   #���ϗʐ��K���z�̕��σx�N�g��
    cov <- sigma^2 * inv_XXV
    
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    w <- mvrnorm(1, beta_mu, cov)
    W[i, ] <- w
  }
  W <- W / matrix(colSums(W), nrow=hh, ncol=k, byrow=T) * hh
  
  ##�A�C�e�������s����T���v�����O
  for(j in 1:item){
    #�����x�N�g���̃p�����[�^
    X <- W[iu_id[[j]], ] * C[it_id[[j]], ]
    inv_XXV <- solve(t(X) %*% X +inv_tau)
    beta_mu <- inv_XXV %*% t(X) %*% y[item_index[[j]]]   #���ϗʐ��K���z�̕��σx�N�g��
    cov <- sigma^2 * inv_XXV
    
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    h <- mvrnorm(1, beta_mu, cov)
    H[j, ] <- h
  }
  H <- H / matrix(colSums(H), nrow=item, ncol=k, byrow=T) * item
  
  ##���Ԃ̓����s����T���v�����O
  for(j in 1:time){
    #�����x�N�g���̃p�����[�^
    X <- W[tu_id[[j]], ] * H[ti_id[[j]], ]
    inv_XXV <- solve(t(X) %*% X +inv_tau)
    beta_mu <- inv_XXV %*% t(X) %*% y[time_index[[j]]]   #���ϗʐ��K���z�̕��σx�N�g��
    cov <- sigma^2 * inv_XXV
    
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    c <- mvrnorm(1, beta_mu, cov)
    C[j, ] <- c
  }
  
  
  ##���f���̕W���΍����T���v�����O
  #���f���̌덷�𐄒�
  mu <- (W[user_id, ] * H[item_id, ] * C[time_id, ]) %*% vec   #���f���̕���
  er <- y - mu   #���f���̌덷
  
  #�t�K���}���z�̃p�����[�^
  s <- as.numeric(t(er) %*% er) + s0
  v <- N + v0
  
  #�t�K���}���z����W���΍����T���v�����O
  sigma <- sqrt(1/rgamma(1, v/2, s/2))
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    W_array[, , mkeep] <- W
    H_array[, , mkeep] <- H
    C_array[, , mkeep] <- C
    SIGMA[mkeep] <- sigma
  }
  
  #�ΐ��ޓx���v�Z
  if(rp%%disp==0){
    LLi <- dnorm(y, mu, sigma, log=TRUE)
    LL <- sum(LLi)
    
    #�T���v�����O���ʂ̕\��
    print(rp)
    print(c(LL, LLst))
    print(round(c(sigma, sigmat), 3))
  }
}

matplot(t(W_array[1, , ]), type="l")
matplot(t(H_array[1, , ]), type="l")
matplot(t(C_array[1, , ]), type="l")
C_array
W
H
C
C_array
round(W_array, 3)

index_na <- which(z_vec==0)
sum((y0[index_na] - rowSums(W[user_id0[index_na], ] * H[item_id0[index_na], ] * C[time_id0[index_na], ]))^2)
sd((y0[index_na] - rowSums(WT0[user_id0[index_na], ] * HT0[item_id0[index_na], ] * CT0[time_id0[index_na], ]))^2)
