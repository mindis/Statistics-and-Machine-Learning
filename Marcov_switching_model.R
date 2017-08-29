#####�}���R�t�؂�ւ����f��#####
library(MASS)
library(MSwM) 
library(reshape2)
library(gtools)
library(dplyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
n <- 1000   #�T���v����
k1 <- 3   #�؂�ւ���
k2 <- 2   #�ϑ��m���̃p�����[�^��


##�����m���̒�`
Pf <- c(0.4, 0.3, 0.3)

##���ڍs��̒�`
pr1 <- c(0.1, 0.7, 0.2)
pr2 <- c(0.2, 0.1, 0.7)
pr3 <- c(0.7, 0.2, 0.1)
Pr <- rbind(pr1, pr2, pr3)

##�ϑ��m���̒�`
P <- matrix(0, nrow=k1, ncol=k2)
for(i in 1:k1){
  alpha <- runif(1, 0.4, 1)
  P[i, ] <- rdirichlet(1, rep(alpha, k2))
}

P <-rbind(c(0.9, 0.1), c(0.6, 0.4), c(0.1, 0.9))

##�����ϐ��̔���
Z <- matrix(0, nrow=n, ncol=k1)
Y <- matrix(0, nrow=n, ncol=k2)

#���ݕϐ��̏����l
Z[1, ] <- rmultinom(1, 1, Pf)
Y[1, ] <- rmultinom(1, 1, P[which.max(Z[1, ]), ])

#2��ڈȍ~�̉����ϐ��𒀎��I�ɔ���������
for(i in 2:n){
  Z[i, ] <- rmultinom(1, 1, Pr[which.max(Z[i-1, ]), ])
  Y[i, ] <- rmultinom(1, 1, P[which.max(Z[i, ]), ])
}


####EM�A���S���Y���Ń}���R�t�؂�ւ����f���𐄒�####
##�����l�̐ݒ�
#�����m���̐ݒ�
rho <- rep(0.25, k1)   

#�}���R�t���ڍs��̏����l�̐ݒ�
A <- matrix(0, nrow=k1, ncol=k1)
for(i in 1:k1){
  p_rand <- runif(k1, 0.1, 1)
  A[i, ] <- p_rand / sum(p_rand)
}

#�ϑ����f���̃p�����[�^
B <- matrix(0, nrow=k1, ncol=k2)
for(i in 1:k1){
  p_rand <- runif(k2, 0.1, 1)
  B[i, ] <- p_rand / sum(p_rand)
}
y <- Y %*% 1:k2

#�f�[�^�̏���
y_delta <- array(0, dim=c(n, k1, k2))
for(i in 1:k2){
  y_delta[, , i] <- matrix(Y[, i], nrow=n, ncol=k1)
}

#�p�����[�^�̊i�[�p�z��
alpha <- matrix(0, nrow=n, ncol=k1)
alpha_s <- matrix(0, nrow=n, ncol=k1)
alpha_mu <- rep(0, n)
B_vec <- matrix(0, nrow=n, ncol=k1)
beta <- matrix(0, nrow=n, ncol=k1)
beta_s <- matrix(0, nrow=n, ncol=k1)

LL1 <- sum(log(apply(Y, 1, function(x) dmultinom(x, 1, B[1, ]))))   #�ΐ��ޓx�̏����l
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 0.01

L <- c()
for(i in 1:n){
  L <- c(L, dmultinom(Y[i, ], 1, P[which.max(Z[i, ]), ]))
}
LT <- sum(log(L))


####EM�A���S���Y��(�o�E���E�F���`�A���S���Y��)�ŉB��}���R�t���f���𐄒�####
while(abs(dl) >= tol){
  A <- Pr
  B <- P
  rho <- Pf
  
  ##�O�����A���S���Y����alpha�𐄒�
  B_vec[1, ] <- B[, y[1]]
  alpha[1, ] <- rho * B_vec[1, ]
  alpha_mu[1] <- 1 / sum(alpha[1, ])
  alpha_s[1, ] <- alpha_mu[1] * alpha[1, ] 
  
  for(i in 2:n){
    B_vec[i, ] <- B[, y[i]]
    alpha[i, ] <- alpha_s[i-1, ] %*% A * B_vec[i, ]
    alpha_mu[i] <- 1 / sum(alpha[i, ])
    alpha_s[i, ] <- alpha[i, ] * alpha_mu[i]
  }
  
  ##�������A���S���Y����beta�𐄒�
  beta[n, ] <- 1
  beta_s[n, ] <- alpha_mu[n]
  
  for(i in n:2){
    beta[i-1, ] <- A %*% (B_vec[i, ] * beta_s[i, ])
    beta_s[i-1, ] <- beta[i-1, ] * alpha_mu[i-1] 
  }
  
  ##�p�����[�^���X�V
  #���ڊm��A���X�V
  for(i in 1:k1){
    A_vec <- matrix(A[i, ], nrow=n-1, ncol=k1, byrow=T)
    a11 <- matrix(alpha_s[1:(n-1), i], n-1, k1) * A_vec * B_vec[2:n, ] * beta_s[2:n, ]
    a12 <- matrix(alpha_s[1:(n-1), i], n-1, k1) * matrix(beta_s[1:(n-1), i], n-1, k1) / alpha_mu[1:(n-1)]
    a <- colSums(a11)/colSums(a12)
    A[i, ] <- a
  }
  
  
  #�ϑ��f�[�^�̃p�����[�^B���X�V
  for(j in 1:k2){
    B[, j] <- colSums(y_delta[, , j] * alpha_s * beta_s / alpha_mu) / colSums(alpha_s * beta_s / alpha_mu)
  }
  
  #�����m��rho���X�V
  rho <- alpha_s[1, ] * beta_s[1, ] / alpha_mu[1]
  
  ##�ΐ��ޓx���v�Z
  LL <- sum(log(1/alpha_mu))
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

round(cbind(B, P), 3)
round(cbind(A, Pr), 3)
LT

