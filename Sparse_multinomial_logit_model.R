####�X�p�[�X�������W�b�g���f��####
library(MASS)
library(ncvreg)
library(glmnet)
library(lars)
library(extraDistr)
library(Matrix)
library(matrixStats)
library(reshape2)
library(plyr)
library(dplyr)

####�f�[�^�̐���####
##�f�[�^�̐ݒ�
hh1 <- 10000   #�w�K�p�f�[�^�̃T���v����
hh2 <- 5000   #���ؗp�f�[�^�̃T���v����
hh <- hh1 + hh2
select <- 10   #�I������
k1 <- 300   #�A���ϐ��̐����ϐ�
k2 <- 200   #���U�ϐ��̐����ϐ�
k <- k1 + k2

##�����ϐ��̐���
#�A���ϐ��̐���
X1 <- mvrnorm(hh, rep(0, k1), diag(1, k1))

#���U�ϐ��̐���
X2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.75)
  X2[, j] <- rbinom(hh, 1, pr)
}

#�f�[�^�̌���
X <- cbind(1, X1, X2)
X1 <- X[1:hh1, ]; X2 <- X[-(1:hh1), ]


##�����ϐ��̐���
#�p�����[�^��ݒ�
theta0 <- t(mvrnorm(select-1, rep(0, k+1), diag(0.05, k+1)))
theta <- thetat <- cbind(ifelse(abs(theta0) <= 0.2, 0, theta0), 0)   #�p�����[�^��0�̃V�������N

#���W�b�g�ƑI���m����ݒ�
logit <- X %*% theta   #���W�b�g
Pr <- exp(logit) / rowSums(exp(logit))   #�I���m��

#�������z����I�����ʂ𐶐�
y <- rmnom(hh, 1, Pr)
y1 <- y[1:hh1, ]; y2 <- y[-(1:hh1), ]
colSums(y)


####����������W�~���@�ɂ��X�p�[�X�������W�b�g���f���𐄒�####
##�������W�b�g���f���̑ΐ��ޓx
loglike <- function(y, X, theta, select){
  #���W�b�g�Ɗm�����`
  logit <- X %*% theta
  Pr <- exp(logit) / rowSums(exp(logit)) 
  
  #�ΐ��ޓx���v�Z
  LL <- sum(log(rowSums(y * Pr)))
  return(LL)
} 

##���W�~���@�̐ݒ�
#�������p�����[�^��ݒ�
lambda <- seq(0.001, 0.01, length=20)  #�������p�����[�^
n_lambda <- hh1*lambda 
X1_sq <- X1^2

#�ΐ��ޓx�̃x�X�g
LLbest <- loglike(y2, X2, thetat, select)

#�p�����[�^�̊i�[�p�z��
LLtest <- c()
THETA <- array(0, dim=c(k+1, select, length(lambda)))


####�X�p�[�X�������W�b�g���f���Ő������p�����[�^���œK��####
for(rp in 1:length(lambda)){
  print(lambda[rp])
  
  ##�A���S���Y���̐ݒ�
  #�p�����[�^�̏����l
  theta <- t(matrix(0, nrow=select, ncol=k+1))
  LL <- loglike(y1, X1, theta, select)
  
  #�A���S���Y���̍X�V�X�e�[�^�X
  LLs <- LL
  iter <- 1
  dl <- -100   #�ΐ��ޓx�̍��̏����l
  tol <- 2.5
  LL1 <- LL   #�ΐ��ޓx�̏����l
  
  ##���z�~���@�Ńp�����[�^���X�V
  while(abs(dl) >= tol){
    
    ##���W�b�g�ƑI���m�����`
    logit <- X1 %*% theta
    Pr <- exp(logit) / rowSums(exp(logit))
    
    ##�p�����[�^���X�V
    #�����W���̃p�����[�^���X�V
    w <- (Pr * (1-Pr))[, -select]
    z <- (logit[, -select] + (y1 - Pr)[, -select] / w)
    
    #�ؕЂ̍X�V
    theta[1, -select] <- colSums(w * (z - logit[, -select])) / colSums(w)
    
    #��A�W���̍X�V
    for(j in 2:(k+1)){
      S <- colSums(w * X1[, j] * (z - X1[, -j] %*% theta[-j, -select]))
      theta[j, -select] <- (S * (abs(S) > n_lambda[rp])) / colSums(w * X1_sq[, j])
    }
    
    ##�A���S���Y�����X�V
    LL <- loglike(y1, X1, theta, select)   #�ΐ��ޓx
    iter <- iter+1
    dl <- LL1 - LL
    LL1 <- LL
    LLs <- c(LLs, LL)
    print(LL)
  }
  
  ##�p�����[�^���i�[
  #�e�X�g�f�[�^�ɑ΂���ΐ��ޓx���X�V
  LLc <- loglike(y2, X2, theta, select)
  
  #�p�����[�^���i�[
  LLtest <- c(LLtest, LLc)
  THETA[, , rp] <- theta
  print(c(LLc, LLbest))
}

####�x�X�g�ȃp�����[�^�ŗ\���덷�����߂�####
##�x�X�g�Ȑ������p�����[�^
lambda_best <- which.max(LLtest)
theta <- THETA[, , lambda_best]   #�x�X�g�ȃp�����[�^
LLtest[lambda_best]   #�x�X�g�ȑΐ��ޓx

##�\���덷���v�Z
round(cbind(theta, thetat), 2)   #�^�l�Ɛ���p�����[�^�̔�r

#�e�X�g�f�[�^�ɑ΂��郍�W�b�g�Ɨ\���m�����v�Z
logit <- X2 %*% theta   #���W�b�g
Pr <- exp(logit) / rowSums(exp(logit))   #�\���m��

#�\�����ʂƊϑ��l�̐�����
pred <- apply(Pr, 1, which.max)   #�\������
round(pred_class <- table(pred, y2 %*% 1:select), 3)   
round(pred_rate <- pred_table / rowSums(pred_table), 3)   
round(sum(diag(pred_class / sum(pred_class))), 3)
