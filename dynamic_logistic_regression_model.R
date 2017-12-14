#####���I���U�I�����f��#####
library(MASS)
library(Matrix)
library(matrixStats)
library(RcppSMC)
library(SMC)
library(dml)
library(KFAS)
library(extraDistr)
library(reshape2)
library(dplyr)

####�f�[�^�̔���####
#set.seed(54390)

##�f�[�^�̐ݒ�
n <- 1000   #�ϑ�����
time <- 1:n   #�ϑ�����id
k <- 4   #�����ϐ���


####�����ϐ��𔭐�####
#���i�𐶐�
PRICE <- runif(n, 0, 1.0)   

#�v�����[�V�����L���𐶐�
DISP <- rbinom(n, 1, 0.4)   #���ʒ�L��
AD <- rbinom(n, 1, 0.4)   #�`���V�f�ڗL��

#�f�[�^������
Data <- cbind(PRICE, DISP, AD)


####�p�����[�^�𐶐����A�����ϐ�������####
##�g�����h�����̐���
b0 <- -0.8   #�g�����h�̏����l
v0 <- 0.015   #�V�X�e�����f���̕��U

#���Ԃ��Ƃɒ����I�Ƀg�����h�����𐶐�
trend <- rep(0, n)
trend[1] <- b0
s <- seq(0.7, 0.3, length=n-1)   
for(i in 2:n){
  diff <- rnorm(5, 0, v0)   #�ω��̌��𐶐�
  sortlist <- sort(diff)   #�����ɕ��ёւ���
  bi <- rbinom(1, 1, s[i-1])   #�ω��̎d��������
  trend[i] <- trend[i-1] + bi*sortlist[4] + (1-bi)*sortlist[2]
}
plot(1:n, trend, type="l", xlab="�ϑ�����")
summary(trend); trend0 <- trend

##�����ϐ��̓��I�p�����[�^�𐶐�
#�����l��ݒ�
beta1 <- beta2 <- beta3 <- rep(0, n)
beta1[1] <- -0.9   #���i�̏����l
beta2[1] <- 0.8   #���ʒ�̏����l
beta3[1] <- 0.7   #�`���V�f�ڗL���̏����l 

#�V�X�e�����f���̕��U
v1 <- 0.015   
v2 <- 0.015
v3 <- 0.015

#���Ԃ��Ƃɒ����I�ɓ��I�p�����[�^�𐶐�
s1 <- seq(0.6, 0.4, length=n-1)
s2 <- seq(0.40, 0.7, length=n-1)
s3 <- seq(0.45, 0.6, length=n-1)
for(i in 2:n){
  diff1 <- rnorm(5, 0, v1); diff2 <- rnorm(5, 0, v2); diff3 <- rnorm(5, 0, v3)
  sortlist1 <- sort(diff1); sortlist2 <- sort(diff2); sortlist3 <- sort(diff3)
  bi1 <- rbinom(1, 1, s1[i-1]); bi2 <- rbinom(1, 1, s2[i-1]); bi3 <- rbinom(1, 1, s3[i-1])
  beta1[i] <- beta1[i-1] + bi1*sortlist1[2] + (1-bi1)*sortlist1[4]
  beta2[i] <- beta2[i-1] + bi2*sortlist2[4] + (1-bi2)*sortlist2[2]
  beta3[i] <- beta3[i-1] + bi3*sortlist3[4] + (1-bi3)*sortlist3[2]
}
plot(1:n, beta1, type="l", xlab="�ϑ�����")
plot(1:n, beta2, type="l", xlab="�ϑ�����")
plot(1:n, beta3, type="l", xlab="�ϑ�����")

#�p�����[�^������
beta0 <- cbind(beta1, beta2, beta3)


##���I���W�X�e�B�b�N��A���f�����牞���ϐ��𐶐�
#���W�b�g�Ɗm���̒�`
logit <- trend + rowSums(Data * beta0)   #���W�b�g
Pr0 <- exp(logit)/(1+exp(logit))   #�m��

#�x���k�[�C���z���牞���ϐ��𐶐�
y <- rbinom(n, 1, Pr0)

#���������������ϐ����v���b�g
plot(1:n, y, xlab="�ϑ�����", ylab="�w���L��", main="�w���L���ƍw���m���̊֘A")
par(new=T)
plot(1:n, Pr0, xlim=c(0, n), ylim=c(0, 1), ylab="", xlab="", type="b", pch=4)


####���q�t�B���^�œ��I���W�X�e�B�b�N��A���f���𐄒�####
##���W�X�e�B�b�N��A���f���̑ΐ��ޓx�֐�(�����l���萔)
fr <- function(beta, X, y){
  
  #���W�b�g�Ɗm�����`
  logit <- X %*% beta
  Pr <- exp(logit)/(1+exp(logit))
  
  #�ΐ��ޓx�֐��̘a
  Li <- y*log(Pr) + (1-y)*log(1-Pr)
  LL <- sum(Li)
  return(LL)
}

####���q�t�B���^�œ��I�p�����[�^�𐄒�####
##���q�t�B���^�̐ݒ�
s <- 10000   #���q��
tau <- diag(0.05^2, k)   #�V�X�e�����U
LL <- rep(0, n)
BETA <- array(0, dim=c(s, k, n))

##�V�X�e�����f���̃p�����[�^�̍X�V
betan <- mvrnorm(s, rep(0, k), diag(0.25, k))

##�ϑ����f���̖ޓx��]��
#���W�b�g�Ɗm�����v�Z
logit <- betan[, 1] + rowSums(matrix(Data[1, ], nrow=s, ncol=k-1, byrow=T) * betan[, -1])
Pr <- exp(logit)/(1+exp(logit))

#�ޓx��]��
Li <- Pr^y[1] * (1-Pr)^(1-y[1])   #���q���Ƃ̖ޓx
LL[1] <- sum(Li)   #�ޓx�̘a

#�ޓx�̕��S���ɉ����ăp�����[�^�����T���v�����O
w <- Li/sum(Li) 
index <- as.numeric(rmnom(1, s, w))
resample <- rep(1:s, index)
BETA[, , 1] <- betan[resample, ]   #���T���v�����O���ꂽ�p�����[�^


##2���ڈȍ~�𗱎q�t�B���^�Œ����I�ɍX�V
for(i in 2:n){
  ##�V�X�e�����f���̃p�����[�^�̍X�V
  betan <- BETA[, , i-1] + mvrnorm(s, rep(0, k), tau)
  
  ##�ϑ����f���̖ޓx��]��
  #���W�b�g�Ɗm�����v�Z
  logit <- betan[, 1] + rowSums(matrix(Data[i, ], nrow=s, ncol=k-1, byrow=T) * betan[, -1])
  Pr <- exp(logit)/(1+exp(logit))
  
  #�ޓx��]��
  Li <- Pr^y[i] * (1-Pr)^(1-y[i])   #���q���Ƃ̖ޓx
  LL[i] <- sum(Li)   #�ޓx�̘a
  
  #�ޓx�̕��S���ɉ����ăp�����[�^�����T���v�����O
  w <- Li/sum(Li) 
  index <- as.numeric(rmnom(1, s, w))
  resample <- rep(1:s, index)
  BETA[, , i] <- betan[resample, ]   #�p�����[�^�����T���v�����O
}
#�ΐ��ޓx�̘a
LLs <- sum(log(LL)) - n*log(s)


##�Œ胉�O�������ɂ��ŏI�I�ȃp�����[�^���m��
L <- 20   #���O��
LAG <- array(0, dim=c(s, k, L))
THETA <- BETA

for(i in L:n){
  print(i)
  
  lag <- i-L+1
  LAG <- THETA[, , lag:i]
  betan <- LAG[, , L-1] + mvrnorm(s, rep(0, k), tau)
  
  ##�ϑ����f���̖ޓx��]��
  #���W�b�g�Ɗm�����v�Z
  logit <- betan[, 1] + rowSums(matrix(Data[i, ], nrow=s, ncol=k-1, byrow=T) * betan[, -1])
  Pr <- exp(logit)/(1+exp(logit))
  
  #�ޓx��]��
  Li <- Pr^y[i] * (1-Pr)^(1-y[i])   #���q���Ƃ̖ޓx
  
  #�ޓx�̕��S���ɉ����ăp�����[�^�����T���v�����O
  w <- Li/sum(Li) 
  index <- as.numeric(rmnom(1, s, w))
  resample <- rep(1:s, index)
  THETA[, , lag:i] <- LAG[resample, , ]   #�p�����[�^�����T���v�����O
}


####���茋�ʂ̊m�F�Ɖ���####
##���㕽�ς��m�F
theta1 <- t(apply(BETA, c(2, 3), mean))
theta2 <- t(apply(THETA, c(2, 3), mean))
round(cbind(theta1, theta2, trend, beta0), 3)

#�t�B���^�����O�A�������A�^�l������
matplot(theta1, type="l", xlab="�ϑ�����", ylab="�p�����[�^")
matplot(theta2, type="l", xlab="�ϑ�����", ylab="�p�����[�^")
matplot(cbind(trend, beta1, beta2, beta3), type="l", xlab="�ϑ�����", ylab="�p�����[�^")

##�m���Ƒΐ��ޓx���Čv�Z
logit <- trend + rowSums(Data * beta0)
Pr1 <- exp(logit)/(1+exp(logit))
Li1 <- Pr1^y * (1-Pr1)^(1-y)   #���q���Ƃ̖ޓx
LL1 <- sum(log(Li1))

logit <- theta1[, 1] + rowSums(Data * theta1[, -1])
Pr2 <- exp(logit)/(1+exp(logit))
Li2 <- Pr2^y * (1-Pr2)^(1-y)   #���q���Ƃ̖ޓx
LL2 <- sum(log(Li2))

logit <- theta2[, 1] + rowSums(Data * theta2[, -1])
Pr3 <- exp(logit)/(1+exp(logit))
Li3 <- Pr3^y * (1-Pr3)^(1-y)   #���q���Ƃ̖ޓx
LL3 <- sum(log(Li3))
c(LL1, LL2, LL3)
round(cbind(y, Pr1, Pr2, Pr3), 3)
c(LL1, LL2, LL3)

##����M�p��Ԃ𐄒�
theta_lower <- matrix(0, nrow=n, ncol=k)
theta_upper <- matrix(0, nrow=n, ncol=k)
for(i in 1:n){
  theta_lower[i, ] <- apply(THETA[, , i], 2, function(x) quantile(x, 0.05))
  theta_upper[i, ] <- apply(THETA[, , i], 2, function(x) quantile(x, 0.95))
}


##�g�����h���������n��Ƀv���b�g
plot(1:n, y, xlab="�ϑ�����", ylab="�w���L��", main="�w���L���ƍw���m���̊֘A")
par(new=T)
plot(1:n, Pr0, xlim=c(0, n), ylim=c(0, 1), ylab="", xlab="", type="b", pch=4)
par(new=T)
plot(1:n, exp(theta2[, 1])/(1+exp(theta2[, 1])), xlim=c(0, n), ylim=c(0, 1), 
     ylab="", xlab="", type="l", lwd=2, col=2)
par(new=T)
plot(1:n, exp(theta_lower[, 1])/(1+exp(theta_lower[, 1])), xlim=c(0, n), ylim=c(0, 1), 
     ylab="", xlab="", type="l", lwd=2, col=3)
par(new=T)
plot(1:n, exp(theta_upper[, 1])/(1+exp(theta_upper[, 1])), xlim=c(0, n), ylim=c(0, 1), 
     ylab="", xlab="", type="l", lwd=2, col=3)
par(new=T)
plot(1:n, exp(theta1[, 1])/(1+exp(theta1[, 1])), xlim=c(0, n), ylim=c(0, 1), 
     ylab="", xlab="", type="l", lwd=2, col=4)
par(new=T)
plot(1:n, exp(trend)/(1+exp(trend)), xlim=c(0, n), ylim=c(0, 1), ylab="", xlab="", 
     type="l", lwd=2, col=5)
