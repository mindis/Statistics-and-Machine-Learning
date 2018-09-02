#####Gibbs Sampler for logistic regression#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(FAdist)
library(bayesm)
library(extraDistr)
library(condMVNorm)
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
N <- 2000
beta <- 0.5   
sigma <- 1.0
x <- rep(1, N)

##�����ϐ��̐���
prob <- exp(beta) / (1 + exp(beta))
y <- rbinom(N, 1, prob)
index_y1 <- which(y==1); index_y0 <- which(y==0)
n1 <- length(index_y1); n0 <- length(index_y0)

####�M�u�X�T���v���[�Ńp�����[�^���T���v�����O####
##MCMC�̐ݒ�
R <- 1000
beta <- 0   #�����l
BETA <- rep(0, R)   #�p�����[�^�̊i�[�p�z��

##�p�����[�^���T���v�����O
for(rp in 1:100000){
  mu <- x * beta
  prob <- exp(mu) / (1 + exp(mu))
  
  r <- rep(0, N)
  r[index_y1] <- runif(n1, 0, prob[index_y1])
  r[index_y0] <- runif(n0, prob[index_y0], 1)
  logit <- log(r / (1-r))
  
  a <- max(logit[index_y1])
  b <- min(logit[index_y0])
  beta <- extraDistr::rtnorm(1, 0, 1, a, b)
  BETA[rp] <- beta
}
beta
plot(1:rp, BETA, type="l")


mean(u[index_y1])
mean(u[index_y0])
hist(u)


beta
plot(1:1000, BETA, type="l")

mean(log(r/(1-r)))



