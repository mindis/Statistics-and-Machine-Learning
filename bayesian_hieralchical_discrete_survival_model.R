#####���E�̐��ݕϐ����܂ފK�w�x�C�Y���ϗʗ��U���ԃn�U�[�h���f��#####
library(MASS)
library(nlme)
library(glmm)
library(suvival)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

set.seed(9483)

####�f�[�^�̔���####
hh <- 1000   #�T���v����
pt <- 36   #�ϑ�����
m <- 9

##ID�̐ݒ�
u.id <- rep(1:hh, rep(pt, hh))
t.id <- rep(1:pt, hh)
ID <- data.frame(no=1:(hh*pt), id=u.id, time=t.id)

####�����ϐ��̔���####
##�K�w���f���̐����ϐ�
cont1 <- 3; bin1 <- 3; multi1 <- 4
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
ZX <- cbind(1, X.cont, X.bin, X.multi)


##�̓����f���̐����ϐ��̔���
#���U�ł���UR�������Ƃɓ���
ur <- matrix(0, nrow=pt, ncol=m)

for(i in 1:pt){
  for(j in 1:1000){
    ur[i, ] <- t(rmultinom(1, 2, rep(1/m, m)))
    if(i==1){ break
    } else {
      if(max(colSums(ur[(i-1):i, ]))==1) break
    }
  }
}

#�S�f�[�^���̃f�[�^�𔭐�
UR <- matrix(as.numeric(t(UR)), nrow=pt*hh, ncol=m, byrow=T)
colnames(UR) <- c("honoka", "kotori", "umi", "rin", "hanayo", "maki", "nico", "nozomi", "eri")
UR <- UR[, -1]   #������o�[���폜


#���U�D�҂̗L���𔭐�
camp <- rbinom(pt, 1, runif(1, 0.4, 0.55))
Camp <- rep(camp, hh)


####�����ϐ��𔭐�####
##�g�����h�����𔭐�
T <- 10000
for(t1 in 1:T){
  tb <- -1.0   #�����l
  trend <- numeric()
  s <- seq(0.85, 0.2, length=t)
  for(i in 1:pt){
    r <- rnorm(5, tb, 0.05)
    sort <- sort(r)
    bi <- rbinom(1, 1, s[i])
    bb <- ifelse(bi == 1, sort[4], sort[2])
    tb <- bb
    trend <- c(trend, bb)
  }
  if(max(pnorm(trend)) < 0.25 & min(pnorm(trend)) > 0.15) break
  print(t1)
}  
plot(pnorm(trend), type="l", lwd=1, xlab="��", ylab="p")
pnorm(trend)


##�K�w���f���̃p�����[�^��ݒ�
#��A�p�����[�^��ݒ�
theta0 <- c(runif(1, -0.5, 0.5), runif(cont1, 0, 0.6), runif(bin1+multi1-1, -0.6, 0.6))
theta1 <- rbind(c(0.8, 0.2, -0.6, -0.4, 1.0, 0.6, -0.8, 0.5), matrix(runif(cont1*(m-1), 0, 0.6), nrow=cont1, ncol=m-1), 
                matrix(runif((bin1+multi1-1)*(m-1), -0.6, 0.7), nrow=bin1+multi1-1, ncol=m-1))
theta2 <- c(runif(1, 0.3, 0.8), runif(cont1, 0, 0.6), runif(bin1+multi1-1, -0.5, 0.8))
theta3 <- c(runif(1, -0.9, -0.4), runif(cont1, 0, 0.6), runif(bin1+multi1-1, -0.5, 0.8))
theta4 <- c(runif(1, 0.2, 0.5), runif(cont1, 0, 0.4), runif(bin1+multi1-1, -0.4, 0.5))
theta_T <- cbind(theta0, theta1, theta2, theta3, theta4)

#���U�����U�s���ݒ�
Cov0 <- diag(runif(ncol(theta_T), 0.4, 0.8))

#���ϗʐ��K���z���p�����[�^�𔭐�
BETA0 <- ZX %*% theta_T + mvrnorm(hh, rep(0, ncol(theta_T)), Cov0)
BETA_T <- BETA0


##�����ϐ��𔭐������Ȃ��玞�Ԉˑ��ϐ��𔭐�������