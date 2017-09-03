#####�������X�N�������ԉ��#####
library(MASS)
library(survival)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

##�f�[�^�̐ݒ�
hh <- 15000
m <- 3  #�������X�N��

####�����ϐ��̐ݒ�####
#�A���ϐ��̔���
k1 <- 4
X1 <- abs(matrix(runif(hh*k1, 0, 1.5), nrow=hh, ncol=k1))


#��l�ϐ��̔���
k2 <- 3
X2 <- matrix(0, hh, k2)
for(i in 1:k2){
  X2[, i] <- rbinom(hh, 1, runif(1, 0.3, 0.7))
}

#���l�ϐ��̔���
k3 <- 4
p.multi <- c(0.3, 0.2, 0.1, 0.4)
X3 <- t(rmultinom(hh, 1, p.multi))
X3 <- X3[, -which.min(colSums(X3))]

#�f�[�^�̌���
X <- cbind(1, X1, X2, X3)

####�����ϐ��̐ݒ�####
##�p�����[�^�̐ݒ�
Y <- matrix(0, nrow=hh, ncol=m)
alpha <- rep(0, m)
b1 <- matrix(0, nrow=m, ncol=ncol(X))
lambda <- matrix(0, nrow=hh, ncol=m)

for(j in 1:m){
  for(i in 1:10000){
    print(i)
    alpha[j] <- runif(1, 0.7, 1.6)   #�`��p�����[�^
    b1[j, ] <- c(runif(1, 0, 1.4), runif(k1, -0.6, 0.8), runif(k2, -0.6, 1.1), runif(k3-1, -0.7, 1.1))   #�Œ���ʂ̃p�����[�^
    
    #���C�u�����z����C�x���g���Ԃ𔭐�
    lambda[, j] <- exp(X %*% b1[j, ])   #���`����
    Y[, j] <- rweibull(nrow(lambda), shape=alpha[j], scale=lambda[, j])
    if(min(Y[, j]) > 0.01 & max(Y[, j]) < 200) break   
  }
}

summary(Y)
alphat <- alpha
betat <- b1

##�ł��؂�̐ݒ�
#�ł��؂�ϐ��̐ݒ�
Z <- matrix(0, nrow=hh, ncol=m)
z <- apply(Y, 1, which.min)
for(i in 1:hh) {Z[i, z[i]] <- 1 }

#�C�x���g���Ԃ����ꂼ��̃C�x���g�̍ŏ����Ԃɍ��킹��
y1 <- rowSums(Y * Z)
D <- round(cbind(Y, y1, Z, lambda), 2)


####�������X�N���f�����Ŗސ���####
##�������X�N���f���̑ΐ��ޓx
llike <- function(theta, y, X, Z, m, k){
  
  #�p�����[�^�̐ݒ�
  a <- exp(theta[1:m])
  beta <- matrix(theta[(m+1):length(theta)], nrow=k, ncol=m)

  #�ΐ��ޓx�̌v�Z
  LLi <- matrix(0, nrow=hh, ncol=m)
  for(i in 1:m){
    lambda <- exp(X %*% beta[, i])   #���`����
    LLi[, i] <- Z[, i]*(log(lambda)+log(a[i])+(a[i]-1)*log(y)) - lambda*y^a[i]   #�ΐ��ޓx���v�Z
  }
  
  LL <- sum(LLi)
  return(LL)
}

##���j���[�g���@�őΐ��ޓx���ő剻
for(i in 1:10000){
  print(i)
  #�����p�����[�^�̐ݒ�
  alpha0 <- runif(m, -0.2, 0.2)
  beta0 <- as.numeric(matrix(c(runif(m, 0.5, 1.2), runif(k1*m, -0.5, 1.0), runif((k2+k3-1)*m, -0.6, 1.0)), 
                             nrow=ncol(X), ncol=m, byrow=T))
  theta0 <- c(alpha0, beta0)

  
  #���j���[�g���@�őΐ��ޓx���ő剻
  res <- try(optim(theta0, llike, y=y1, X=X, Z=Z, m=m, k=ncol(X), method="BFGS", 
                   hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #�G���[����
}

####���ʂ̊m�F�Ɨv��#### 
##���肳�ꂽ�p�����[�^
alpha <- exp(res$par[1:3])   #�`��p�����[�^�̐���l
theta <- matrix(res$par[(m+1):length(res$par)], nrow=m, ncol=ncol(X), byrow=T)   #��A�p�����[�^�̐���l

#�^�̃p�����[�^�Ƃ̔�r
colSums(Z)
round(rbind(alpha, alphat), 3)
round(rbind(matrix(res$par[(m+1):length(res$par)], nrow=m, ncol=ncol(X), byrow=T), betat), 3)

##���v�ʂ�AIC
round(res$value, 3)   #�ő�ΐ��ޓx
round(tval <- res$par/sqrt(-diag(solve(res$hessian))), 3)   #t�l
round(AIC <- -2*res$value + 2*length(res$par), 3)   #AIC
round(BIC <- -2*res$value + log(hh)*length(res$par), 3)   #BIC


