#####�x�C�W�A���������������W�b�g���f��#####
library(MASS)
library(matrixStats)
library(flexmix)
library(glmnet)
library(FAdist)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
hh <- 3000
select <- 6
k <- 300   #�����ϐ���

##�����ϐ��̔���
freq <- rpois(hh, 200)   #�|�A�\�����z����p�x�𔭐�
p <- rdirichlet(hh, runif(k, 0.2, 1.0))   #�f�B���N�����z����o���m���𔭐�
X <- t(apply(cbind(freq, p), 1, function(x) rmultinom(1, x[1], x[-1])))   #�������z��������ϐ��𔭐�


#�����ϐ����x�N�g����
#ID�̃x�N�g����
u.id <- rep(1:hh, rep(select, hh))
i.id <- rep(1:select, hh)
ID <- data.frame(no=1:(hh*select), u.id=u.id, i.id=i.id)

#�ؕЂ̃x�N�g����
BP <- matrix(diag(select), nrow=hh*select, ncol=select, byrow=T)[, -select]

#�����ϐ��̃x�N�g����
X_vec <- matrix(0, nrow=hh*select, ncol=ncol(X)*(select-1))

for(i in 1:hh){
  x_diag0 <- c()
  for(j in 1:ncol(X)){
    x_diag0 <- cbind(x_diag0, diag(X[i, j], select-1))
  }
  X_vec[ID$u.id==i, ] <- rbind(x_diag0, 0)
}
XM_vec <- cbind(BP, X_vec)


##�����ϐ��̔���
for(i in 1:1000){
  print(i)
  
  #�p�����[�^�̐ݒ�
  b00 <- runif(select-1, -1.0, 1.0)
  bo1 <- matrix(runif((select-1)*k, -1, 1), nrow=k, ncol=select-1) * matrix(rbinom(k, 1, 0.15), nrow=k, ncol=select-1)
  b01 <- as.numeric(t(ifelse(abs(bo1) > 0.25, bo1, 0)))
  b0 <- c(b00, b01)
  
  #���W�b�g�Ɗm���̌v�Z
  logit <- matrix(XM_vec %*% b0, nrow=hh, ncol=select, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #�������z���牞���ϐ��𔭐�
  y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
  if(min(colMeans(y)) > 0.1 & max(colMeans(y)) < 0.4) break
}

####�}���R�t�A�������e�J�����@��L1�������������W�b�g���f���𐄒�####
fr <- function(theta, lambda, y, X, hh, select){
  
  #���W�b�g�Ɗm���̌v�Z
  logit <- matrix(X %*% theta, nrow=hh, ncol=select, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #�����t���ΐ��ޓx���`
  LLi <- rowSums(y*log(Pr)) - lambda*sum(abs(theta))
  LL <- sum(LLi)
  return(LL)
}

##�ΐ��ޓx���ő剻
theta <- runif((select-1)*(k+1), -0.2, 0.2)
lambda <- 0.001
res <- try(optim(theta[1:15], fr, gr=NULL, lambda, y, XM_vec[, 1:15], hh, select, method="BFGS", hessian=FALSE, 
                 control=list(fnscale=-1, trace=TRUE)), silent=TRUE)
res$par
b0
fr(b0, lambda, y, XM_vec, hh, select)

