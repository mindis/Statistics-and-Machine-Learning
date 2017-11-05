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
hh <- 6000
select <- 10
k <- 30   #�����ϐ���

##�����ϐ��̔���
freq <- rpois(hh, 200)   #�|�A�\�����z����p�x�𔭐�
p <- rdirichlet(hh, runif(k, 0.2, 1.0))   #�f�B���N�����z����o���m���𔭐�
X <- scale(t(apply(cbind(freq, p), 1, function(x) rmultinom(1, x[1], x[-1]))))   #�������z��������ϐ��𔭐�


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
for(i in 1:5000){
  print(i)
  
  #�p�����[�^�̐ݒ�
  b00 <- runif(select-1, -0.7, 0.7)
  bo1 <- matrix(runif((select-1)*k, 0, 0.9), nrow=k, ncol=select-1) * matrix(rbinom(k, 1, 0.3), nrow=k, ncol=select-1)
  b01 <- as.numeric(t(ifelse(abs(bo1) > 0.2, bo1, 0)))
  b0 <- c(b00, b01)
  
  #���W�b�g�Ɗm���̌v�Z
  logit <- matrix(XM_vec %*% b0, nrow=hh, ncol=select, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #�������z���牞���ϐ��𔭐�
  y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
  if(min(colMeans(y)) > 0.05 & max(colMeans(y)) < 0.5) break
}
round(b0, 3)
colMeans(y)

####��������W�œK���ɂ��L1�������������W�b�g���f���𐄒�####
loglike <- function(b1, b2, lambda, y, X1, X2, hh, select){

  #���W�b�g�Ɗm���̌v�Z
  logit <- matrix(X1 %*% b1 + X2 %*% b2, nrow=hh, ncol=select, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=select)
  
  #�����t���ΐ��ޓx���`
  LLi <- rowSums(y*log(Pr)) - lambda*sum(abs(b1)) - lambda*sum(abs(b2))
  LL <- sum(LLi)
  return(LL)
}

#�^�̑ΐ��ޓx
LLt <- loglike(b1=b0[1], b2=b0[-1], lambda=lambda, y=y, X1=as.matrix(XM_vec[, 1]), X2=XM_vec[, -1], hh=hh, select=select)


##�ΐ��ޓx���ő剻
##�����l�̐ݒ�
lambda <- 0.0025
theta <- rep(0, (select-1)*(k+1))

##�A���S���Y���̐ݒ�
max.iter <- 30
iter <- 1
tol <- 10
diff <- 100
L1 <- 0

##L1�������������W�b�g���f���𐄒�
while(diff >= tol & iter <= max.iter){
  for(i in 1:length(b0)){
    #�p�����[�^�̐ݒ�
    b1 <- theta[i]
    b2 <- theta[-i]

    #�����t���ΐ��ޓx���ő剻
    res <- optim(b1, loglike, gr=NULL, b2, lambda, y, as.matrix(XM_vec[, i]), XM_vec[, -i], hh, select, method="Nelder-Mead",
                 hessian=FALSE, control=list(fnscale=-1))
  
    #�p�����[�^�̍X�V
    theta[i] <- res$par
    print(c(i, LLt, res$value))
  }
  
  #�A���S���Y���̃p�����[�^�̍X�V
  iter <- iter+1
  LL <- res$value
  diff <- abs(LL-L1)
  L1 <- LL
}

round(cbind(theta, b0), 3)
colMeans(y)
loglike()

sum(abs(theta) >= 0.1)