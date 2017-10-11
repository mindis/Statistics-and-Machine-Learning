#####�����t���l�X�e�b�h���W�b�g���f��#####
library(MASS)
library(mlogit)
library(flexmix)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(8645)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 10000   #����Ґ�
pt <- 50   #�ϑ�����
hhpt <- hh*pt   #���T���v����
choise <- 6   #�u�����h��

##ID�̐ݒ�
id <- rep(1:hh, rep(pt, hh))
time <- rep(1:pt, hh)
ID <- data.frame(no=1:hhpt, id=id, t=time)

####�����ϐ��̔���####
##�u�����h�w�����f���̐����ϐ�
#�O���̏���x�o�̔���
Consume <- rep(0, hhpt)
for(i in 1:hh){
  x <- runif(1, 0.2, 0.7)
  Consume[ID$id==i] <- rbinom(pt, 1, x)
}

#�J�e�S�����C�����e�B�̔���
ROYL <- rep(0, hhpt)
for(i in 1:hh){
  x <- rnorm(1)
  roy <- x + rnorm(pt, 0, runif(1, 0.1, 0.4))
  ROYL[ID$id==i] <- roy 
}

#�Ƒ��l���̔���
FAMILY <- rep(0, hhpt)
for(i in 1:hh){
  x <- round(rgamma(1, 5.0, 2), 0)
  x[x==0] <- 1
  FAMILY[ID$id==i] <- x
}

#�f�[�^�̌���
X1 <- data.frame(r=1, Cons=Consume, Last=0, Royl=ROYL, logsum=0, Family=FAMILY)
XM1 <- as.matrix(X1)


##�u�����h�I�����f���̐����ϐ�
#�ʏ퉿�i�̔���
PRICE <- matrix(runif(hhpt*choise, 0.6, 1), nrow=hhpt, ncol=choise) - 1   

#�f�B�X�J�E���g���̔���
DISC <- matrix(runif(hhpt*choise, 0, 0.4), nrow=hhpt, ncol=choise) 

#���ʒ�̔��� 
DISP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hhpt, 1, r)
}

#���ʃL�����y�[���̔���
CAMP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hhpt, 1, r)
}

#�u�����h���C�����e�B�̔���
BRAND <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  x_mean <- runif(choise, -0.7, 0.7)
  x_cov <- runif(choise, 0.1^2, 0.3^2)
  BRAND[ID$id==i, ] <- mvrnorm(pt, x_mean, diag(x_cov))
}


##�����ϐ����x�N�g����
#�ؕЂ̃x�N�g����
BP <- matrix(as.numeric(diag(choise)), nrow=hhpt*choise, ncol=choise, byrow=T)[, -choise]

#�J�e�S�����C�����e�B�̃x�N�g����
Royl_vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
for(i in 1:hhpt){
  r <- ((i-1)*choise+1):((i-1)*choise+choise)
  Royl_vec[r, ] <- diag(ROYL[i], choise)
}
Royl_vec <- Royl_vec[, -choise]

#���̑��̕ϐ��̃x�N�g����
Price_vec <- as.numeric(t(PRICE))
Disc_vec <- as.numeric(t(DISC))
Disp_vec <- as.numeric(t(DISP))
Camp_vec <- as.numeric(t(CAMP))
Brand_vec <- as.numeric(t(BRAND))

#ID���x�N�g����
id_vec <- rep(1:hh, rep(pt*choise, hh))
time_vec <- rep(rep(1:pt, rep(choise, pt)), hh)
c_vec <- rep(1:choise, hhpt)
ID_vec <- data.frame(no=1:(hhpt*choise), id=id_vec, time=time_vec, brand=c_vec)

#�f�[�^������
X2 <- data.frame(BP=BP, Price=Price_vec, Disc=Disc_vec, Disp=Disp_vec, Camp=Camp_vec, Brand=Brand_vec, Royl=Royl_vec)
XM2 <- as.matrix(X2)


####�����ϐ��̔���####
##�p�����[�^�̐ݒ�
#�u�����h�w�����f���̃p�����[�^
alpha00 <- -1.4   #�ؕ�
alpha01 <- runif(1, 0.7, 1.1)   #�O���̏���L���̃p�����[�^
alpha02 <- runif(1, 0.12, 0.18)   #�O��w������̌o�ߎ��Ԃ̃p�����[�^
alpha03 <- runif(1, 0.6, 0.9)   #�J�e�S�����C�����e�B�̃p�����[�^
alpha04 <- runif(1, 0.3, 0.7)   #���O�T���ϐ��̃p�����[�^
alpha05 <- runif(1, 0.05, 0.15)   #�Ƒ��l���̃p�����[�^
alpha0 <- c(alpha00, alpha01, alpha02, alpha03, alpha04, alpha05)

#�u�����h�I�����f���̃p�����[�^
beta00 <- c(1.2, 0.9, 2.0, 1.6, 0.6)
beta01 <- runif(1, -4.0, -3.6)   #���i�̃p�����[�^
beta02 <- runif(1, -4.5, -3.9)   #�������̃p�����[�^
beta03 <- runif(1, 2.4, 2.9)   #���ʒ�̃p�����[�^
beta04 <- runif(1, 2.1, 2.7)   #����L�����y�[���̃p�����[�^
beta05 <- runif(1, 0.6, 0.9)   #�u�����h���C�����e�B�̃p�����[�^
beta06 <- c(0.5, -0.2, 0.6, 0.9, -0.5)   #�J�e�S�����C�����e�B�̃p�����[�^
beta0 <- c(beta00, beta01, beta02, beta03, beta04, beta05, beta06)

rho1 <- 0.3   #�N���X�^�[1(�u�����h1�A2�A3)�̑��փp�����[�^
rho2 <- 0.5   #�N���X�^�[2(�u�����h4�A5)�̑��փp�����[�^
rho0 <- c(rho1, rho2, 1)

##���O�T���ϐ��̍쐬
nest <- cbind(c(1, 1, 1, 0, 0, 0), c(0, 0, 0, 1, 1, 0), c(rep(0, choise-1), 1))   #�l�X�g�\�����` 
nest1 <- matrix(c(1, 1, 1, 0, 0, 0), nrow=hhpt, ncol=choise, byrow=T)

#���p���`
logit <- matrix(XM2 %*% beta0, nrow=hhpt, ncol=choise, byrow=T)

#�u�����h�I�����f���̃��O�T���ϐ��̒�`
nest_list <- list()
logsum02 <- matrix(0, nrow=hhpt, ncol=length(rho0))
logsum01 <- rep(0, hhpt)
Pr02 <- matrix(0, nrow=hhpt, ncol=choise)

for(i in 1:ncol(nest)){
  nest_list[[i]] <- matrix(nest[, i], nrow=hhpt, ncol=choise, byrow=T)
  U <- exp(logit * nest_list[[i]] / rho0[i]) * nest_list[[i]]
  Pr02[, nest[, i]==1] <- U[, nest[, i]==1] / rowSums(U)   #�ŉ��w�̏����t���m��
  logsum02[, i] <- log(rowSums(U))   #���O�T���ϐ�
}

#�u�����h�I���m�����v�Z
Pr2 <- matrix(0, nrow=hhpt, ncol=choise)
V <- exp(logsum02 * matrix(rho0, nrow=hhpt, ncol=length(rho0), byrow=T))
CL <- V / rowSums(V)   #�l�X�g�̑I���m��

#�ŏI�I�ȃu�����h�I�����v�Z
for(i in 1:ncol(nest)){
  Pr2[, nest[, i]==1] <- matrix(CL[, i], nrow=hhpt, ncol=sum(nest[, i])) * Pr02[, nest[, i]==1]
}


#�u�����h�w�����f���̃��O�T���ϐ��̒�`
logsum01 <- log(rowSums(exp(logsum02)))
CV <- logsum01 / mean(logsum01)   #�l���傫���̂ŕ��ς�����
X1$logsum <- CV


##�O��w������̌o�ߎ��Ԃ̏����l��ݒ�
week <- rpois(hh, 3.2)
X1$Last[ID$t==1] <- ifelse(week==0, 1, week)
XM1 <- as.matrix(X1)


##�����I�ɍw���L���ƃu�����h�I���𔭐�������
y1 <- rep(0, hhpt)
y2 <- rep(0, hhpt)
Y2 <- matrix(0, nrow=hhpt, ncol=choise)

for(j in 1:pt){
  
  ##�u�����h�w���L���𔭐�
  #���p�Ɗm�����v�Z
  logit1 <- XM1[ID$t==j, ] %*% alpha0
  Pr1 <- exp(logit1) / (1+exp(logit1))
  
  #�x���k�[�C���z����w���L���𔭐�
  y1[ID$t==j] <- rbinom(hh, 1, Pr1)
  
  #�w���Ԋu���X�V
  X1$Last[ID$t==j+1] <- ifelse(y1[ID$t==j]==1, 1, XM1[ID$t==j, "Last"]+1)
  XM1 <- as.matrix(X1)
  
  ##�u�����h�w�����������ꍇ�����u�����h�I���𔭐�
  Y2[ID$t==j, ] <- t(apply(Pr2[ID$t==j, ], 1, function(x) rmultinom(1, 1, x))) * matrix(y1[ID$t==j], nrow=hh, ncol=choise)
  y2[ID$t==j] %*% Y2[ID$t==j, ] %*% 1:choise
}

##�����������ϐ��̏W�v
colSums(Y2)
colMeans(Y2[y1==1, ])
barplot(colSums(Y2))


####�Ŗޖ@�ŏ����t���l�X�e�b�h���W�b�g���f���𐄒�####
##�����t���l�X�e�b�h���W�b�g���f���̑ΐ��ޓx�֐���ݒ�
loglike <- function(x, y1, y2, X1, X2, nest, index1, index2, index3, hhpt, choise){
  #�p�����[�^�̐ݒ�
  alpha <- x[index1]
  beta <- x[index2]
  rho <- c(x[index3], 1)
  
  ##�u�����h�I�����f���̊m���̒�`
  logit2 <- matrix(X2 %*% beta, nrow=hhpt, ncol=choise, byrow=T)
  
  #�u�����h�I�����f���̃��O�T���ϐ��̒�`
  nest_list <- list()
  logsum02 <- matrix(0, nrow=hhpt, ncol=length(rho))
  logsum01 <- rep(0, hhpt)
  Pr02 <- matrix(0, nrow=hhpt, ncol=choise)
  
  for(i in 1:ncol(nest)){
    nest_list[[i]] <- matrix(nest[, i], nrow=hhpt, ncol=choise, byrow=T)
    U <- exp(logit2 * nest_list[[i]] / rho[i]) * nest_list[[i]]
    Pr02[, nest[, i]==1] <- U[, nest[, i]==1] / rowSums(U)   #�ŉ��w�̏����t���m��
    logsum02[, i] <- log(rowSums(U))   #���O�T���ϐ�
  }
  
  #�u�����h�I���m�����v�Z
  Pr2 <- matrix(0, nrow=hhpt, ncol=choise)
  V <- exp(logsum02 * matrix(rho, nrow=hhpt, ncol=length(rho), byrow=T))
  CL <- V / rowSums(V)   #�l�X�g�̑I���m��
  
  #�ŏI�I�ȃu�����h�I�����v�Z
  for(i in 1:ncol(nest)){
    Pr2[, nest[, i]==1] <- matrix(CL[, i], nrow=hhpt, ncol=sum(nest[, i])) * Pr02[, nest[, i]==1]
  }
  
  ##�w�����f���̊m���̒�`
  #�u�����h�w�����f���̃��O�T���ϐ��̒�`
  logsum01 <- log(rowSums(exp(logsum02)))
  CV <- logsum01 / mean(logsum01)   #�l���傫���̂ŕ��ς�����
  X1[, "logsum"] <- CV
  
  #���W�b�g�Ɗm���̒�`
  logit1 <- X1 %*% alpha
  Pr1 <- exp(logit1)/(1+exp(logit1))
  
  ##�ΐ��ޓx���`
  LL1 <- sum(y2 * log(Pr2))
  LL2 <- sum(y1*log(Pr1) + (1-y1)*log(1-Pr1))
  LL <- LL1 + LL2
  return(LL)
}


##�p�����[�^�̃C���f�b�N�X�̐ݒ�
index1 <- 1:ncol(XM1)
index2 <- (length(index1)+1):(length(index1)+ncol(XM2))
index3 <- (index2[length(index2)]+1):(index2[length(index2)]+2)

##���j���[�g���@�ŏ����t���l�X�e�b�h���W�b�g���f���̃p�����[�^�𐄒�
for(i in 1:1000){
  x <- c(rep(0, index2[length(index2)]), runif(2, 0.4, 0.7))
  res <- try(optim(x, loglike, gr=NULL, y1=y1, y2=Y2, X1=XM1, X2=XM2, nest=nest, index1=index1, index2=index2, index3=index3,
                   hhpt=hhpt, choise=choise, method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE)), silent=FALSE)
  if(class(res)=="try-error") {next} else {break}   #�G���[����
}

##���茋�ʂƓK���x
alpha <- res$par[index1]
beta <- res$par[index2]
rho <- res$par[index3]
LL <- res$value

#���肳�ꂽ�p�����[�^�Ɛ^�̃p�����[�^���r
round(rbind(alpha, alpha0), 3)
round(rbind(beta, beta0), 3)
round(rbind(rho, rho0=rho0[1:2]), 3)

#�K���x�̌v�Z
round(LL, 3)   #�ΐ��ޓx
round(tval <- res$par/sqrt(-diag(solve(res$hessian))), 3)   #t�l
round(AIC <- -2*res$value + 2*length(res$par), 3)   #AIC
round(AIC <- -2*res$value + log(hhpt)*length(res$par), 3)   #BIC


round(XM2, 3
      )