#####�ł��؂�f�[�^�̃��f�����O#####
library(MASS)
library(reshape2)
library(plyr)

####�Б��ł��؂胂�f��####
####�f�[�^�̔���####
n <- 10000   #�T���v����
p <- 15   #�����ϐ���
b <- runif(p, -1.5, 4.5)   #��A�W��
b0 <- 8.4   #�ؕ�
sigma <- 8   #�W���΍�
X <- matrix(runif(n*p, -1.0, 5.0), n, p)   #�����ϐ�

betaT1 <- c(b, b0, sigma)

#�^�̃f�[�^�𔭐�
D <- trunc(X %*% b + b0 + rnorm(n, 0, sigma))   #�^�̎��v�֐�
S <- trunc(X %*% b + b0 + runif(n, 0, 2.5))   #�^�̋����֐�

#�w���f�[�^�𔭐�(���v�������������Ă���ꍇ�������w���f�[�^�Ƃ���)
B <- ifelse(D > S, S, D)
(BDS <- data.frame(B, D, S))

#�ł��؂�f�[�^�̎w���ϐ����쐬
z1 <- subset(1:n, BDS$D < BDS$S)   #���v����������Ă���f�[�^
z2 <- subset(1:n, BDS$D > BDS$S)   #���v�������������Ă���f�[�^
length(z1); length(z2)

####�ł��؂�f�[�^���f���𐄒�#####
##�ΐ��ޓx�֐��̒�`
fr <- function(theta, D, B, X, z1, z2, p){
  beta <- theta[1:p]
  beta0 <- theta[p+1]
  sigma <- exp(theta[p+2])   #�񕉐���
  Xb <- beta0 + as.matrix(X) %*% as.vector(beta)   #���σx�N�g��
  
  #��ł��؂�f�[�^�̖ޓx
  L1 <- sum(-log(sigma^2) - ((D - Xb)[z1])^2 / sigma^2)
  
  #�ł��؂�f�[�^�̖ޓx
  Lt <- 1-pnorm((B - Xb)[z2] / sigma)
  if(sum(Lt==0)!=0){i <- subset(1:length(Lt), Lt==0); Lt[i] <- 10^-100}
  L2 <- sum(log(Lt))   
  
  #�ΐ��ޓx�����v
  LL <- sum(L1 + L2)
  return(LL)
}

##�����l�̐ݒ�
fitf <- lm(B ~ X)
betaf <- fitf$coef[2:16]
betaf0 <- fitf$coef[1]
betaff <- as.numeric(c(betaf, betaf0, 1))

##�ΐ��ޓx�̍ő剻
fit <- optim(betaff, fn=fr, gr=NULL, D, B, X, z1, z2, p, 
             method="BFGS", hessian=T, control=list(fnscale=-1))

##���ʂƓ��v��
round(b <- c(fit$par[1:16], exp(fit$par[17])), 3)   #���肳�ꂽ�p�����[�^
round(betaT1, 3)   #�^�̌W��


c(b[1:16], log(b[17]))/sqrt(-diag(solve(fit$hessian)))   #t�l
(AIC <- -2*fit$value + 2*length(fit$par))   #AIC
(BIC <- -2*fit$value + log(nrow(X))*length(fit$par))   #BIC



####�����ł��؂胂�f��####
####�f�[�^�̔���####
n <- 300000   #�T���v����
p <- 20   #�����ϐ���

##�ő�l��9�ȉ��A�ŏ��l��-4�ȏ�ɂȂ�悤�ɕϐ��Ɖ�A�W�����쐬
T <- 10000
for(t in 1:T){
  b <- c(rnorm(12, 0.24, 0.18), rnorm(8, -0.21, 0.13))   #��A�W��
  b0 <- 0.6   #�ؕ�
  sigma <- 0.5   #�W���΍�
  X1 <- matrix(trunc(rnorm(n*p, 3, 1)), n, p)   #�����ϐ�
  X2 <- ifelse(X1 > 5, 5, X1)   #5�ȏ�̐��l��5�ɂ���
  X <- ifelse(X2 < 1, 1, X2)   #1�ȉ��̐��l��1�ɂ���
  
  #�^�̃f�[�^�̔���
  S <- b0 + X %*% b + rnorm(n, 0, sigma)   #�^�̃X�R�A
  print(min(S))
  print(max(S))
  if(max(S) < 15 && max(S) > 8 && min(S) < -2 && min(S) > -9) break
}
score_true <- round(S, 0)   #�^�̃X�R�A 
cbind(S, score_true)   
betaT2 <- c(b, b0, sigma)   #�^�̉�A�W��


##1�ȉ������5�ȏ�̃X�R�A�Ƀt���O�����āA�ϑ��f�[�^�̃X�R�A�ɂ���
#�X�R�A��5�ȏ�
upper.z <- ifelse(score_true >= 5, 1, 0)
table(upper.z)   
yupper <- ifelse(score_true >= 5, 5, score_true)

#�X�R�A��1�ȉ�
lower.z <- ifelse(score_true <= 1, 1, 0)
table(lower.z)
yobs <- ifelse(score_true <= 1, 1, yupper)

table(yobs)   #�X�R�A�̕��z������

####�ł��؂�f�[�^���f���𐄒�#####
##�ΐ��ޓx�֐��̒�`
fr <- function(theta, y, X, upper, lower, p){
  beta <- theta[1:p]
  beta0 <- theta[p+1]
  sigma <- exp(theta[p+2])   #�񕉐���
  Xb <- beta0 + as.matrix(X) %*% as.vector(beta)   #���σx�N�g��
  z <- abs(upper + lower - 1)   #��ł��؂�f�[�^�̎w���ϐ�
  
  #��ł��؂�f�[�^�̖ޓx
  L1 <- sum(-log(sigma^2) - (y[z==1] - Xb[z==1])^2 / sigma^2)
  
  #�㑤�ł��؂�f�[�^�̖ޓx
  Lupper <- 1-pnorm((y[upper==1] - Xb[upper==1]) / sigma);
  if(sum(Lupper==0)!=0){i <- subset(1:length(Lupper), Lupper==0); Lupper[i] <- 10^-100}
  L2 <- sum(log(Lupper))   
  
  #�����ł��؂�f�[�^�̖ޓx
  Llower <- pnorm((y[lower==1] - Xb[lower==1]) / sigma);
  if(sum(Llower==0)!=0){i <- subset(1:length(Llower), Llower==0); Llower[i] <- 10^-100}
  L3 <- sum(log(Llower))   
  
  #�ΐ��ޓx�����v
  LL <- sum(L1 + L2 + L3)
  return(LL)
}

##�����l�̐ݒ�
fitf <- lm(yobs ~ X)
betaf <- fitf$coef[2:(p+1)]
betaf0 <- fitf$coef[1]
betaff <- as.numeric(c(betaf, betaf0, 1))

##�ΐ��ޓx�̍ő剻
fit <- optim(betaff, fn=fr, gr=NULL, yobs, X, upper.z, lower.z, p, 
             method="BFGS", hessian=T, control=list(fnscale=-1))


####���ʂƓ��v��####
round(b <- c(fit$par[1:(p+1)], exp(fit$par[length(fit$par)])), 3)   #���肳�ꂽ�p�����[�^
round(betaT2, 3)   #�^�̌W��
round(betaff, 3)   #�ŏ����@�ł̌W��

#�^�̕��ύ\��
round(ym <- X %*% b[1:p] + b[(p+1)], 3)
cbind(score=round(ym, 0), score_true)   #�^�̃X�R�A�Ƃ̔�r
table(round(ym, 0))   #���茋�ʂ̃X�R�A�̕��z
table(score_true)   #�^�̃X�R�A�̕��z

c(b[1:(p+1)], log(b[length(b)]))/sqrt(-diag(solve(fit$hessian)))   #t�l
(AIC <- -2*fit$value + 2*length(fit$par))   #AIC
(BIC <- -2*fit$value + log(nrow(X))*length(fit$par))   #BIC