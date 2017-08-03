#####Paired combinatorial logit model#####
library(MASS)
library(mlogit)
library(nnet)
library(flexmix)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####���ϗʐ��K�����𔭐�������֐�####
##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}


####�f�[�^�̔���####
####�f�[�^�̐ݒ�####
hh <- 10000   #�T���v����
k <- 5   #�I��

####�����ϐ��̔���####
#�����t���̐����ϐ��̔���
X1.cont <- matrix(rnorm(2*(hh*k), 0, 1), nrow=hh, ncol=k*2)

X1.bin <- matrix(0, nrow=hh, ncol=k*2)
for(i in 1:(k*2)){
  X1.bin[, i]  <- rbinom(hh, 1, runif(1, 0.35, 0.6))
}

#�����^�̐����ϐ��̔���
X2.cont <- matrix(rnorm(hh, 0, 1), nrow=hh, ncol=1)
X2.bin <- matrix(rbinom(hh, 1, runif(1, 0.35, 0.7)), nrow=hh, ncol=1)


##�����ϐ����x�N�g���`���̃f�[�^�t�H�[�}�b�g�ɕύX
#ID��ݒ�
id <- rep(1:hh, rep(k, hh))
choise <- rep(1:k, hh)
ID <- data.frame(no=1:length(id), id=id, choise=choise)

#�ؕЂ̐ݒ�
p <- c(1, rep(0, k))
Pop <- matrix(p, nrow=hh*length(p), ncol=k, byrow=T)
Pop <- subset(Pop, rowSums(Pop) > 0)[, -k]


#�����^�����ϐ����x�N�g���`���ɐݒ�
X2v.cont <- matrix(0, hh*k, ncol=k)
X2v.bin <- matrix(0, hh*k, ncol=k)

for(i in 1:hh){
  index.v <- ((i-1)*k+1):((i-1)*k+k)
  v.cont <- diag(X2.cont[i, ], k)
  v.bin <- diag(X2.bin[i, ], k)
  X2v.cont[index.v, ] <- v.cont 
  X2v.bin[index.v, ] <- v.bin
}
X2v.cont <- X2v.cont[, -k]
X2v.bin <- X2v.bin[, -k]

#�����t�������ϐ����x�N�g���`���ɐݒ�
X1v.cont <- matrix(0, nrow=hh*k, ncol=2)
X1v.bin <- matrix(0, nrow=hh*k, ncol=2)

for(i in 1:2){
  index.r <- ((i-1)*k+1):((i-1)*k+k)
  X1v.cont[, i] <- as.numeric(t(X1.cont[, index.r]))
  X1v.bin[, i] <- as.numeric(t(X1.bin[, index.r]))
}

##�f�[�^������
X <- data.frame(pop=Pop, c1=X1v.cont, b1=X1v.bin, c2=X2v.cont, b2=X2v.bin)
round(XM <- as.matrix(X), 3)


####PCL���f���Ɋ�Â������ϐ��𔭐�####
##�p�����[�^�̐ݒ�
#��A�p�����[�^�̐ݒ�
b0 <- runif(k-1, -1.2, 2.0)
b1 <- runif(2, 0, 1.3)
b2 <- runif(2, -1.2, 1.4)
b3 <- runif(k-1, 0, 1.5)
b4 <- runif(k-1, -1.3, 1.5)
b <- c(b0, b1, b2, b3, b4)
beta.t <- b

#�ގ��x�p�����[�^�̐ݒ�
tau <- matrix(0, nrow=k, ncol=k-1) 
Cov <- corrM(col=k, lower=0, upper=0.8)   #�ގ��x�p�����[�^�𔭐�

for(i in 1:k){
  r <- 1:k
  tau[i, ] <- Cov[i, r[-i]]
}

#PCL���f���Ɋ�Â��m���̌v�Z
logit <- matrix(XM %*% b, nrow=hh, ncol=k, byrow=T)   #���p�֐��̐ݒ�
P1 <- matrix(0, nrow=hh, ncol=k)   #�m���̕��q�����̃p�����[�^�̊i�[�p�z��
P2 <- matrix(0, nrow=hh, ncol=k-1)   #�m���̕��ꕔ���̃p�����[�^�̊i�[�p�z��

#�m���v�Z�̕K�v�������v�Z
for(i in 1:k){
  #���q�����̌v�Z
  tau.m1 <- matrix(1-tau[i, ], nrow=hh, ncol=k-1, byrow=T)
  tau.m2 <- matrix(-tau[i, ], nrow=hh, ncol=k-1, byrow=T)
  logit1 <- matrix(logit[, i], nrow=hh, ncol=k-1)
  logit2 <- logit[, (1:k)[-i]]
  P1[, i] <- rowSums(tau.m1 * (exp(logit1/tau.m1) + exp(logit2/tau.m1))^tau.m2 * exp(logit1/tau.m1))
  
  if(i==k) {break}
  r <- (1:k)[-(1:i)]
  tau.m3 <- matrix(cbind(0, 1-tau)[i, r], nrow=hh, ncol=length(r), byrow=T)
  logit3 <- matrix(logit[, i], nrow=hh, ncol=length(r))
  logit4 <- logit[, ((i+1):k)]
  P2[, i] <- rowSums(tau.m3 * (exp(logit3/tau.m3) + exp(logit4/tau.m3))^tau.m3)
}

#�m���̌v�Z
Pr <- P1 / matrix(rowSums(P2), nrow=hh, ncol=k)

Pr1 <- matrix(exp(XM %*% b), hh, k, byrow=T) / rowSums(matrix(exp(XM %*% b), hh, k, byrow=T))
round(cbind(Pr, Pr1), 3)
round(Pr, 3)

##�����ϐ��𔭐�
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
colSums(Y)
round(data.frame(Y=Y, P=Pr), 3)

####Paired combinatorial logit model���Ŗސ���####
loglike <- function(x, Y, X, k, P1, P2){
  
  #�p�����[�^�̐ݒ�
  b <- x[1:ncol(X)]
  tau.v <- abs(x[(ncol(X)+1):length(x)])
  
  #�ގ��x�p�����[�^�̐ݒ�
  tau1 <- matrix(0, nrow=k, ncol=k-1)
  K <- diag(k)
  K[lower.tri(K)] <- tau.v
  K_t <- t(K) + K - diag(k)
  
  for(i in 1:k){
    r <- 1:k
    tau[i, ] <- K_t[i, r[-i]]
  }
  
  #PCL���f���Ɋ�Â��m���̌v�Z
  logit <- matrix(X %*% b, nrow=hh, ncol=k, byrow=T)   #���p�֐��̐ݒ�
  
  #�m���v�Z�̕K�v�������v�Z
  for(i in 1:k){
    #���q�����̌v�Z
    tau.m1 <- matrix(1-tau[i, ], nrow=hh, ncol=k-1, byrow=T)
    tau.m2 <- matrix(-tau[i, ], nrow=hh, ncol=k-1, byrow=T)
    logit1 <- matrix(logit[, i], nrow=hh, ncol=k-1)
    logit2 <- logit[, (1:k)[-i]]
    P1[, i] <- rowSums(tau.m1 * (exp(logit1/tau.m1) + exp(logit2/tau.m1))^tau.m2 * exp(logit1/tau.m1))
    
    if(i==k) {break}
    r <- (1:k)[-(1:i)]
    tau.m3 <- matrix(cbind(0, 1-tau)[i, r], nrow=hh, ncol=length(r), byrow=T)
    logit3 <- matrix(logit[, i], nrow=hh, ncol=length(r))
    logit4 <- logit[, ((i+1):k)]
    P2[, i] <- rowSums(tau.m3 * (exp(logit3/tau.m3) + exp(logit4/tau.m3))^tau.m3)
  }
  
  #�ΐ��ޓx���v�Z
  LLl <- rowSums(Y*log(P1)) - log(rowSums(P2))
  LL <- sum(LLl)
  return(LL)
}

##�������W�b�g���f���̑ΐ��ޓx�֐�
LL_logit <- function(x, X, Y, hh, k){
  #�p�����[�^�̐ݒ�
  theta <- x
  
  #���p�֐��̐ݒ�
  U <- matrix(X %*% theta, nrow=hh, ncol=k, byrow=T)
  
  #�ΐ��ޓx�̌v�Z
  d <- rowSums(exp(U))
  LLl <- rowSums(Y * U) - log(d)
  LL <- sum(LLl)
  return(LL)
}


##����t�����j���[�g���@��PCL���f�����Ŗސ���
#�p�����[�^�̏����l�����l�𑽍����W�b�g���f���Ō���
#���j���[�g���@�ōŖސ���
theta <- runif(ncol(XM), -1, 1)
res.z <- optim(theta, LL_logit, gr=NULL,  X=XM, Y=Y, hh=hh, k=k, method="BFGS", hessian=FALSE,
               control=list(fnscale=-1))   
theta <- res.z$par   #�p�����[�^�̏����l

P1 <- matrix(0, nrow=hh, ncol=k)   #�m���̕��q�����̃p�����[�^�̊i�[�p�z��
P2 <- matrix(0, nrow=hh, ncol=k-1)   #�m���̕��ꕔ���̃p�����[�^�̊i�[�p�z��

for(i in 1:1000){
  print(i)
  
  #�ގ��x�p�����[�^�̏����l
  tau.first1  <- corrM(col=k, lower=0, upper=0.7)   #�ގ��x�p�����[�^�𔭐�
  tau.first <- tau.first1[lower.tri(tau.first1)]

  #�����l���x�N�g���Ɍ���
  x <- c(theta, tau.first)
  
  #�p�����[�^�̏���l�Ɖ����l��ݒ�
  l <- c(theta+2, rep(0, sum(1:(k-1))))
  u <- c(theta-2, rep(1, sum(1:(k-1))))

  #���j���[�g���@�ōŖސ���
  res <- try(optim(x, loglike, gr=NULL, X=XM, Y=Y, k=k, P1=P1, P2=P2, method="SANN", hessian=TRUE, control=list(fnscale=-1)), 
             silent=FALSE)
  if(class(res)=="try-error") {next} else {break}   #�G���[����
}

####�p�����[�^���茋�ʂƗv��####
##���肳�ꂽ�p�����[�^�Ɛ^�̃p�����[�^�̔�r
round(beta <- res$par, 3)
round(c(beta.t, Cov[lower.tri(Cov)]), 3)

#�ގ��x�p�����[�^�𑊊փp�����[�^�ɕϊ�
rho.m <- diag(k)
rho <- beta[(ncol(XM)+1):length(beta)]
rho.m[lower.tri(rho.m)] <- rho
round(rho.M <- abs(t(rho.m) + rho.m - diag(k)), 3)
round(Cov, 3)

#�ő剻���ꂽ�ΐ��ޓx��AIC
round(res$value, 3)

round(tval <- beta/sqrt(-diag(solve(res$hessian))), 3)   #t�l
round(AIC <- -2*res$value + 2*length(beta), 3)   #AIC
round(BIC <- -2*res$value + log(hh)*length(beta), 3)   #BIC
