#####����q�^�g�[�r�b�g���f��#####
library(MASS)
library(reshape2)
library(plyr)
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

##���֍s�񂩂番�U�����U�s����쐬����֐����`
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #�ŗL�l�����ŋ����I�ɐ���l�s��ɏC������
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}


####�f�[�^�̔���####
##�f�[�^�̐ݒ�
N <- 4000   #�T���v����
k <- 12    #�����ϐ���
k.cont <- 5    #�A���ϐ���
k.bin <- 3   #��l�ϐ���
k.multi <- 4    #���l�ϐ���

##ID�̐ݒ�
ID <- 1:N

##�����ϐ��̔���
#�A���ϐ�
X.cont <- matrix(rnorm(N*k.cont), nrow=N, ncol=k.cont)

#��l�ϐ�
X.bin <- matrix(0, nrow=N, ncol=k.bin)
for(i in 1:k.bin){
  p.bin <- runif(1, 0.3, 0.6)
  X.bin[, i] <- rbinom(N, 1, p.bin)
}

#���l�ϐ�
p.multi <- runif(k.multi)
X.multi <- t(rmultinom(N, 1, p.multi))
X.multi <- X.multi[, -which.min(colSums(X.multi))]

X <- data.frame(i=1, cont=X.cont, bin=X.bin, multi=X.multi)
XM <- as.matrix(X)

##�p�����[�^�̐ݒ�
#�v���r�b�g���f���̃p�����[�^�̐ݒ�
for(t in 1:1000){
  print(t)
  b1 <- c(-0.5, runif(ncol(X)-1, -0.7, 1.0))
  b2 <- c(1.6, runif(ncol(X)-1, -3.3, 4.4))
  
  #�����U�s��̐ݒ�
  Cov <- matrix(c(1, 1.2, 1.2, 4.5), nrow=2, ncol=2)
  
  ##2�ϗʉ�A���f��������ݕϐ��𔭐�
  er <- mvrnorm(N, rep(0, 2), Cov)   #�덷�𔭐������� 
  Z <- XM %*% cbind(b1, b2) + er   #���ݕϐ�Z�𔭐�
  
  #���ݕϐ��������ϐ�Y�ɕϊ�
  Y1 <- ifelse(Z[, 1] >= 0, 1, 0) 
  Y2 <- ifelse(Z[, 2] >= 0, Z[,2 ], 0) * Y1
  Y3 <- Y1 * ifelse(Y2 > 0, 0, 1)
  Y4 <- ifelse(Y2 > 0, 1, 0)
  
  agg <- colSums(cbind(Y1=(1-Y1), Y3, Y4))/N   #�����p�^�[���ʂ̔䗦
  if(max(agg) < 0.6 & min(agg) > 0.2) {break} else {next}
}

round(cbind(Y1, Y2, Y3, Y4), 3)   #�����ϐ����m�F
colSums(cbind(Y1=(1-Y1), Y3, Y4))/N   #�����p�^�[���ʂ̔䗦

####�Ŗސ���œ���q�^�g�[�r�b�g���f���𐄒�####
##����q�^�g�[�r�b�g���f���̑ΐ��ޓx�֐�
loglike <- function(b, Y1, Y2, Y3, Y4, X){
  
  #�p�����[�^��ݒ�
  beta1 <- b[1:ncol(X)]
  beta2 <- b[(ncol(X)+1):(2*ncol(X))]
  sigma1 <- 1
  sigma_sq <- b[(2*ncol(X))+1]
  rho <- b[(2*ncol(X))+2]
  sigma2 <- sqrt(sigma_sq)
  
  #���`���f�����v�Z
  bx <- X %*% beta1
  cw <- X %*% beta2
  
  #���`���f����W����
  z1 <- (Y1-bx) / sigma1
  z2 <- (Y2-cw) / sigma2
  
  #�ΐ��ޓx���v�Z
  d10 <- log(1-pnorm(bx/sigma1))
  d20 <- log(1/sigma1 *dnorm(z1) * (1-pnorm((1/(1-rho^2)^(1/2)) * (cw/sigma2 + rho*z1))))
  d21 <- 1/2*log(2*pi) - log((1/sigma1)*(1/sigma2)*(1/(1-rho^2)^(1/2))) + 
              1/(1-rho^2) * (z1^2+z2^2 - 2*rho*z1*z2)/2
  
  #nan��0�ɒu��������
  d10[is.infinite(d10)] <- -0.01
  d20[is.infinite(d20)] <- -0.01
  d21[is.infinite(d21)] <- 0.01
  
  #�ΐ��ޓx�̑��a���v�Z
  LL <- sum((1-Y1)*d10 + Y3*d20 - Y4*d21)
  return(LL)
}


##���j���[�g���@�œ���q�^�g�[�r�b�g���f���𐄒�
#�����p�����[�^��ݒ�
res <- list()
for(j in 1:20){
  print(j)
  for(t in 1:1000){
    b0 <- c(runif(ncol(X), -1.5, 1.5), runif(ncol(X), -2.2, 3.0), runif(1, 2, 4), runif(1, 0.2, 0.6))
    
    #���j���[�g���@�őΐ��ޓx���ő剻
    res[[j]] <- try(optim(b0, loglike, gr=NULL, Y1=Y1, Y2=Y2, Y3=Y3, Y4=Y4, X=XM, method="BFGS", 
                          hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)
    if(class(res[[j]]) == "try-error") {next} else {break}   #�G���[����
  }
}


####���肳�ꂽ�p�����[�^�̗v��Ɠ��v�ʂ̌v�Z####
##�x�X�g�ȃp�����[�^��ΐ��ޓx�Ō���
LL_res <- c()
for(i in 1:length(res)){
  LL_res <- c(LL_res, res[[i]]$value)
}
opt <- which.max(LL_res)

##�x�X�g�ȑΐ��ޓx�̃p�����[�^
beta <- res[[opt]]$par  

#���肳�ꂽ��A�W���Ɛ^�̉�A�W���̔�r
round(rbind(beta[1:ncol(X)], beta[(ncol(X)+1):(2*ncol(X))]), 2)
round(rbind(b1, b2), 2)

#���肳�ꂽ���U�p�����[�^�Ɛ^�̕��U�p�����[�^�̔�r
round(beta[(length(beta)-1):length(beta)], 3)
round(c(Cov[2, 2], cov2cor(Cov)[1, 2]), 3)

#�ő剻���ꂽ�ΐ��ޓx��AIC
round(res[[opt]]$value, 3)

round(tval <- beta/sqrt(-diag(solve(res[[opt]]$hessian))), 3)   #t�l
round(AIC <- -2*res[[opt]]$value + 2*length(beta), 3)   #AIC
round(BIC <- -2*res[[opt]]$value + log(N)*length(beta), 3)   #BIC


