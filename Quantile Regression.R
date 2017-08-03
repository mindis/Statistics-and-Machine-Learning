#####���ʓ_��A����#####
#####�K�w�x�C�Y�������W�b�g���f��#####
library(quantreg)
library(MASS)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#�f�[�^�̐ݒ�
n <- 2000
k1 <- 4
k2 <- 2

##�����ϐ��̔���
#�A���ϐ�
X_cont <- matrix(runif(n*k1, 0, 2), nrow=n, ncol=k1)

#��l�ϐ�
X_bin <- matrix(0, nrow=n, ncol=k2)
for(i in 1:k2){
  X_bin[, i] <- rbinom(n, 1, runif(1, 0.3, 0.7))
}

#�����ϐ�������
X <- cbind(X_cont, X_bin)
Xi <- cbind(1, X)

##�p�����[�^��ݒ�
a <- runif(1, 4.0, 6.0)
b <- runif(ncol(X), 2.0, 5.0)


##�����ϐ��𔭐�
Y <- c()
b_qual <- matrix(0, nrow=n, ncol=k1)

for(i in 1:n){
  b_qual[i, ] <- b[1:k1] + mvrnorm(1, rep(0, 4), diag(X[i, 1:k1]))
  Y <- c(Y, a + X[i, 1:k1] %*% b_qual[i, ] + X[i, (k1+1):ncol(X)] %*% b[(k1+1):ncol(X)] + rnorm(1, 0, 1))
}

#�U�z�}
plot(X[, 1], Y, xlab="X�̒l")

#�ŏ������ŗ\��
round(b_sq <- as.numeric(solve(t(Xi) %*% Xi) %*% t(Xi) %*% Y), 3)
round(b, 3)
Y_sq <- Xi %*% b_sq   #�\���l
round(data.frame(Y, Y_sq), 3)


####���ʓ_��A�𐄒�####
##�ŏ�������֐���ݒ�(��Ό덷)
qual_reg <- function(beta, Y, X, tau){
  er <- Y - X %*% beta
  rho <- (abs(er) + (2*tau-1)*er)/2
  return(sum(rho))
}

##���ʓ_���Ƃɕ��ʓ_��A�𓖂Ă͂߂�
qual <- c(0.1, 0.25, 0.5, 0.75, 0.9)
res <- list()
beta <- matrix(0, nrow=length(qual), ncol=ncol(Xi))
er <- matrix(0, nrow=length(qual), ncol=ncol(Xi)) 
res_func <- list()
beta_func <- matrix(0, nrow=length(qual), ncol=ncol(Xi))

for(i in 1:length(qual)){
  b1 <- runif(ncol(Xi), 0, 2)
  res[[i]] <- optim(b1, qual_reg, Y=Y, X=as.matrix(Xi), tau=qual[i], method="CG")
  beta[i, ] <- res[[i]]$par   #���肳�ꂽ�p�����[�^
  er[i, ] <- res[[i]]$value   #�ŏ������ꂽ��Ό덷
  
  #�֐��ŕ��ʓ_��A
  res_func[[i]]  <- rq(Y ~ X, tau=qual[i])
  beta_func[i, ] <- as.numeric(res_func[[i]]$coefficients)
}

#���ʂ��r
round(beta, 3)
round(beta_func, 3)
round(b_sq, 3)


