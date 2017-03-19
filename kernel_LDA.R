#####�J�[�l���������ʂ���уJ�[�l���听������#####
library(kernlab)
library(plyr)
library(MASS)
set.seed(555)
####�V�~�����[�V�����f�[�^�̔���####
##���̓��̓f�[�^���쐬
x1 <- matrix(rnorm(50000, 0, 1), 1000, 50)
v51 <- as.vector(rpois(1000, 8))
x3 <- t(rmultinom(1000, 1, c(0.2, 0.2, 0.3, 0.15, 0.1, 0.05)))
x4 <- t(rmultinom(1000, 1, c(0.1, 0.15, 0.15, 0.05, 0.2, 0.1, 0.2, 0.05)))
Xv <- cbind(x1, v51, x3, x4)

#�璷�ȕϐ����폜����
Xv <- Xv[, -65]
Xv <- Xv[, -57]

#�f�[�^�t���[���ɕϊ�
X <- as.data.frame(Xv)

#�������̕ϐ����J�e�S���J���ϐ��ɕϊ�
X[, 1:20] <- X[, 1:20] > runif(20, -0.7, 0.7)
Xnew <- as.data.frame(as.matrix(X))
round(head(Xnew, 5), 3)
summary(Xnew)

##�O�����s��ɕϊ�����
#�^�̃J�[�l���֐��̐ݒ�
n <- nrow(Xnew)
gm <- matrix(0, n, n)

##�O�����s����쐬����
#�������J�[�l��
gram1 <- (3 + as.matrix(Xnew) %*% t(as.matrix(Xnew))) + (3 + as.matrix(Xnew) %*% t(as.matrix(Xnew)))^2
round(gram1[1:15, 1:15], 3)
round(gram1[985:1000, 985:1000], 3)
round(eigen(gram1)$value, 3)   #������l���ǂ����m�F

#�K�E�X�J�[�l��
sigma <- 1/2
kf_gauss <- function(x1, sigma){
  x1 <- as.matrix(x1)
  g1 <- matrix(t(x1), nrow(x1)^2, ncol(x1), byrow=T)
  g2 <- matrix(rep(x1, c(rep(nrow(x1), nrow(x1)*ncol(x1)))), nrow(x1)^2, ncol(x1))
  g3 <- (g2 - g1)^2
  gram <- exp(-sigma*sqrt(matrix(apply(g3, 1, sum), nrow(x1), nrow(x1), byrow=T)))
  return(gram)
}
gram2 <- kf_gauss(x1=Xnew, sigma=sigma)
round(gram2[1:15, 1:15], 3)
round(gram2[985:1000, 985:1000], 3)
round(eigen(gram2)$value, 3)   #������l���ǂ����m�F

#�֐����g��
L <- kernelMatrix(rbfdot(sigma=1/2), Xnew)   #�K�E�X�J�[�l���ŕϊ�

#forloop�ŃO�����s����쐬����ƌv�Z���Ԃ�����
#for(i in 1:n){
#  for(j in 1:n){
#    r <- sum(Xnew[i, ] - Xnew[j, ])
#    gm[i, j] <- r
#  }
#  print(i)
#}

##�O�����s����g���ĕ��ނ��쐬����
#�W���f�[�^�̍쐬
a <- rnorm(1000, 0, 0.5)
#����̐��l�ȏゾ���W���Ƃ���
rand <- abs(rnorm(1000, 0, 3.5))
tcnt <- abs(a1) >= rand
table(tcnt)
alpha <- ifelse(abs(a1) >= rand, a1, rnorm(1000, 0, 0.1))   
round(alpha, 3)
z <- gram2 %*% as.vector(alpha)
plot(z)

library(mlbench)
train <- mlbench.2dnormals(500, cl=20, sd=0.5)   #�w�K�f�[�^
plot(train$x)