#####���̑������z���f��#####
library(MASS)
library(vcd)
library(gtools)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(534)
N <- 5000
k <- 20

#�K���}���z����p�����[�^�ɂ𔭐�������
theta.t <- matrix(0, nrow=k, ncol=2)
Y <- matrix(0, nrow=N, ncol=k)
lambda <- matrix(0, nrow=N, ncol=k)

for(i in 1:k){
  #�K���}���z����lambda�𔭐�
  theta.t[i, ] <- c(runif(1, 0.001, 3), 1.5)
  lambda[, i] <- rgamma(N, shape=theta.t[i, 1], scale=theta.t[i, 2])
  
  #lambda����|�A�\�������𔭐�������
  poison <- rpois(N, lambda[, i])  
  Y[, i] <- poison
  
}
#�����������f�[�^���m�F
round(data.frame(Y=Y, lambda=lambda), 1)
summary(Y)
round(apply(Y, 2, sd), 2)

#�f�[�^�̕��z
hist(Y[, 1], breaks=25, col="grey", xlab="�p�x", main="���̓񍀕��z")   #Y�̕��z
hist(lambda[, 1], breaks=25, col="grey", xlab="�p�x", main="lambda�̕��z")   #lambda�̕��z


##�Z�����ƕp�x���v�Z
#�Z�����̌v�Z
perm_rate <- c()
for(i in 1:k){
 perm_rate <- c(perm_rate, mean(ifelse(Y[, i]==0, 0, 1)))
}

#�p�x�̌v�Z
summary(rowSums(Y))
freq.all <- colSums(Y)/N   #�S�̂ł̕p�x
freq.cond <- colSums(Y)/colSums(apply(Y, 2, function(x) ifelse(x > 0, 1, 0)))

#�Z�����ƕp�x�̌������Ĕ�r
round(data.frame(perm_rate, freq.all, freq.cond), 3)   
plot(perm_rate, freq.cond, xlab="�Z����", ylab="�w���p�x")


####���̑������z���f�����Ŗސ���####
fr <- function(b0, Y, N, k){
  #�p�����[�^�̐ݒ�
  alpha <- exp(b0[1:k])
  scale <- exp(b0[k+1])
  
  #�ΐ��ޓx���v�Z
  LL <- sum(
            Y * log(scale/(1+scale)) +
            matrix(alpha*log(1/(1+scale)), nrow=N, ncol=k, byrow=T) +
            log(gamma(Y+matrix(alpha, nrow=N, ncol=k, byrow=T))) -
            lfactorial(Y) - log(matrix(gamma(alpha), nrow=N, ncol=k, byrow=T)) 
            )
  return(LL)
}

##�ΐ��ޓx���ő剻
b0 <- c(runif(k, 0.1, 2), 1.0)   #�����p�����[�^�̐ݒ�
res <- optim(b0, fr, gr=NULL, Y, N, k, method="BFGS", hessian=TRUE, control=list(fnscale=-1))   #���j���[�g���@�ōœK��

####���茋�ʂƉ���####
##���肳�ꂽ�p�����[�^
b <- res$par
round(exp(b), 2)   #���茋��
round(c(theta.t[, 1], theta.t[1, 2]), 2)   #�^�̃p�����[�^

alpha <- exp(b[1:k])
scale <- exp(b[k+1])

##�K���x��AIC
res$value   #�ő剻���ꂽ�ΐ��ޓx
(tval <- b/sqrt(-diag(solve(res$hessian))))   #t�l
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(length(poison))*length(b))   #BIC


##���茋�ʂ���Z�����ƍw���p�x���v�Z
#�Z�����̌v�Z
rate.nmd <- (1-(scale+1)^-alpha) / (1-(scale+1)^-sum(alpha))  
round(data.frame(rate.nmd, perm_rate), 3)

#�w���p�x�̌v�Z
freq.nmd <- (alpha*scale) / ((1-(scale+1)^-sum(alpha))*rate.nmd)   #�w���p�x
round(data.frame(freq.nmd, freq.cond), 2)


##�Z�����ƍw���p�x�̊֌W������
#�Z�����ƍw���p�x�̊֌W���֐��X���[�W���O�Ő���
lo <- loess(freq.nmd ~ rate.nmd)
x <-  seq(min(rate.nmd), max(rate.nmd), length=500)
pred.lo <- predict(lo, x)

#�֐��X���[�W���O(���)�ƐZ�����ƍw���p�x�̎����l���v���b�g
plot(perm_rate, freq.cond, xlab="�ϑ����ꂽ�Z����", ylab="�ϑ����ꂽ�w���p�x", main="�Z�����ƍw���p�x�̊֌W", 
     pch=3, cex=1.25)
lines(x, pred.lo, type="l", col=2)
