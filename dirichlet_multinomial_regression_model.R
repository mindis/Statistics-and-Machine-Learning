#####�f�B�N����������A���f��#####
library(MASS)
detach("package:bayesm", unload=TRUE)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

set.seed(7853)

####�f�[�^�̔���####
hh <- 5000   #���[�U�[��
select <- 6   #�I������
st <- 6   #��I����

##�����ϐ��̔���
#�u�����h���C�����e�B
ROYL <- matrix(runif(hh*select), nrow=hh, ncol=select)

#����L��
PRE <- matrix(0, nrow=hh, ncol=select)
for(j in 1:select){
  par <- runif(1, 0.2, 0.4)
  PRE[, j] <- rbinom(hh, 1, par)
}

#�����̐ݒ�
DIST <- matrix(runif(hh*select), nrow=hh, ncol=select)


#�Ƒ��\��
FAMI0 <- rpois(hh, 2.4)
FAMI <- ifelse(FAMI0==0, 1, FAMI0)

#�L���ڐG���x��
PROM0 <- round(rnorm(hh, 3, 1))
PROM <- ifelse(PROM0 < 1, 1, ifelse(PROM0 > 5, 5, PROM0))

##�����ϐ����x�N�g���ϊ�
#ID��ݒ�
u_id <- rep(1:hh, rep(select, hh))
i_id <- rep(1:select, hh)
ID_vec <- data.frame(no=1:length(u_id), id=u_id, item=i_id)

#�ؕЂ̐ݒ�
BP <- matrix(diag(1, select), nrow=hh*select, ncol=select, byrow=T)[, -st]

#�I�������Ƃɐ����ϐ����ς������ϐ����x�N�g���ϊ�
ROYL_vec <- as.numeric(t(ROYL))
DIST_vec <- as.numeric(t(DIST))
PRE_vec <- as.numeric(t(PRE))

#�I���d���ɐ����ϐ����ς��Ȃ������ϐ����x�N�g���ϊ�
FAMI_vec0 <- matrix(0, nrow=hh*select, ncol=select)
PROM_vec0 <- matrix(0, nrow=hh*select, ncol=select)

for(i in 1:hh){
  FAMI_vec0[ID_vec$id==i, ] <- diag(FAMI[i], select)
  PROM_vec0[ID_vec$id==i, ] <- diag(PROM[i], select)
}

FAMI_vec <- FAMI_vec0[, -st]
PROM_vec <- PROM_vec0[, -st]

#�f�[�^�̌���
X <- cbind(b=BP, roy=ROYL_vec, pre=PRE_vec, dist=DIST_vec, fam=FAMI_vec, prom=PROM_vec)


####�f�B�N�������z���牞���ϐ��𔭐�####
##�p�����[�^��ݒ�
beta00 <- runif(select-1, -0.7, 0.9)
beta01 <- c(runif(2, 0, 0.8), runif(1, -1.1, -0.6), runif(select-1, 0.075, 0.125), runif(select-1, 0.1, 0.2))
theta0 <- c(beta00, beta01)

##�f�B�N�������z���牞���ϐ��𔭐�
alpha0 <- matrix(exp(X %*% theta0), nrow=hh, ncol=select, byrow=T)   #�f�B�N�������z�̕��ύ\��
y <- t(apply(alpha0, 1, function(x) rdirichlet(1, x)))   #�f�B�N�������z����I���m���𔭐�
summary(y)   #�W�v


####�Ŗސ���Ńf�B�N����������A���Ŗސ���####
##�f�B�N����������A���f���̑ΐ��ޓx�֐����`
fr <- function(theta, y, X, hh, select){
  #�p�����[�^��ݒ�
  beta <- theta
  
  #�f�B�N�������z�̃p�����[�^�̕��ύ\����ݒ�
  alpha <- matrix(exp(X %*% beta), nrow=hh, ncol=select, byrow=T)
  
  #�f�B�N�������z�̑ΐ��ޓx��a
  LLi <- log(ddirichlet(y, alpha))
  LL <- sum(LLi)
  return(LL)
}

##���j���[�g���@�őΐ��ޓx���ő剻
#�p�����[�^�̏����l
beta0 <- runif(select-1, -0.5, 0.5)
beta1 <- c(runif(2, 0, 0.5), runif(1, -0.6, -0.1), runif(select-1, 0, 0.1), runif(select-1, 0, 0.1))
theta <- c(beta0, beta1)

#�ΐ��ޓx���ő剻����
res <- optim(theta, fr, gr=NULL, y, X, hh, select, method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))

##���茋�ʂ̊m�F�ƓK���x
theta <- res$par
LL <- res$value

#���肳�ꂽ�p�����[�^�Ɛ^�̃p�����[�^�̔�r
round(rbind(theta, theta0), 3)

#�K���x�̌v�Z
round(res$value, 3)   #�ő�ΐ��ޓx
round(tval <- res$par/sqrt(-diag(solve(res$hessian))), 3)   #t�l
round(AIC <- -2*res$value + 2*length(res$par), 3)   #AIC
round(BIC <- -2*res$value + log(hh)*length(res$par), 3) #BIC
