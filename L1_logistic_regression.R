#####L1���������W�X�e�B�b�N��A���f��#####
library(MASS)
library(ncvreg)
library(glmnet)
library(lars)
library(reshape2)
library(plyr)

####�f�[�^�̔���####
#set.seed(319)
p <- 100   #�����ϐ��̌�
pm <- 10   #�Ӗ��̂���ϐ��̌�
n <- 2000   #�T���v����
X <- matrix(rnorm(n*p), nrow=n, ncol=p)   #�f�[�^�̔���
b <- c(1.1, runif(pm, -1.0, 1.0), rep(0, length=(p-pm)))   #�^�̉�A�W��

Xb <- b[1] + X %*% b[2:length(b)]   #���`����
pr <- exp(Xb) / (1+exp(Xb))   #�^�̊m��
y <- rbinom(n, 1, pr)   #�ϑ��f�[�^�̐���


####L1���������W�X�e�B�b�N��A���f���𐄒�####
##�ΐ��ޓx�֐��̐ݒ�
#�ΐ��ޓx�̐ݒ�
fr <- function(b1, b2, lambda, x1, x2, y){
  #�p�����[�^�̐ݒ�
  alpha <- b1[1]
  beta <- b1[2]
  
  #�ޓx���`���č��v����
  Xb <- alpha + as.matrix(x1) %*% beta + as.matrix(x2) %*% b2 
  p <- exp(Xb) / (1 + exp(Xb))
  LLS <- y*log(p) + (1-y)*log(1-p) - lambda*abs(sum(b2)) - lambda*abs(beta)
  LL <- sum(LLS)
  return(LL)
}

##�����l�̐ݒ�
b0 <- c(0.5, runif(p, -1, 1))
lambda <- 0.005 

##�A���S���Y���̐ݒ�
max.iter <- 30
iter <- 1
tol <- 1
diff <- 100
L1 <- 0

##L1���������W�X�e�B�b�N��A���f���𐄒�
while(diff >= tol & iter <= max.iter){
  for(i in 1:p){
    b1 <- b0[c(1, (i+1))]   #�œK������p�����[�^
    b2 <- b0[c(-1,-(i+1))]   #�Œ肷��p�����[�^
    x1 <- X[, i]   #�œK������p�����[�^�ɑΉ���������ϐ�
    x2 <- X[, -i]   #�Œ肷��p�����[�^�ɑΉ���������ϐ�
      
    #�����l��Nelder-Mead�@�Ō��肷��
    res1 <- optim(b1, fr, gr=NULL, b2, lambda, x1, x2, y,
                  method="Nelder-Mead", hessian=FALSE, control=list(fnscale=-1))
    
    #�X�V�����p�����[�^�̊i�[
    b0[1] <- res1$par[1]
    b0[(i+1)] <- res1$par[2]
  }
  #�A���S���Y���̃p�����[�^�̍X�V
  iter <- iter+1
  LL <- res1$value
  diff <- abs(LL-L1)
  L1 <- LL
  
  #�ΐ��ޓx�̕\��
  print(diff)
}

res1$value   #�ΐ��ޓx
round(b0, 2)   #���肳�ꂽ�p�����[�^
round(b, 2)   #�^�̃p�����[�^


####�N���X�o���e�[�V�����ōœK��lambda�����߂�####
##�A���S���Y���ƃp�����[�^�̐ݒ�
#�p�����[�^�̐ݒ�
b0 <- c(0.5, runif(p, -1, 1))   #�����l
lambdaE <- seq(0.001, 0.03, length=12)   #lambda�̌��

##�A���S���Y���̐ݒ�
spl <- 5
len <- nrow(X)/spl   #�T���v����5����
max.iter <- 30   #�ő�J��Ԃ���
iter <- 1   #�J��Ԃ����̏����l
tol <- 1   #�ޓx�̍��̂������l
diff <- 100   #�ޓx�̍��̏����l
L1 <- 0   #�ޓx�̏����l
CR <- c()   #�������̊i�[�p

##5�����N���X�o���f�[�V�����ɂ��œK��lambda�̑I��
for(lam in 1:length(lambdaE)){
  lambda <- lambdaE[lam]
  cr <- c()   
  
  for(k in 1:spl){
    l <- ((k-1)*len+1):(k*len)
    index <- subset(l, l <= nrow(X))   #�f�U�C���s��̍s���ȏ�̃C���f�b�N�X�͍폜
    x.cv <- X[-index, ]
    y.cv <- y[-index]
    diff <- 100   #diff�̏�����
    iter <- 30   #iter�̏�����
    
    ##coordinate Descent�@�ɂ�鐄��
    while(diff >= tol & iter <= max.iter){
      for(i in 1:p){
        b1 <- b0[c(1, (i+1))]   #�œK������p�����[�^
        b2 <- b0[c(-1,-(i+1))]   #�Œ肷��p�����[�^
        x1 <- x.cv[, i]   #�œK������p�����[�^�ɑΉ���������ϐ�
        x2 <- x.cv[, -i]   #�Œ肷��p�����[�^�ɑΉ���������ϐ�
        
        #�����l��Nelder-Mead�@�Ō��肷��
        res1 <- optim(b1, fr, gr=NULL, b2, lambda, x1, x2, y.cv,
                      method="Nelder-Mead", hessian=FALSE, control=list(fnscale=-1))
        
        #�X�V�����p�����[�^�̊i�[
        b0[1] <- res1$par[1]
        b0[(i+1)] <- res1$par[2]
      }
      #�A���S���Y���̃p�����[�^�̍X�V
      iter <- iter+1
      LL <- res1$value
      diff <- abs(LL-L1)
      L1 <- LL
    }
    
    ##���肳�ꂽ�p�����[�^�̓K���x�̕]��
    #�p�����[�^�̊i�[
    beta0 <- b0[1]   #�ؕЂ̐���l
    beta1 <- b0[2:length(b0)]   #��A�p�����[�^�̐���l
    
    #�m���̌v�Z
    logi <- beta0 + X[index, ] %*% beta1   #�e�X�g�f�[�^�̃��W�b�g
    tf <- as.numeric(exp(logi) / (1 + exp(logi)) > 0.5)   #�N���X����
    
    #�������̌v�Z
    cr <- c(cr, sum(diag(table(y[index], tf))) / sum(table(y[index], tf)))
    print(cr)
  }
  CR <- c(CR, mean(cr))
  cat("��*'��')�Ă������[ ����lambda�� \n", 
      round(lambdaE[lam], 5), "�Ō��ʂ� \n", 
      CR , "���� \n")
}

####�œK��lambda��L1���������W�X�e�B�b�N��A���f���𐄒�####
#�œK�Ȑ������p�����[�^��I��
plot(lambdaE, CR, type="l", lwd=2)
round(lambdaE, 3); round(CR, 3)   #lambda�Ɛ������̕\��
(opt.lambda <- lambdaE[which.max(CR)])   #�œK��lambda��I��
b0 <- b0 + runif(length(b0), -0.2, 0.2)   #�����p�����[�^


##�A���S���Y���̐ݒ�
max.iter <- 30
iter <- 1
tol <- 1
diff <- 100
L1 <- 0

##L1���������W�X�e�B�b�N��A���f���𐄒�
while(diff >= tol & iter <= max.iter){
  for(i in 1:p){
    b1 <- b0[c(1, (i+1))]   #�œK������p�����[�^
    b2 <- b0[c(-1,-(i+1))]   #�Œ肷��p�����[�^
    x1 <- X[, i]   #�œK������p�����[�^�ɑΉ���������ϐ�
    x2 <- X[, -i]   #�Œ肷��p�����[�^�ɑΉ���������ϐ�
    
    #�����l��Nelder-Mead�@�Ō��肷��
    res1 <- optim(b1, fr, gr=NULL, b2, opt.lambda, x1, x2, y,
                  method="Nelder-Mead", hessian=FALSE, control=list(fnscale=-1))
    
    #�X�V�����p�����[�^�̊i�[
    b0[1] <- res1$par[1]
    b0[(i+1)] <- res1$par[2]
  }
  #�A���S���Y���̃p�����[�^�̍X�V
  iter <- iter+1
  LL <- res1$value
  diff <- abs(LL-L1)
  L1 <- LL
  
  #�ΐ��ޓx�̕\��
  print(diff)
}

##�֐����g�����ꍇ�Ő���
fit <- glmnet(x=X, y=y, lambda=opt.lambda, family="binomial", alpha=1)
b.func <- c(as.numeric(fit$a0), as.numeric(fit$beta))

##���肳�ꂽ���ʂƐ^�̃p�����[�^���r
res1$value   #�ΐ��ޓx
round(b0, 2)   #���肳�ꂽ�p�����[�^
round(b.func, 2)   #�֐��Ő��肳�ꂽ�p�����[�^
round(b, 2)   #�^�̃p�����[�^

