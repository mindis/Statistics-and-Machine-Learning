#####���E�̐��ݕϐ��𗣎U���Ԑ������f��#####
library(MASS)
library(nlme)
library(glmm)
library(survival)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(9483)

####�f�[�^�̔���####
hh <- 1000   #�T���v����
pt <- 36   #�ϑ�����
m <- 9

##ID�̐ݒ�
u.id <- rep(1:hh, rep(pt, hh))
t.id <- rep(1:pt, hh)
ID <- data.frame(no=1:(hh*pt), id=u.id, time=t.id)

####�����ϐ��̔���####
cont1 <- 3; bin1 <- 3; multi1 <- 4
X.cont <- matrix(rnorm(hh*cont1*pt), nrow=hh*pt, ncol=cont1)
X.bin <- matrix(0, nrow=hh*pt, ncol=bin1)
X.multi <- matrix(0, nrow=hh*pt, ncol=multi1)

#��l�����ϐ���ݒ�
for(i in 1:bin1){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(hh*pt, 1, p)
}

#���l�����ϐ���ݒ�
p <- runif(multi1)
x.multi <- t(rmultinom(hh, 1, p))

for(i in 1:hh){
  X.multi[ID$id==i, ] <- matrix(x.multi[i, ], nrow=pt, ncol=multi1, byrow=T)   #�璷�ȕϐ��͍폜
}
X.multi <- X.multi[, -which.min(colSums(X.multi))] #�璷�ȕϐ��͍폜

#�f�[�^������
X <- cbind(1, sqrt(ID$time), X.cont, X.bin, X.multi)


####���ݕϐ��Ɖ����ϐ��𔭐�####
##�����ϐ��̔���
#�p�����[�^�̐ݒ�
beta0 <- c(runif(1, -1.2, -1.0), runif(1, -0.07, 0.07), runif(cont1, 0, 0.6), runif(bin1+multi1-1, -1.0, 1.4))

#�w���̃��W�b�g�Ɗm���̌v�Z
logit <- X %*% beta0
P <- exp(logit)/(1+exp(logit))

#�w���L����񍀕��z���甭��
y <- rbinom(hh*pt, 1, P)

##���ݕϐ��̔���
#�p�����[�^�̐ݒ�
alpha0 <- runif(1, 1.8, 2.8)
alpha1 <- rep(0, hh*pt)
for(i in 1:hh) {alpha1[ID$id==i] <- rnorm(1, 0, 0.75)}
alpha2 <- beta0[(length(beta0)-(multi1-2)):length(beta0)] + runif(multi1-1, 0.7, 1.0)

#���E�̃��W�b�g�Ɗm���̌v�Z
logit0 <- alpha0 + alpha1 + X.multi %*% alpha2 
P0 <- exp(logit0)/(1 + exp(logit0))
summary(P0); hist(P0)

#���E�L����񍀕��z���甭��
z0 <- rbinom(hh*pt, 1, P0)

#���E���Ă�����A����ȍ~�����Ɨ��E
for(i in 1:hh){
  w <- z0[ID$id==i]
  if(min(w)==1){
    print("���E�Ȃ�")
    next
  } else {
    w[which.min(w):length(w)] <- 0
    z0[ID$id==i] <- w
  }
}   
surv_rate <- tapply(z0, ID$time, mean)   #������

##�ϑ�����鉞���ϐ��̐ݒ�
Y <- z0 * y
YZX <- round(data.frame(Y, y, P, z0, P0, ID, X), 3)


####�}���R�t�A�������e�J���@�Ő��ݕϐ����U���Ԑ������f���𐄒�####
##���W�X�e�B�b�N��A���f���̑ΐ��ޓx��ݒ�
loglike <- function(beta, y, X){
  logit <- X %*% beta    #���W�b�g�̌v�Z
  p <- exp(logit)/(1+exp(logit))   #�m���̌v�Z
  
  #�ΐ��ޓx�̌v�Z
  LLs <- y*log(p) + (1-y)*log(1-p)
  LL <- sum(LLs)
  return(LL)
}

##�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4
len <- length(z)
par <- ncol(X)

##���ݕϐ��̏����l�̐ݒ�
z <- rep(1, hh*pt)
for(i in 1:hh){
  print(i)
  y_ind <- Y[ID$id==i]
  
  if(sum(y_ind)==0){
    z[ID$id==i] <- 0
  } else {
    index <- max(subset(1:length(y_ind), y_ind==1))
    z[ID$id==i][index:pt] <- 0
    z[ID$id==i][index] <- 1
  }
}
z2 <- z   #�p�����[�^�X�V�p�̐��ݕϐ�

#ID�Ǝ��Ԃ̃C���f�b�N�X���X�g���쐬
time_list <- list()
id_list <- list()
for(i in 1:pt){
  time_list[[i]] <- subset(1:nrow(ID), ID$time==i)
}
for(i in 1:hh){
  id_list[[i]] <- subset(1:nrow(ID), ID$id==i)
}

len <- nrow(X)


##���O���z�̐ݒ�
betas <- rep(0, ncol(X))  #��A�W���̎��O���z
rootBi <- 0.01*diag(ncol(X))

##�����l�̐ݒ�
oldbeta <- rep(0, ncol(X))   #��A�W���̕��ς̏����l
r <- tapply(z2, ID$time, mean)   #���E���̎��O���z�̏����l

##�����_���E�H�[�N�̕��U��ݒ�
#�ΐ��ޓx���ő剻
b0 <- c(rep(0, ncol(X)))   #�����p�����[�^�̐ݒ�
res <- optim(b0, loglike, gr=NULL, y=Y[z==1], X=X[z==1, ], method="BFGS", hessian=TRUE, control=list(fnscale=-1))
rw <- solve(-res$hessian)   #�����_���E�H�[�N�̕��U

##�p�����[�^�̕ۑ��p�z��
BETA <- matrix(0, nrow=R/keep, ncol=ncol(X))
Z <- matrix(0, nrow=R/keep, ncol=nrow(X))
mix_rate <- matrix(0, nrow=R/keep, ncol=pt)
z_time <- list()
z_surv <- list()
index_surv <- list()


####MCMC�@�Ő��ݕϐ����U���Ԑ������f���𐄒�####
for(rp in 1:R){

  ##���g���|���X�w�C�X�e�B���O�A���S���Y���ŉ�A�W��beta���T���v�����O
  #z=1�̕ϐ������o��
  index_z1 <- subset(1:len, z==1)
  XZ1 <- X[index_z1, ]
  YZ1 <- Y[index_z1]
  
  #MH�@��beta���T���v�����O
  betad <- oldbeta
  betan <- betad + 0.5*mvrnorm(1, rep(0, par), rw)
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew <- loglike(betan, y=YZ1, X=XZ1)
  logold <- loglike(betad, y=YZ1, X=XZ1)
  logpnew <- lndMvn(betan, betas, rootBi)
  logpold <- lndMvn(betad, betas, rootBi)
  
  #�T���v�����O���ꂽbeta���̑����邩����
  #MH�T���v�����O
  alpha <- min(1, exp(lognew + logpnew - logold - logpold))
  if(alpha == "NAN") alpha <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V����beta���̑�
  if(u < alpha){
    oldbeta <- betan
    logl <- lognew
    
    #�����łȂ��Ȃ�beta���X�V���Ȃ�
  } else {
    oldbeta <- betad
  }
  
  ##���ݕϐ�z���T���v�����O
  index_z2 <- subset(1:len, z2==0)
  XZ2 <- X[index_z2, ]
  ID_z <- ID[index_z2, ]
  r_rate <- r[ID_z$time]
  
  #���W�b�g�Ɗm�����v�Z
  logit <- XZ2 %*% oldbeta
  Pr <- as.numeric(exp(logit)/(1+exp(logit)))

  #���_t�܂łɍw�����Ă��Ȃ��m�����v�Z
  pr_prod <- unlist(tapply(Pr, ID_z$id, cumprod))

  #���E���Ă���m�����v�Z
  z_rate <- ((1-r_rate) * (1-pr_prod))/((r_rate * pr_prod) + ((1-r_rate) * (1-pr_prod)))
  z1 <- rbinom(length(z_rate), 1, z_rate)

  #���E���Ă�����A����ȍ~�͂��ׂė��E������
  for(i in 1:hh){
    index_id <- subset(1:nrow(ID_z), ID_z$id==i)
    z_ind <- z1[index_id]
    index <- which.max(z_ind)

    if(max(z_ind)==0 | length(z_ind)==0){
      next
    } else {
      z_ind[index:length(z_ind)] <- 1
      z1[index_id] <- z_ind
    }
  }
  
  ##���ݕϐ�z�Ə����t�����������X�V
  #���ݕϐ�z���X�V
  z[index_z2] <- 1-z1

  #�J�v�����}�C���[�@�ŏ����t��������r���X�V
  z_time[[1]] <- z[time_list[[1]]]
  z_surv[[1]] <- z_time[[1]]
  index_surv[[1]] <- subset(1:length(z_time[[1]]), z_time[[1]]==1)
  r[1] <- mean(z_surv[[1]])
 
  for(i in 2:pt){ 
    z_time[[i]] <- z[time_list[[i]]]
    z_surv[[i]] <- z_time[[i]][index_surv[[i-1]]]
    index_surv[[i]] <- subset(1:length(z_time[[i]]), z_time[[i]]==1)
    r[i] <- mean(z_surv[[i]])
  }
  
  
  ##�T���v�����O���ʂ�ۑ�
  if(rp%%keep==0){
    mkeep <- rp/keep
    surv_rate1 <- tapply(z, ID$time, mean)
    BETA[mkeep, ] <- oldbeta
    Z[mkeep, ] <- z
    mix_rate[mkeep, ] <- surv_rate1
    
    print(rp)
    print(logl)
    print(round(rbind(oldbeta, beta0), 3))
    print(round(rbind(surv_rate1, surv_rate2=surv_rate), 3))
  }
}

####���茋�ʂ̊m�F�Ɨv��####
burnin <- 5000/keep   #�o�[���C������

logit <- X %*% oldbeta
Pr <- exp(logit)/(1+exp(logit))
res <- round(cbind(ID, z, z0, z2, Y, Pr, P), 3)


##�T���v�����O���ʂ̃v���b�g
matplot(BETA[, 1:4], type="l")
matplot(BETA[, 5:8], type="l")
matplot(BETA[, 9:11], type="l")
matplot(mix_rate[, 1:5], type="l")
matplot(mix_rate[, 6:10], type="l")
matplot(mix_rate[, 11:15], type="l")
matplot(mix_rate[, 16:20], type="l")

mean(z[ID$time==pt])
mean(z0[ID$time==pt])

matplot(BETA[, 9:10], type="l")
matplot(mix_rate[, 33:36], type="l")
ncol(BETA)