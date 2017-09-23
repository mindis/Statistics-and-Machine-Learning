#####���֍\���̂���K�w�x�C�YRF���f��#####
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

#set.seed(431278)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 1000   #�T���v����
pt <- rpois(hh, 15.0)
pt <- ifelse(pt==0, 1, pt)
hhpt <- sum(pt)
dt <- 150

##ID�̐ݒ�
id <- rep(1:hh, pt)
time <- c()
for(i in 1:hh){time <- c(time, 1:pt[i])}
ID <- data.frame(no=1:hhpt, id=id, time=time)

####�����ϐ��̔���####
##�K�w���f���̐����ϐ��̔���
k1 <- 5
cont1 <- 2; bin1 <- 3
X1 <- matrix(0, nrow=hh, ncol=k1)
for(i in 1:hh){
  X1[i, 1:cont1] <- rnorm(2, 0, 1)
  for(j in 1:bin1){
    X1[i, (cont1+j)] <- rbinom(1, 1, runif(1, 0.4, 0.6))
  }
}

##���ς̐����ϐ�
Price <- scale(rbeta(hhpt, 7.5, 1.5))   #����̕��ϒl����
Disc <- scale(rbeta(hhpt, 2.0, 6.5))   #����̒l�����Q�[���̊���
Prom <- rbinom(hhpt, 1, 0.4)   #�v�����[�V�����L��

#�V��Q�[���{���ƃW�������̔���
g <- 6
Genre <- matrix(0, nrow=hhpt, ncol=g)
pr_gen <- runif(g, 0.3, 1.0)

#�T���v�����Ƃɑ������z����Q�[���W�������𔭐�������
for(i in 1:hhpt){
  new <- rpois(1, 6.5)
  if(new==0){
    next
  } else {
    Genre[i, ] <- t(rmultinom(1, new, pr_gen))
  }
}

##�f�[�^������
ZX <- X1
colnames(ZX) <- c("cont1", "cont2", "bin1", "bin2", "bin3")

XM <- cbind(1, Price, Disc, Prom, Genre)
colnames(XM) <- c("�ؕ�", "Price", "Disc", "Prom", "genre1", "genre2", "genre3", "genre4", "genre5", "genre6")

##�ϗʌ��ʂ̃f�U�C���s��̐ݒ�
Z <- matrix(0, nrow=hhpt, ncol=hh)
for(i in 1:hh){
  Z[ID$id==i, i] <- 1 
}

####�����ϐ��̔���####
for(i in 1:1000){
  ##�p�����[�^�̐ݒ�
  #�K�w���f���̃p�����[�^�̐ݒ�
  theta01 <- c(runif(cont1, 0, 0.4), runif(bin1, -0.4, 0.8))
  theta02 <- c(runif(cont1, 0, 0.4), runif(bin1, -0.3, 0.2))
  theta0 <- cbind(theta01, theta02)
  
  ##�ϗʌ��ʂ̃p�����[�^�𑽕ϗʐ��K���z���甭��
  #���U�����U�p�����[�^��ݒ�
  Cov0 <- matrix(c(0.6, -0.35, -0.35, 0.5), nrow=2, ncol=2)
  
  #�ϗʌ��ʂ̃p�����[�^�𔭐�
  random <- mvrnorm(hh, rep(0, 2), Cov0)
  theta <- ZX %*% theta0 + random
 
  ##�Œ���ʂ̃p�����[�^�𔭐�
  #�������f���̃p�����[�^�̐ݒ�
  alpha0 <- runif(1, 0.8, 1.8)   #�`��p�����[�^
  beta0 <- c(runif(1, 0, 5.0), runif(1, 0.1, 0.4), runif(1, -0.4, -0.1), runif(1, -0.7, -0.2), runif(g, -0.4, 0))
  
  #�p�x���f���̃p�����[�^�̐ݒ�
  gamma0 <- c(runif(1, 0, 0.5), runif(1, -0.25, 0), runif(1, 0, 0.25), runif(1, 0, 0.35), runif(g, 0, 0.2)) 
  
  ##�������f���ƕp�x���f���̉����ϐ��𔭐�
  #���ύ\����ݒ�
  scale0 <- exp(XM %*% beta0 + Z %*% theta[, 1])
  lambda0 <- exp(XM %*% gamma0 + Z %*% theta[, 2])
  
  #���C�u�����z�ƃ|�A�\�����z���w���Ԋu�ƍw�����𔭐�
  y1 <- rweibull(hhpt, alpha0, scale0)
  y2 <- rpois(hhpt, lambda0)
  
  print(round(c(min(y1), max(y1), min(y2), max(y2)), 3))
  if(min(y1) > 0.01 & max(y2) < 50) break
}

Y <- cbind(y1, y2)

#���������������ϐ�������
hist(y1[y1 < 100], xlab="�w���Ԋu", main="�Q�[���X�ւ̖K��Ԋu", col="grey")
hist(y2[y2 < 20], xlab="�w����", main="�Q�[���̍w����", col="grey")

####�ł��؂�̐ݒ�####
##���[�U�[���Ƃ�T = 150�܂Ŋϑ�
#�ϐ��̊i�[�p���X�g
ID.list <- list()
y1.list <- list()
y2.list <- list()
X.list <- list()
Z.list <- list()
z.list <- list()

#�l���Ƃɑł��؂�ϐ���ݒ�
for(i in 1:hh){
  print(i)
  y1_ind <- Y[ID$id==i, 1]
  y2_ind <- Y[ID$id==i, 2]
  z <- rep(0, length(y1_ind))
  c_sum <- cumsum(y1_ind)
  
  #�ݐώ��Ԃ�dt�ȏ�̃C�x���g�͑ł��؂�
  index1 <- subset(1:length(c_sum), c_sum <= dt)
  
  if(max(c_sum) <= dt){
    index2 <- index1
  } else {
    index2 <- c(index1, length(index1)+1)
  }
  
  #�����ϐ��̑ł��؂��ݒ�
  if(max(c_sum) > dt & length(index1) > 0){
    print(1)
    y_vec <- c(y1_ind[index1], dt-c_sum[length(index1)])
    z[length(y_vec)] <- 1
  } else if(max(c_sum) > dt & length(index1)==0) {
    print(2)
    y_vec <- dt
    z <- 1
  } else {
    print(3)
    y_vec <- y1_ind[index2]
  }

  #�ł��؂�ꂽ�ϐ����i�[
  y1.list[[i]] <- y_vec[index2]
  y2.list[[i]] <- y2_ind[index2]
  ID.list[[i]] <- ID[ID$id==i, ][index2, ]
  if(class(XM[ID$id==i, ])=="numeric"){
    X.list[[i]] <- XM[ID$id==i, ]
    Z.list[[i]] <- Z[ID$id==i, ]
  } else {
    X.list[[i]] <- XM[ID$id==i, ][index2, ]
    Z.list[[i]] <- Z[ID$id==i, ][index2, ]
  }
  z.list[[i]] <- z[index2]
}

#���X�g���x�N�g�����邢�͍s��
y1 <- unlist(y1.list)
y2 <- unlist(y2.list) 
ID <- do.call(rbind, ID.list)
XM <- do.call(rbind, X.list)
Z <- do.call(rbind, Z.list)
z <- 1-unlist(z.list)
hhpt <- nrow(ID)

#�f�[�^�̊m�F�Ɖ���
round(cbind(ID, y1, y2, z, XM), 3)
hist(y1, xlab="�w���Ԋu", main="�Q�[���X�ւ̖K��Ԋu", col="grey")
hist(y2, xlab="�w����", main="�Q�[���̍w����", col="grey")

cbind(ID, y1, z)

####�}���R�t�A�������e�J�����@�ŊK�w�x�C�YRF���f���𐄒�####
##RF���f���̑ΐ��ޓx�֐��̐ݒ�
loglike <- function(alpha, beta, gamma, theta, y1, y2, z, X, Z){
  
  #���ύ\���̐ݒ�
  scale <- exp(X %*% beta + Z %*% theta[, 1])
  lambda <- exp(X %*% gamma + Z %*% theta[, 2])
  
  #�ΐ��ޓx���v�Z
  LL1 <- z*(log(scale)+log(alpha)+(alpha-1)*log(y1)) - scale*y1^alpha   #���C�u�����f���̑ΐ��ޓx
  LL2 <- y2*log(lambda)-lambda - lfactorial(y2)   #�|�A�\�����f���̑ΐ��ޓx
  LLi <- LL1 + LL2
  LL <- sum(LL1 + LL2)
  LL_val <- list(LLi=LLi, LL=LL)
  return(LL_val)
}

##�ΐ��ޓx�֐��̍ő剻�p�֐��̐ݒ�
rf_model <- function(theta, y1, y2, z, X, index_par1, index_par2){
  
  #�p�����[�^�̐ݒ�
  alpha <- theta[1]
  beta <- theta[index_par1]
  gamma <- theta[index_par2]
  
  #���ύ\���̐ݒ�
  scale <- exp(X %*% beta)
  lambda <- exp(X %*% gamma)
  
  #�ΐ��ޓx���v�Z
  LL1 <- z*(log(scale)+log(alpha)+(alpha-1)*log(y1)) - scale*y1^alpha   #���C�u�����f���̑ΐ��ޓx
  LL2 <- y2*log(lambda)-lambda - lfactorial(y2)   #�|�A�\�����f���̑ΐ��ޓx
  LL <- sum(LL1) + sum(LL2)
  return(LL)
}

##�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4
len <- nrow(X)
par <- ncol(X)

##���O���z�̐ݒ�
#�`��p�����[�^�̎��O���z
alpha_mu <- 0
alpha_sigma <- 2.5

#�Œ���ʂ̃p�����[�^�̎��O���z
betas <- rep(0, ncol(XM)*2)   #��A�W���̕��ς̎��O���z
sigma <- diag(rep(0.01), ncol(XM)*2)   #��A�W���̎��O���z�̕��U

#�K�w���f���̎��O���z
Deltabar <- matrix(0, nrow=ncol(ZX), ncol=2)
Adelta <- 0.01 * diag(ncol(ZX))
nu <- 2 + ncol(ZX)
V <- nu * diag(2)

##�T���v�����O���ʂ̕ۑ��p�z��
ALPHA <- rep(0, R/keep)
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XM)*2)
THETA0 <- matrix(0, nrow=R/keep, ncol=ncol(ZX)*2)
THETA1 <- matrix(0, nrow=R/keep, ncol=hhpt)
THETA2 <- matrix(0, nrow=R/keep, ncol=hhpt)
SIGMA <- matrix(0, nrow=R/keep, ncol=2^2)


##���p���Ƒΐ��ޓx�̕ۑ��p�z��
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##�����l�̐ݒ�
#�K�w���f���̃p�����[�^�̐ݒ�
oldtheta01 <- -c(runif(cont1, 0, 0.5), runif(bin1, -0.8, 0.8))
oldtheta02 <- c(runif(cont1, 0, 0.6), runif(bin1, -0.3, 0.3))
oldtheta0 <- cbind(oldtheta01, oldtheta02)

##�ϗʌ��ʂ̃p�����[�^�𑽕ϗʐ��K���z���甭��
#���U�����U�p�����[�^��ݒ�
oldcov <- matrix(c(0.5, -0.2, -0.2, 0.5), nrow=2, ncol=2)
cov_inv <- solve(oldcov)

#�ϗʌ��ʂ̃p�����[�^�𔭐�
oldrandom <- mvrnorm(hh, rep(0, 2), Cov0)
theta_mu <- ZX %*% oldtheta0
oldtheta <- ZX %*% oldtheta0 + oldrandom

##�Œ���ʂ̃p�����[�^�𔭐�
##�p�����[�^�̏����l�ƃ����_���E�H�[�N�̕��U��ݒ肷�邽�߂�RF���f���̑ΐ��ޓx���ő剻
#�����l�̐ݒ�
oldalpha <- runif(1, 0.8, 1.8)   #�`��p�����[�^
oldbeta <- -c(runif(1, 1, 5.5), runif(1, 0.4, 1.0), runif(1, -0.8, -0.3), runif(1, -0.7, -0.2), runif(g, -0.4, 0))
oldgamma <- c(runif(1, 0, 0.6), runif(1, -0.5, -0.2), runif(1, 0, 0.3), runif(1, 0, 0.35), runif(g, 0, 0.2)) 
par <- c(oldalpha, oldbeta, oldgamma)

#�p�����[�^�̃C���f�b�N�X
index1 <- 2:(1+ncol(XM))
index2 <- (2+ncol(XM)):length(par)

#���j���[�g���@�őΐ��ޓx���ő剻
res <- try(optim(par, rf_model, y1=y1, y2=y2, z=z, X=XM, index_par1=index1, index_par2=index2, method="BFGS", 
                 hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)

#���茋��
par <- res$par

#�����_���E�H�[�N�̕��U��ݒ�
rw <- diag(solve(-res$hessian))
rw1 <- rw[1]
rw2 <- diag(rw[2:length(par)])

#�������f���̃p�����[�^�̏����l��ݒ�
oldalpha <- res$par[1]   #�`��p�����[�^�̏����l
oldbeta <- par[index1]

#�p�x���f���̃p�����[�^�̐ݒ�
oldgamma <- par[index2]

#�p�����[�^������
olddelta <- c(oldbeta, oldgamma)
index_par1 <- 1:length(oldbeta)
index_par2 <- (length(oldbeta)+1):length(olddelta)


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##MH�@�ŌŒ���ʃp�����[�^���T���v�����O
  deltad <- olddelta 
  deltan <- deltad + 0.15 * mvrnorm(1, rep(0, length(deltad)), rw2)
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew1 <- loglike(oldalpha, deltan[index_par1], deltan[index_par2], oldtheta, y1, y2, z, XM, Z)$LL
  logold1 <- loglike(oldalpha, deltad[index_par1], deltad[index_par2], oldtheta, y1, y2, z, XM, Z)$LL
  logpnew1 <- lndMvn(deltan, betas, sigma)
  logpold1 <- lndMvn(deltad, betas, sigma)
  
  #MH�T���v�����O
  alpha1 <- min(1, exp(lognew1 + logpnew1 - logold1 - logpold1))
  if(alpha1 == "NAN") alpha2 <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V�����Œ����beta���̑�
  if(u < alpha1){
    olddelta <- deltan
    logl <- lognew1
    
    #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
  } else {
    olddelta <- deltad
    logl <- logold1
  }  
  
  ##MH�@�Ō`��p�����[�^���T���v�����O
  alphad <- abs(oldalpha)
  alphan <- abs(alphad + rnorm(1, 0, sqrt(rw1)))

  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew2 <- loglike(alphan, olddelta[index_par1], olddelta[index_par2], oldtheta, y1, y2, z, XM, Z)$LL
  logold2 <- loglike(alphad, olddelta[index_par1], olddelta[index_par2], oldtheta, y1, y2, z, XM, Z)$LL
  logpnew2 <- -1/2 * alpha_sigma^(-1) * (log(alphan) - alpha_mu)^2
  logpold2 <- -1/2 * alpha_sigma^(-1) * (log(alphad) - alpha_mu)^2
  
  #MH�T���v�����O
  alpha2 <- min(1, exp(lognew2 + logpnew2 - logold2 - logpold2))
  if(alpha2 == "NAN") alpha2 <- -1

  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V�����Œ����beta���̑�
  if(u < alpha2){
    oldalpha <- alphan
    
    #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
  } else {
    oldalpha <- alphad
  }  
  
  
  ##MH�@�ŕϗʌ��ʂ��T���v�����O
  oldthetad <- oldtheta
  oldthetan <- oldthetad + 0.5 * mvrnorm(hh, rep(0, 2), oldcov)
  
  #���O���z�̌덷���v�Z
  er_new <- oldthetan - theta_mu
  er_old <- oldthetad - theta_mu
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew3 <- loglike(oldalpha, olddelta[index_par1], olddelta[index_par2], oldthetan, y1, y2, z, XM, Z)$LLi
  logold3 <- loglike(oldalpha, olddelta[index_par1], olddelta[index_par2], oldthetad, y1, y2, z, XM, Z)$LLi
  logpnew3 <- apply(er_new, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  logpold3 <- apply(er_old, 1, function(x) -0.5 * (x %*% cov_inv %*% x))
  
  #ID�ʂɑΐ��ޓx�̘a�����
  logl.ind <- as.matrix(data.frame(logl1=lognew3, logl2=logold3, id=ID$id) %>%
                          dplyr::group_by(id) %>%
                          dplyr::summarize(new=sum(logl1), old=sum(logl2)))[, 2:3]
  
  ##MH�T���v�����O
  #�T���v�����O���̑����邩�ǂ���������
  rand <- runif(hh)   #��l�������痐���𔭐�
  LLind.diff <- exp(logl.ind[, 1] + logpnew3 - logl.ind[, 2] - logpold3)   #�̑𗦂��v�Z
  alpha <- ifelse(LLind.diff > 1, 1, LLind.diff)
  alpha <- matrix(alpha, nrow=hh, ncol=2)
  
  #alpha�Ɋ�Â�beta���̑�
  oldtheta.r <- ifelse(alpha > rand, oldthetan, oldthetad)   #alpha��rand�������Ă�����̑�
  adopt <- sum(oldtheta[, 1]!=oldtheta.r[, 1])/hh   #�̑�
  oldtheta <- oldtheta.r   #�p�����[�^���X�V
  
  
  ##���ϗʉ�A���f���ɂ��K�w���f���̃T���v�����O
  out <- rmultireg(Y=oldtheta, X=ZX, Bbar=Deltabar, A=Adelta, nu=nu, V=V)
  oldtheta0 <- out$B
  oldcov <- out$Sigma
  
  #�K�w���f���̃p�����[�^���X�V
  cov_inv <- solve(oldcov)
  theta_mu <- ZX %*% oldtheta0
  

  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    ALPHA[mkeep] <- oldalpha
    BETA[mkeep, ] <- olddelta
    THETA0[mkeep, ] <- as.numeric(oldtheta0)
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    colnames(oldtheta0) <- c("theta01", "theta02")
   
    print(rp)
    print(round(c(logl, res$value), 2))
    print(round(c(oldalpha, alpha0), 3))
    print(round(rbind(olddelta, delta_t=c(beta0, gamma0), ML_par=par[2:length(par)]), 3))
    print(round(rbind(theta=t(oldtheta0), t(theta0)), 3))
    print(round(cbind(oldcov, Cov0), 3))
    print(round(cbind(cov2cor(oldcov), cov2cor(Cov0)), 3))
    print(round(adopt, 3))
  }
}

matplot(BETA[, 6:10], type="l")
matplot(THETA0[, 1:5], type="l")
matplot(SIGMA, type="l")
XM