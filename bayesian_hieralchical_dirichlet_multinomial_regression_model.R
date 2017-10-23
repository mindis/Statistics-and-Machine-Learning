#####�K�w�x�C�Y�f�B�N����������A���f��#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(extraDistr)
library(matrixStats)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(6812)

####�f�[�^�̔���####
hh <- 1000   #���[�U�[��
pt <- 10   #�ϑ�����
hhpt <- hh*pt   #���T���v����
select <- 10   #�I������
st <- 10   #��I����

##ID�̐ݒ�
id <- rep(1:hh, rep(pt, hh))
time <- rep(1:pt, hh)
ID <- data.frame(no=1:hhpt, id, time)

#�x�N�g���^ID�̐ݒ�
u.id <- rep(1:hh, rep(pt*select, hh))
i.id <- rep(1:select, pt*hh)
time <- rep(rep(1:pt, rep(select, pt)), hh)
ID_vec <- data.frame(no=1:length(u.id), id=u.id, time=time, item=i.id) 

####�����ϐ��̔���####
##���U�I�����f���̐����ϐ��̔���
##�����t���̐����ϐ��̔���
X1.cont <- matrix(rnorm(hhpt*select*2, 0, 1), nrow=hhpt*select, ncol=2)

X1.bin <- matrix(0, nrow=hhpt*select, ncol=3)
for(i in 1:3){
  X1.bin[, i]  <- rbinom(hhpt*select, 1, runif(1, 0.35, 0.6))
}

X1.multi <- t(rmultinom(hhpt*select, 1, c(0.2, 0.2, 0.3, 0.3)))[, -4]

#�ؕЂ̐ݒ�
Pop0 <- matrix(diag(1, select), nrow=hhpt*select, ncol=select, byrow=T)
Pop <- Pop0[, -st]

#�f�[�^�̌���
X <- data.frame(bp=Pop, c=X1.cont, b=X1.bin, m=X1.multi)
XM <- as.matrix(X)
round(XM, 2)


##�����^�̐����ϐ��̔���
X2.cont <- rnorm(hh, 0, 1)
X2.pois <- rpois(hh, 2)
X2.bin <- matrix(rbinom(hh*2, 1, 0.4), nrow=hh, ncol=2)

#�f�[�^�̌���
Z <- data.frame(i=1, c=X2.cont, p=X2.pois, b=X2.bin)
ZX <- as.matrix(Z)


####�����ϐ��̔���####
##�w�����ʂ𔭐�
freq <- rep(0, hhpt)
for(i in 1:hh){
  par <- rgamma(hh, 25, 1.25)
  freq[ID$id==i] <- rpois(pt, par)   #�w����
}

##�K�w��A���f������p�����[�^�𔭐�
#��A�p�����[�^��ݒ�
gamma00 <- c(runif(select-1, -0.3, 0.4), runif(ncol(X1.cont), 0, 0.25), runif(ncol(X1.bin)+ncol(X1.multi), -0.25, 0.25))
gamma01 <- runif(ncol(XM), 0, 0.25)
gamma02 <- runif(ncol(XM), -0.15, 0.15)
gamma03 <- matrix(runif(ncol(XM*2), -0.2, 0.25), nrow=2, ncol=ncol(XM))
gamma0 <- rbind(gamma00, gamma01, gamma02, gamma03)
rownames(gamma0) <- rep("", ncol(ZX))

#���U�p�����[�^��ݒ�
Cov0 <- diag(runif(ncol(XM), 0.15, 0.4))

#���ϗʐ��K���z����f�B�N����������A���f���̉�A�p�����[�^�𔭐�
beta_mu <- ZX %*% gamma0   #���ύ\��
er <- mvrnorm(hh, rep(0, ncol(XM)), Cov0)   #�ϗʌ��ʂ̌덷
beta0 <- beta_mu + er   #�f�B�N����������A���f���̃p�����[�^

##�f�B�N�������z����m���𔭐�
Alpha_m <- matrix(0, nrow=hhpt, ncol=select)
Prob <- matrix(0, nrow=hhpt, ncol=select)

for(i in 1:hh){
  Alpha_m[ID$id==i, ] <- matrix(exp(XM[ID$id==i, ] %*% beta0[i, ]), nrow=pt, ncol=select, byrow=T)
  Prob[ID$id==i, ] <- t(apply(Alpha_m[ID$id==i, ], 1, function(x) rdirichlet(1, x)))
}

##�������z���牞���ϐ��𔭐�
Y <- t(apply(cbind(freq, Prob), 1, function(x) rmultinom(1, x[1], x[-1])))
y <- as.numeric(t(Y))
Pr <- Y / matrix(rowSums(Y), nrow=hhpt, ncol=select)

####�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y�f�B�N����������A���f���𐄒�####
##�f�B�N����������A���f���̑ΐ��ޓx�֐����`
loglike <- function(theta, y, X, hh, select){

  #�p�����[�^��ݒ�
  beta <- theta
  
  #�f�B�N�������z�̃p�����[�^�̕��ύ\����ݒ�
  alpha <- matrix(exp(X %*% beta), nrow=hh, ncol=select, byrow=T)
  
  #�f�B�N�������z�̑ΐ��ޓx��a
  LLi <- ddirichlet(y, alpha, log=TRUE)
  LLi[is.infinite(LLi)] <- 0
  LL <- sum(LLi)
  return(LL)
}

##�f�B�N����������A���f���̖ޓx�֐����`
fr <- function(theta, y, X, hh, select){
  #�p�����[�^��ݒ�
  beta <- theta
  
  
  
  #�f�B�N�������z�̃p�����[�^�̕��ύ\����ݒ�
  alpha <- matrix(exp(X %*% beta), nrow=hh, ncol=select, byrow=T)
  
  #�f�B�N�������z�̑ΐ��ޓx��a
  LLi <- ddirichlet(y, alpha)
  LL <- prod(LLi)
  return(LL)
}

##�A���S���Y���̐ݒ�
R <- 20000   #�T���v�����O��
keep <- 4   #2���1��̊����ŃT���v�����O���ʂ��i�[
iter <- 0

##���O���z�̐ݒ�
Deltabar <- matrix(0, nrow=ncol(ZX), ncol=ncol(XM))
Adelta <- 0.01 * diag(rep(1, ncol(ZX)))   #�K�w���f���̉�A�W���̕��U�̎��O���z
nu <- (ncol(XM)) + select   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(ncol(XM))

##�T���v�����O���ʂ̕ۑ��p�z��
PROB <- matrix(0, hhpt, select)
BETA <- array(0, dim=c(hh, ncol(XM), R/keep))
GAMMA <- matrix(0, nrow=R/keep, ncol=length(gamma0))
SIGMA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
gc(); gc()

#���p���Ƒΐ��ޓx�̕ۑ��p�z��
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##�C���f�b�N�X���쐬
index_vec <- matrix(1:nrow(ID_vec), nrow=hh, ncol=pt*select, byrow=T)
index_id <- matrix(1:nrow(ID), nrow=hh, ncol=pt, byrow=T)

##�p�����[�^�̊i�[�p�z��
lognew <- rep(0, hh)
logold <- rep(0, hh)
logpnew <- rep(0, hh)
logpold <- rep(0, hh) 

##�����l�̐ݒ�
alpha_mu <- matrix(0, nrow=hh, ncol=select)
for(j in 1:ncol(Y)){
  alpha_mu[, j] <- tapply(Y[, j] + 1, ID$id, sum)
}
freq_mu <- tapply(rowSums(Y) + 1, ID$id, mean)
Pr <- alpha_mu/rowSums(alpha_mu)


#�f�B�N�������f���̃p�����[�^�̏����l
#�����ϐ��̐ݒ�
Pr <- (Y + 1)/rowSums(Y + 1)

#�p�����[�^�̏����l
phi0 <- runif(select-1, -0.5, 0.5)
phi1 <- c(runif(2, 0, 0.5), runif(3, -0.4, 0.4), runif(3, -0.3, 0.3))
phi <- c(phi0, phi1)

#�ΐ��ޓx���ő剻����
res <- optim(phi, loglike, gr=NULL, Pr, XM, hhpt, select, method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))
beta <- res$par
LL <- res$value
rw <- -diag(diag(solve(res$hessian)))
oldbeta <- matrix(beta, nrow=hh, ncol=ncol(XM), byrow=T) + matrix(rnorm(hh*ncol(XM), 0, 0.15), nrow=hh, ncol=ncol(XM))

#�p�����[�^alpha��ݒ�
oldalpha <- matrix(0, nrow=nrow(Y), ncol=select)
for(i in 1:hh){
  oldalpha[index_id[i, ], ] <- matrix(exp(XM[index_vec[i, ], ] %*% oldbeta[i, ]), nrow=pt, ncol=select, byrow=T)
}

#�K�w���f���̃p�����[�^�̏����l
oldcov <- diag(rep(0.2, ncol(XM)))
inv_cov <- solve(oldcov)
oldgamma <- matrix(runif(ncol(XM)*ncol(ZX), -0.2, 0.2), nrow=ncol(ZX), ncol=ncol(XM))

####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O
for(rp in 1:R){

  ##�����f�B�N�������z����p�����[�^Prob���T���v�����O
  d_par <- Y + oldalpha
  oldprob <- t(apply(d_par, 1, function(x) rdirichlet(1, x)))
  
  ##MH�@�Ń��[�U�[���ƂɃf�B�N����������A���f���̃p�����[�^beta���T���v�����O
  #�V�����p�����[�^���T���v�����O
  betad <- oldbeta
  betan <- betad + 2 * mvrnorm(hh, rep(0, ncol(XM)), rw)
  
  #�p�����[�^�̎��O���z�ƌ덷
  mu <- ZX %*% oldgamma
  er_new <- betan - mu 
  er_old <- betad - mu
  
  for(i in 1:hh){
    lognew[i] <- loglike(betan[i, ], oldprob[index_id[i, ], ], XM[index_vec[i, ], ], pt, select)
    logold[i] <- loglike(betad[i, ], oldprob[index_id[i, ], ], XM[index_vec[i, ], ], pt, select)
    logpnew[i] <- -0.5 * (er_new[i, ] %*% inv_cov %*% er_new[i, ])
    logpold[i] <- -0.5 * (er_old[i, ] %*% inv_cov %*% er_old[i, ])
  }

  #���g���|���X�w�C�X�e�B���O�@�Ńp�����[�^�̍̑�������
  rand <- matrix(runif(hh), nrow=hh, ncol=ncol(betan))   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew + logpnew - logold - logpold)   #�̑𗦂��v�Z
  omega <- matrix(ifelse(LLind_diff > 1, 1, LLind_diff), nrow=hh, ncol=ncol(betan))
  
  #omega�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  oldbeta <- ifelse(omega > rand, betan, betad)   #omega��rand�������Ă�����̑�
  
  #alpha���X�V
  for(i in 1:hh){
    oldalpha[index_id[i, ], ] <- matrix(exp(XM[index_vec[i, ], ] %*% oldbeta[i, ]), nrow=pt, ncol=select, byrow=T)
  }   
  
  ##���ϗʉ�A���f���ɂ��K�w���f���̃T���v�����O
  out <- rmultireg(Y=oldbeta, X=ZX, Bbar=Deltabar, A=Adelta, nu=nu, V=V)
  oldgamma <- out$B
  oldcov <- diag(diag(out$Sigma))
  inv_cov <- solve(oldcov)
  
  solve(t(ZX) %*% ZX) %*% t(ZX) 
  sum(is.na(oldbeta))
  oldbeta
  t(ZX) %*% oldbeta
  
  ##�T���v�����O���ʂ�ۑ�
  if(rp%%keep==0){
    cat("��*'��')�� <�������[ ������Ƃ܂��Ă�", paste(rp/R*100, "��"), fill=TRUE)
    mkeep <- rp/keep
    PROB <- PROB + oldprob
    BETA[, , mkeep] <- oldbeta 
    GAMMA[mkeep, ] <- as.numeric(oldgamma)
    SIGMA[mkeep, ] <- diag(oldcov)
    print(rp)
    print(c(sum(ifelse(omega[, 1] > rand[, 1], lognew, logold)), LL))
    print(round(rbind(oldgamma, gamma0), 3))
  }
}

matplot(GAMMA[, 1:5], type="l")
matplot(t(BETA[1:5, 1, ]), type="l")
t(BETA[1:5, 1, 1:1000])

GAMMA

round(cbind(oldprob, Prob), 3)
