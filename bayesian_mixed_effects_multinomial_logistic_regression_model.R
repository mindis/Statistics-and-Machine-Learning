#####�ϗʌ��ʍ������W�b�g���f��#####
library(MASS)
library(mlogit)
library(MCMCpack)
library(bayesm)
library(caret)
library(reshape2)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)


####�f�[�^�̔���####
#set.seed(8437)
##�f�[�^�̐ݒ�
hh <- 250   #�T���v����
pt <- rpois(hh, 15); pt <- ifelse(pt==0, 1, pt)   #�w���@��(�w���@���0�Ȃ�1�ɒu������)
hhpt <- sum(pt)
choise <- 5   #�I���\��
st <- 5   #��u�����h
k <- 5   #�����ϐ��̐�

##�����ϐ��̔���
#ID�̐ݒ�
id <- rep(1:hh, pt)
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}
ID <- cbind(no=1:hhpt, id, t)

#�ʏ퉿�i�̔���
PRICE <- matrix(runif(hhpt*choise, 0.6, 1), nrow=hhpt, ncol=choise, byrow=T)   

#�f�B�X�J�E���g���̔���
DISC <- matrix(runif(hhpt*choise, 0, 0.5), nrow=hhpt, ncol=choise, byrow=T)

#���ʒ�̔���
DISP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hhpt, 1, r)
}

#���ʃL�����y�[���̔���
CAMP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hhpt, 1, r)
}

#�J�e�S���[���C�����e�B
ROYL <- matrix(runif(hhpt, 0, 1), nrow=hhpt, ncol=1)

##�p�����[�^�̐ݒ�
beta1 <- -5.8   #���i�̃p�����[�^
beta2 <- 5.5   #�������̃p�����[�^
beta3 <- 2.0   #���ʒ�̃p�����[�^
beta4 <- 1.8   #�L�����y�[���̃p�����[�^
b1 <- c(1.1, 0.6, -0.7, -0.3)   #�J�e�S���[���C�����e�B�̃p�����[�^
b0 <- c(0.5, 0.9, 1.4, 2.1)   #�u�����h1�`4�̑��΃x�[�X�̔���
betat <- c(b0, beta1, beta2, beta3, beta4)

##�ϗʌ��ʂ̐ݒ�
k.random <- 5   #�ϗʌ��ʂ̕ϐ���
b0.random <- matrix(b0, nrow=hh, ncol=choise-1, byrow=T) + 
  matrix(rnorm(hh*(choise-1), 0, 1), nrow=hh, ncol=choise-1, byrow=T)
beta1.random <- matrix(beta1, nrow=hh, ncol=1, byrow=T) + 
  matrix(rnorm(hh, 0, sqrt(3.5)), nrow=hh, ncol=1, byrow=T)


##���p�𔭐������A�I�����ꂽ�u�����h������
#���W�b�g�̔���
logit <- matrix(0, nrow=hhpt, ncol=st)
for(i in 1:hh){
  r <- subset(1:hhpt, ID[, 2]==i)
  for(j in 1:(st-1)){
    logit[r, j] <- b0.random[i, j] + beta1.random[i]*PRICE[r, j]  + beta2*DISC[r, j] + beta3*DISP[r, j] + beta4*CAMP[r, j] 
  }
}

#��ϐ��̃��W�b�g���v�Z
for(i in 1:hh){
  r <- subset(1:hhpt, ID[, 2]==i)
  logit[r, st] <- beta1.random[i]*PRICE[r, st] + beta2*DISC[r, st] + beta3*DISP[r, st] + beta4*CAMP[r, st] 
}

##�������������W�b�g����I���u�����h������
#�u�����h�I���m�����v�Z
Pr <- exp(logit)/rowSums(exp(logit))
colMeans(Pr); apply(Pr, 2, summary)

#�I���u�����h�𔭐�
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
colMeans(Y); apply(Y, 2, table)

round(cbind(Y %*% 1:choise, Pr), 3)   #�I�����ʂƑI���m��

####�}���R�t�A�������e�J�����@�ŕϗʌ��ʃ��W�X�e�B�b�N��A���f���𐄒�####
##��A���f���𐄒肷�邽�߂ɐ����ϐ����x�N�g���`���ɕύX�ݒ�
#id��ݒ�
id.v <- c()
for(i in 1:hh){
  id.v <- c(id.v, rep(ID[ID[, 2]==i, 2], choise))
}

#�ؕЂ̐ݒ�
p <- c(1, rep(0, choise))
bp <- matrix(p, nrow=hhpt*(choise+1), ncol=choise, byrow=T)
BP <- subset(bp, rowSums(bp) > 0)
BP <- BP[, -st] 

#�J�e�S�����C�����e�B�̐ݒ�
index.royl <- rep(1:hhpt, rep(choise, hhpt))
ROYL.v <- matrix(0, nrow=hhpt*choise, ncol=choise)

for(i in 1:hhpt){
  ROYL.v[index.royl==i, ] <- diag(c(rep(ROYL[i, ], choise-1), 0))
}
ROYL.v <- ROYL.v[, -st]

#�����ϐ��̐ݒ�
PRICE.v <- as.numeric(t(PRICE))
DISC.v <- as.numeric(t(DISC))
DISP.v <- as.numeric(t(DISP))
CAMP.v <- as.numeric(t(CAMP))

round(X <- data.frame(b=BP, PRICE=PRICE.v, DISC=DISC.v, DISP=DISP.v, CAMP=CAMP.v), 2)   #�f�[�^�̌���
XM <- as.matrix(X)


#Z�̐ݒ�
Z <- matrix(0, nrow=nrow(XM), ncol=k.random*hh)
for(i in 1:hh){
  r <- subset(1:nrow(XM), id.v==i)
  c <- ((i-1)*k.random+1):((i-1)*k.random+k.random)
  Z[r, c] <- XM[id.v==i, 1:k.random]
}


##MCMC�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- c()   #�ΐ��ޓx�̕ۑ��p

##���O���z�̐ݒ�
#�Œ���ʂ̎��O���z
betas.fix <- rep(0, ncol(XM))  #��A�W���̕��ς̎��O���z
sigma.fix <- diag(rep(0.01, ncol(XM)))   #��A�W���̎��O���z�̕��U

#�ϗʌ��ʂ̎��O���z
Deltabar <- rep(0, hh*choise)
Adelta <- 0.01*diag(k.random)
nu <- k.random   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(rep(1, k.random))
beta.random <- matrix(0, nrow=hh, ncol=k.random)   #�ϗʌ��ʂ̎��O���z�̕��ς�0�ɌŒ�


##�T���v�����O���ʂ̕ۑ��p
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
B.random <- array(0, dim=c(hh, k.random, R/keep))
SIGMA <- matrix(0, nrow=R/keep, ncol=k.random^2)

##�����l�̐ݒ�
#��A�W���̏����l
oldbeta.f <- c(c(runif(choise-1, 0, 3)), -5.0, 5.0, runif(2, 1, 2))  

#�ϗʌ��ʂ̏����l
cov.random <- diag(runif(k.random, 0.5, 2))
oldbeta.r <- mvrnorm(hh, c(rep(0, k.random-1), -5), cov.random)   #�ϗʌ��ʂ̏����l


#�K�w���f���̏����l
beta.random <- matrix(0, nrow=hh, ncol=k.random)

####�}���R�t�A�������e�J�����@�Ő���####
##mixed logit���f���̑ΐ��ޓx
LLike <- function(beta, b, X, Z, Y, hh, choise){
  d <- rowSums(exp(matrix(X %*% beta + Z %*% b, nrow=hh, ncol=choise, byrow=T)))
  LLl <- rowSums(Y * matrix(X %*% beta + Z %*% b, nrow=hh, ncol=choise, byrow=T)) - log(d)
  LL <- sum(LLl)
  LL.val<- list(LLl=LLl, LL=LL)
  return(LL.val)
}

##�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O
foreach(rp = 1:R) %do% {
  ##MH�T���v�����O�ŌŒ����beta�̃T���v�����O
  oldbeta.rv <- as.numeric(t(oldbeta.r))
  betad.f <- oldbeta.f
  betan.f <- betad.f + rnorm(length(betad.f), 0, 0.025)
  #�����_���E�H�[�N�T���v�����O
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew.f <- LLike(beta=betan.f, b=oldbeta.rv, X=XM, Z=Z, Y=Y, hh=hhpt, choise=choise)$LL
  logold.f <- LLike(beta=betad.f, b=oldbeta.rv, X=XM, Z=Z, Y=Y, hh=hhpt, choise=choise)$LL
  logpnew.f <- lndMvn(betan.f, betas.fix, sigma.fix)
  logpold.f <- lndMvn(betad.f, betas.fix, sigma.fix)
  
  #MH�T���v�����O
  alpha.f <- min(1, exp(lognew.f + logpnew.f - logold.f - logpold.f))
  if(alpha.f == "NAN") alpha.f <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V�����Œ����beta���̑�
  if(u < alpha.f){
    oldbeta.f <- betan.f
    logl.f <- lognew.f
    
    #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
  } else {
    logl.f <- logold.f
  }
  
  
  ##MH�T���v�����O�Ōl�ʂɕϗʌ���beta���T���v�����O
  betad.random <- oldbeta.r 
  rw <- t(0.5 * chol(cov.random) %*% t(matrix(rnorm(hh*(k.random)), nrow=hh, ncol=k.random)))
  
  betan.random <- betad.random + rw
  betad.r <- as.numeric(t(betad.random))
  betan.r <- as.numeric(t(betan.random))
  inv.cov <- solve(cov.random)
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew.r <- LLike(beta=oldbeta.f, b=betan.r, X=XM, Z=Z, Y=Y, hh=hhpt, choise=choise)$LLl
  logold.r <- LLike(beta=oldbeta.f, b=betad.r, X=XM, Z=Z, Y=Y, hh=hhpt, choise=choise)$LLl
  logpnew.r <- apply((betan.random - beta.random), 1, function(x) -0.5 * x %*% inv.cov %*% x)
  logpold.r <- apply((betad.random - beta.random), 1, function(x) -0.5 * x %*% inv.cov %*% x)
  
  #ID�ʂɑΐ��ޓx�̘a�����
  lognew.rind <- data.frame(logl=lognew.r, id=ID[, 2]) %>%
                  dplyr::group_by(id) %>%
                  dplyr::summarize(sum=sum(logl))
  
  logold.rind <- data.frame(logl=logold.r, id=ID[, 2]) %>%
                  dplyr::group_by(id) %>%
                  dplyr::summarize(sum=sum(logl))
  
  #MH�T���v�����O
  rand <- matrix(runif(hh), nrow=hh, ncol=k.random)
  LLind.diff <- exp(lognew.rind$sum + logpnew.r - logold.rind$sum - logpold.r)   #���p�����v�Z
  alpha <- matrix(ifelse(LLind.diff > 1, 1, LLind.diff), nrow=hh, ncol=k.random)      
  
  oldbeta.r <- ifelse(alpha > rand, oldbeta.r <- betan.random, oldbeta.r <- betad.random)   #alpha��rand�������Ă�����̑�
  logl <- ifelse(alpha[, 1] > rand[, 1], logl <- lognew.r, logl <- logold.r)
  
  
  ##�t�E�B�V���[�g���z����sigma���T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^
  V <- var(oldbeta.r)
  VK <- k.random * diag(k.random) + hh * V
  nu1 <- hh + nu - 1 
  
  #�t�E�B�V���[�g���z���番�U�����U�s��𔭐�
  cov.random <- rwishart(nu1, solve(VK))$IW   
  

  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta.f
    B.random[, , mkeep] <- oldbeta.r
    SIGMA[mkeep, ] <- as.numeric(cov.random)
    
    print(sum(logl))
    print(rp)
    print(round(mean(alpha), 3)); print(round(alpha.f, 3))
    print(round(rbind(oldbeta.f, betat), 3))
    print(round(cov.random, 3))
  }
}

####���茋�ʂƗv��####
burnin <- 2500
i <- 6

matplot(BETA[, 1:4], type="l")
matplot(BETA[, 5:8], type="l")
matplot(SIGMA[, c(1, 7, 13)], type="l")
matplot(SIGMA[, c(19, 25)], type="l")
matplot(t(B.random[i, 1:2, ]), type="l")
matplot(t(B.random[i, 3:4, ]), type="l")
plot(1:5000, t(B.random[i, 5, ]), type="l")

#�ϗʌ��ʂ̕��U�̎��㕽��
colMeans(SIGMA[burnin:nrow(SIGMA), c(1, 7, 13, 19, 25)]) 
colMeans(t(B.random[i, , burnin:nrow(SIGMA)])) 
c(b0.random[i, ]-b0, beta1.random[i, ]-beta1)