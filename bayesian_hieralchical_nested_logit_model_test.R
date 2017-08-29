#####�K�w�x�C�Y�l�X�e�b�h���W�b�g���f��#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(extraDistr)
library(gtools)
library(mlogit)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####�C�ӂ̕��U�����U�s����쐬������֐�####
##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
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
hh <- 2000
member <- 9   #�����o�[��
c_num <- 8   #�ߑ��p�^�[��
hhpt <- hh*member*c_num   #�S�ϐ���


####�����ϐ��̔���####
##�̓������ϐ��̔���
#�����o�[�̐����ϐ��̐ݒ�
Mus <- matrix(as.numeric(table(1:hhpt, rep(rep(1:member, rep(c_num, member)), hh))), nrow=hhpt, ncol=member)
colnames(Mus) <- c("hono", "koto", "umi", "rin", "hana", "maki", "nico", "eri", "nozo")
Mus <- Mus[, -ncol(Mus)]

#�ߑ��̐����ϐ��̐ݒ�
ct <- c(1, rep(0, c_num))
cloth <- matrix(ct, nrow=hh*member*(c_num+1), ncol=c_num, byrow=T)
CLOTH <- subset(cloth, rowSums(cloth) > 0)[, -c_num]
colnames(CLOTH) <- c("A", "B", "C", "D", "F", "G", "H")

#�J�[�h���
type <- 3   #��ސ�
card <- t(rmultinom(hhpt, 1, c(3/member, 3/member, 3/member)))
colnames(card) <- c("UR", "SSR", "SR")
CARD <- card[, -type]

#�v�����[�V�����ڐG��
Prom <- scale(rpois(hhpt, 5))

#�f�[�^�̌���
X <- data.frame(Mus, CLOTH, CARD, Prom)
XM1 <- as.matrix(X)
XM2 <- XM1[, -ncol(XM1)]

#�p�����[�^��
k1 <- 2 + ncol(XM1) + ncol(XM2) - (member-1)
k2 <- 2
k11 <- 1:(ncol(XM1)+1)
k12 <- c(ncol(XM1)+2, 2:member, (ncol(XM1)+3):k1)
k13 <- k12[1:(length(k12)-2)]

##ID�̐ݒ�
id <- rep(1:hh, rep(member*c_num, hh))
pt <- rep(1:(member*c_num), hh)
ID <- data.frame(no=1:hhpt, id=id, pt=pt)


##�̊Ԑ����ϐ��̔���
#�A���ϐ�
cont.h <- 3
Z.cont <- matrix(runif(hh*cont.h, 0, 1), nrow=hh, ncol=cont.h) 

#��l�ϐ�
bin.h <- 3
Z.bin <- matrix(0, nrow=hh, ncol=bin.h)
for(i in 1:bin.h){
  p.bin <- runif(1, 0.2, 0.8)
  Z.bin[, i] <- rbinom(hh, 1, p.bin)  
}

#���l�ϐ�
multi.h <- 3
p.multi <- runif(multi.h)
Z.multi <- t(rmultinom(hh, 1, p.multi))
freq.min <- which.min(colSums(Z.multi))
Z.multi <- Z.multi[, -freq.min]

#�f�[�^�̌���
Z <- data.frame(cont=Z.cont, bin=Z.bin, multi=Z.multi)
ZM <- as.matrix(Z)
k3 <- ncol(ZM)

####�����ϐ��̔���####
##�̊ԃp�����[�^�̐ݒ�
#�̊ԕ��U�����U�s���ݒ�
Cov <- corrM(col=k1, lower=-0.55, upper=0.9, eigen_lower=0.025, eigen_upper=0.35)   

#�̊ԉ�A�p�����[�^�̐ݒ�
theta1 <- c(runif(1, -1.1, -0.55), runif(member-1, -0.2, 0.85), runif(c_num-1, -0.5, 0.7), 
            runif(1, -0.7, -0.5), runif(1, -0.5, -0.3), runif(1, 0.05, 0.15), runif(1, -1.3, -1.0), 
            runif(c_num-1, -0.5, 0.7), runif(1, -0.8, -0.5), runif(1, -0.55, -0.3))

theta2 <- matrix(c(runif(k3, -0.4, 0.5), runif((member-1)*k3, -0.4, 0.4), runif((c_num-1)*k3, -0.5, 0.5), 
                   runif(k3, -0.4, 0.4), runif(k3, -0.4, 0.4), runif(k3, -0.1, 0.15), runif(k3, -0.4, 0.4), 
                   runif((c_num-1)*k3, -0.45, 0.6), runif(k3, -0.4, 0.3), runif(k3, -0.3, 0.3)), nrow=k3, ncol=k1, byrow=T)

#�p�����[�^�̌���
theta <- rbind(theta1, theta2)
rownames(theta) <- 1:nrow(theta) 


##�̓���A�p�����[�^�̐ݒ�
#�l�ʂ̉�A�p�����[�^
beta <- cbind(1, ZM) %*% theta + mvrnorm(hh, rep(0, k1), Cov/5)   #��A�W���̃p�����[�^

#���O�T���ϐ��̃p�����[�^
rho.par <- c(1.0, 0)
rho <- exp(rho.par)/(1+exp(rho.par))

#��A�W���̐ݒ�
beta1 <- beta[, k11]
beta2 <- beta[, k12]
beta3 <- beta[, k13]

##�����ϐ��̔���
#���W�b�g�ƃ��O�T���ϐ����`
logit1.list <- list()
logit2.list <- list()
logsum.list <- list()

#�S���[�U�[�̌l���Ƃ̃��O�T���ϐ��ƃ��W�b�g���v�Z
for(i in 1:hh){
  print(i)
  #���O�T���ϐ��̒�`
  logsum.list[[i]] <- log(1 + exp(cbind(1, XM2[ID$id==i, 1:(ncol(XM2)-2)]) %*% beta3[i, ]))
  logsum.ind <- matrix(logsum.list[[i]], nrow=member*c_num, ncol=k2)*CARD[ID$id==i, ]
  
  #���W�b�g�̒�`
  logit1.list[[i]] <- cbind(1, XM1[ID$id==i, ]) %*% beta1[i, ] + logsum.ind %*% rho
  logit2.list[[i]] <- cbind(1, XM2)[ID$id==i, ] %*% beta2[i, ]
}

#���X�g�𐔒l�^�ɕύX
logsum <- unlist(logsum.list)
logit1 <- unlist(logit1.list)
logit2 <- unlist(logit2.list)

#�l�X�e�b�h���W�b�g���f���ɂɂ��m�����v�Z���A�x���k�[�C���z��艞���ϐ��𔭐�
#�J�[�h�������Ă��邩�ǂ���
Pr1 <- exp(logit1)/(1+exp(logit1))
y1 <- rbinom(hhpt, 1, Pr1)

#�J�[�h���o�����Ă��邩�ǂ���
Pr2 <- exp(logit2)/(1+exp(logit2))
y2 <- rbinom(hhpt, 1, Pr2)
y2[y1==0] <- NA   #�J�[�h�������Ă���ꍇ�̂݊o���L�����`����


####�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y�l�X�e�b�h���W�b�g���f���𐄒�####
####�}���R�t�A�������e�J�����@�̐���̂��߂̏���####
##�l�X�e�b�h���W�b�g���f���̑ΐ��ޓx�֐�
loglike <- function(x, y1, y2, XM1, XM2, CARD, type, member, c_num, hhpt, k11, k12, k13, k2){
  
  #�p�����[�^�̐ݒ�
  beta1 <- x[k11]
  beta2 <- x[k12]
  beta3 <- x[k13]
  rho <- exp(x[(length(x)-1):length(x)])/(1+exp(x[(length(x)-1):length(x)]))   #�p�����[�^��0�`1
  
  #���O�T���ϐ��̒�`
  logsum <- log(1 + exp(cbind(1, XM2[, 1:(ncol(XM2)-2)]) %*% beta3))  #���O�T���ϐ�
  logsum.card <- matrix(logsum, nrow=hhpt, ncol=k2)*CARD
  
  #���W�b�g�̒�`
  logit1 <- cbind(1, XM1) %*% beta1 + logsum.card %*% rho    #�J�[�h���L�L���̃��W�b�g
  logit2 <- cbind(1, XM2) %*% beta2   #�o���L���̃��W�b�g
  
  
  #�ΐ��ޓx���`����
  #�J�[�h���L�L���̑ΐ��ޓx
  Pr.l <- exp(logit1) / (1 + exp(logit1))
  LLs.l <- y1*log(Pr.l) + (1-y1)*log(1-Pr.l)  
  LL.l <- sum(LLs.l)
  
  #�o���L���̑ΐ��ޓx
  Pr.b <- exp(logit2[y1==1]) / (1 + exp(logit2[y1==1]))
  LLs.b <- y2[y1==1]*log(Pr.b) + (1-y2[y1==1])*log(1-Pr.b)  
  LL.b <- sum(LLs.b)
  
  #�ΐ��ޓx�����v
  LL <- LL.l + LL.b
  return(LL)
}

##���j���[�g���@�őΐ��ޓx���ő剻����
#���肳�ꂽ�p�����[�^�������l�̎Q�l�ɂ���@
x <- c(runif(k1, -0.5, 0.5), 0.8, 0.5)
res <- optim(x, loglike, y1=y1, y2=y2, XM1=XM1, XM2=XM2, CARD=CARD, type=type, member=member, c_num=c_num, hhpt=hhpt,
             k11=k11, k12=k12, k13=k13, k2=k2, method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))

#���茋��
x1 <- res$par
H <- sqrt(-diag(solve(res$hessian)))


##MCMC�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4
ZMi <- cbind(1, ZM)

#�f�[�^�̐ݒ�
n <- member*c_num
XM11 <- cbind(1, XM1)
XM12 <- cbind(1, XM2)


##�C���f�b�N�X���쐬
#�p�����[�^�̃C���f�b�N�X���쐬
index.list <- list()
for(i in 1:hh){
  index.list[[i]] <- ID$no[ID$id==i]
}
user.count <- as.numeric(table(ID$id))
index.ones <- subset(1:hhpt, y1==1) 

index.par1 <- k11
index.par2 <- k12
index.logsum <- k12[1:(length(k12)-2)]


#�ΐ��ޓx�̕ۑ��p�z��
logpnew1 <- matrix(0, nrow=hh, ncol=1)
logpold1 <- matrix(0, nrow=hh, ncol=1)


##���O���z�̐ݒ�
Deltabar <- matrix(0, nrow=ncol(cbind(1, ZM)), ncol=k1)   #��A�p�����[�^�̊K�w���f���̉�A�W���̕��ς̎��O���z
Adelta <- 0.01 * diag(rep(1, ncol(ZMi)))   #��A�p�����[�^�̊K�w���f���̉�A�W���̕��U�̎��O���z
nu <- (ncol(ZM)+1)+k1   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(k1)   #�t�E�B�V���[�g���z�̃p�����[�^
mu <- rep(-3, 2)
sigma <- 0.01 * diag(2)

##�T���v�����O���ʂ̕ۑ��p�z��
#�p�����[�^�̕ۑ��p�z��
BETA <- array(0, dim=c(hh, k1, R/keep))
RHO <- matrix(0, nrow=R/keep, ncol=k2)
THETA <- matrix(0, nrow=R/keep, ncol=ncol(cbind(1, ZM))*k1)
SIGMA <- matrix(0, nrow=R/keep, ncol=k1*k1)

#���p���Ƒΐ��ޓx�̕ۑ��p�z��
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##�����l�̐ݒ�
#��A�y�����[�^�̏����l
tau1 <- mvrnorm(hh, rep(0, k1), diag(0.3, k1))
oldbetas <- matrix(x1[1:k1], nrow=hh, ncol=k1, byrow=T) + tau1

oldBetas <- matrix(0, nrow=hhpt, ncol=k1) 
for(i in 1:hh){
  oldBetas[index.list[[i]], ] <- matrix(oldbetas[i, ], nrow=user.count[i], ncol=k1, byrow=T)
}

#���O�T���ϐ��̏����l
oldrho <- c(0, 0)

#���O���z�̃p�����[�^�̏����l
oldDelta <- ginv(t(ZMi) %*% ZMi) %*% t(ZMi) %*% oldbetas
oldVbeta <- 1/hh * (t(oldbetas - ZMi %*% oldDelta) %*% (oldbetas - ZMi %*% oldDelta))
oldVbeta_inv <- solve(oldVbeta)


##�l�X�e�b�h���W�b�g���f���̑ΐ��ޓx�֐�
LL_nest <- function(Betas, rho, ID, y1, y2, XM1, XM2, CARD, hhpt, index.par1, index.par2, index.logsum, index.ones, k2, z){
  
  #�p�����[�^�̐ݒ�
  Beta1 <- Betas[, index.par1]
  Beta2 <- Betas[, index.par2]
  rho1 <- exp(rho)/(1+exp(rho))   #�p�����[�^��0�`1
  
  #���O�T���ϐ��̒�`
  logsum <- log(1 + exp(rowSums(XM2[, 1:(ncol(XM2)-2)] * Betas[, index.logsum])))  #���O�T���ϐ�
  logsum.mx <- matrix(logsum, nrow=hhpt, ncol=k2) * CARD
  
  #���W�b�g�̒�`
  logit1 <- rowSums(XM1 * Beta1) + logsum.mx %*% rho1    #�J�[�h���L�L���̃��W�b�g
  logit2 <- rowSums(XM2 * Beta2)   #�o���L���̃��W�b�g
  
  #�ΐ��ޓx���`����
  #�J�[�h���L�L���̑ΐ��ޓx
  Pr.l <- exp(logit1) / (1 + exp(logit1))
  LL.l <- y1*log(Pr.l) + (1-y1)*log(1-Pr.l)  
  
  #�o���L���̑ΐ��ޓx
  LL.b <- matrix(0, nrow=hhpt, ncol=1)
  Pr.b <- exp(logit2[index.ones]) / (1 + exp(logit2[index.ones]))
  LLl.b <- y2[index.ones]*log(Pr.b) + (1-y2[index.ones])*log(1-Pr.b)  
  LL.b[index.ones, ] <- LLl.b
  
  #�p�^�[���ɂ���đΐ��ޓx�̘a�̎�����ς���
  if(z==1){
    #�Œ���ʂ̑ΐ��ޓx�̘a
    LL <- sum(LL.l) + sum(LL.b)
  } else {
    #�ϗʌ��ʂ̑ΐ��ޓx�̘a
    LL <- as.matrix(data.frame(logl=LL.l+LL.b, id=ID$id) %>%
                      dplyr::group_by(id) %>%
                      dplyr::summarize(sum=sum(logl)))[, 2]
  }
  return(LL)
}


#�^�̃p�����[�^�ł̑ΐ��ޓx
Betat <- matrix(0, nrow=hhpt, ncol=k1)
for(i in 1:hh){
  Betat[index.list[[i]], ] <- matrix(beta[i, ], nrow=user.count[i], ncol=k1, byrow=T)
}

#�ΐ��ޓx�̘a
LLind <- LL_nest(Betas=Betat, rho=rho, ID=ID, y1=y1, y2=y2, XM1=XM11, XM2=XM12, CARD=CARD, hhpt=hhpt, 
                 index.par1=index.par1, index.par2=index.par2, index.logsum=index.logsum, index.ones=index.ones, k2=k2, z=2)
LL_t <- sum(LLind)


####�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y�l�X�e�b�h���W�b�g���f���̃p�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##���g���|���X-�w�C�X�e�B���O�@�Ōl�ʂ�beta���T���v�����O
  #�����_���E�H�[�N�T���v�����O
  rw1 <- mvrnorm(hh, rep(0, k1), diag(H[1:k1]))
  betad <- oldbetas.r <- oldbetas
  betan <- oldbetas + rw1
  
  #�p�����[�^���f�[�^�s��ƍ��킹��
  Betad <- oldBetas
  Betan <- matrix(0, nrow=hhpt, ncol=k1) 
  
  #�p�����[�^�̎��O���z�ƌ덷���v�Z
  er_new <- betan - ZMi %*% oldDelta
  er_old <- betad - ZMi %*% oldDelta
  
  ##ID�ʂɑΐ��ޓx�Ƒΐ����O���z���v�Z
  #�ΐ����O���z���v�Z
  for(i in 1:hh){
    Betan[index.list[[i]], ] <- matrix(betan[i, ], nrow=user.count[i], ncol=k1, byrow=T)   #�p�����[�^���f�[�^�s��
    logpnew1[i, ] <- -0.5 * (er_new[i, ] %*% oldVbeta_inv %*% er_new[i, ])
    logpold1[i, ] <- -0.5 * (er_old[i, ] %*% oldVbeta_inv %*% er_old[i, ])
  }
  
  #�ΐ��ޓx���v�Z
  lognew1 <- LL_nest(Betas=Betan, rho=rho, ID=ID, y1=y1, y2=y2, XM1=XM11, XM2=XM12, CARD=CARD, hhpt=hhpt, 
                     index.par1=index.par1, index.par2=index.par2, index.logsum=index.logsum, index.ones=index.ones, 
                     k2=k2, z=2)
  logold1 <- LL_nest(Betas=Betad, rho=rho, ID=ID, y1=y1, y2=y2, XM1=XM11, XM2=XM12, CARD=CARD, hhpt=hhpt, 
                     index.par1=index.par1, index.par2=index.par2, index.logsum=index.logsum,  index.ones=index.ones, 
                     k2=k2, z=2)
  
  ##�T���v�����O���̑����邩�ǂ���������
  #��l�������痐���𔭐�
  rand1 <- runif(hh)  
  
  #�̑𗦂��v�Z
  LLind.diff1 <- exp(lognew1 + logpnew1 - logold1 - logpold1)
  alpha1 <- ifelse(LLind.diff1 > 1, 1, LLind.diff1)
  
  #alpha�Ɋ�Â�beta��rho���̑����邩�ǂ�������
  index.adopt1 <- subset(1:hh, alpha1 > rand1)
  oldbetas.r[index.adopt1, ] <- betan[index.adopt1, ]
  
  for(i in 1:hh){
    oldBetas[index.list[[i]], ] <- matrix(oldbetas.r[i, ], nrow=user.count[i], ncol=k1, byrow=T)   
  }
  
  #�̑𗦂Ƒΐ��ޓx���v�Z
  adopt <- sum(oldbetas[, 1] != oldbetas.r[, 1])/hh
  LLho <- sum(lognew1[index.adopt1]) + sum(logold1[-index.adopt1])   #�ΐ��ޓx
  
  #�p�����[�^���X�V
  oldbetas <- oldbetas.r
  
  
  ##���g���|���X-�w�C�X�e�B���O�@��rho���T���v�����O
  #�����_���E�H�[�N�T���v�����O
  rw2 <- mvrnorm(1, rep(0, k2), diag(0.01, k2))
  rhod <- oldrho.r <- oldrho
  rhon <- oldrho + rw2
  
  ##�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew2 <- LL_nest(Betas=oldBetas, rho=rhon, ID=ID, y1=y1, y2=y2, XM1=XM11, XM2=XM12, CARD=CARD, hhpt=hhpt, 
                     index.par1=index.par1, index.par2=index.par2, index.logsum=index.logsum, index.ones=index.ones, 
                     k2=k2, z=1)
  logold2 <- LL_nest(Betas=oldBetas, rho=rhod, ID=ID, y1=y1, y2=y2, XM1=XM11, XM2=XM12, CARD=CARD, hhpt=hhpt, 
                     index.par1=index.par1, index.par2=index.par2, index.logsum=index.logsum, index.ones=index.ones, 
                     k2=k2, z=1)
  logpnew2 <- lndMvn(rhon, mu, sigma)
  logpold2 <- lndMvn(rhod, mu, sigma)
  
  ##�T���v�����O���̑����邩�ǂ���������
  alpha2 <- min(1, exp(lognew2 + logpnew2 - logold2 - logpold2))
  if(alpha2 == "NAN") alpha2 <- -1
  
  #��l�����𔭐�
  rand2 <- runif(1)
  
  #rand2 < alpha�Ȃ�V����rho���̑�
  if(rand2 < alpha2){
    oldrho.r <- rhon
    
    #�����łȂ��Ȃ�rho���X�V���Ȃ�
  } else {
    oldrho.r <- rhod
  }
  
  #�p�����[�^���X�V
  oldrho <- oldrho.r   
  
  
  ##���ϗʉ�A���f���ɂ��K�w���f���̃T���v�����O
  #��A�p�����[�^�̑��ϗʉ�A
  out <- rmultireg(Y=oldbetas, X=ZMi, Bbar=Deltabar, A=Adelta, nu=nu, V=V)
  oldDelta <- out$B
  oldVbeta <- out$Sigma
  oldVbeta_inv <- solve(oldVbeta)
  
  
  ##�T���v�����O���ʂ�ۑ�
  if(rp%%keep==0){
    cat("��*'��')�� <�������[ ������Ƃ܂��Ă�", paste(rp/R*100, "��"), fill=TRUE)
    mkeep <- rp/keep
    BETA[, , mkeep] <- oldbetas
    RHO[mkeep, ] <- oldrho
    THETA[mkeep, ] <- as.numeric(oldDelta)
    SIGMA[mkeep, ] <- as.numeric(oldVbeta)
    print(round(c(adopt, LLho, LL_t, res$value), 2))
    print(round(rbind(oldbetas[1:3, 1:22], beta[1:3, 1:22]), 2))
    print(round(rbind(colMeans(oldbetas), colMeans(beta)), 2))
    print(round(c(exp(oldrho)/(1+exp(oldrho)), rho), 3))
    print(round(cbind(oldDelta[, c(k11[1:5], k12[1:5])], theta[, c(k11[1:5], k12[1:5])]), 2))
    print(rp)
  }
}