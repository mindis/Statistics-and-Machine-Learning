#####�x�C�W�A���K�w�L�����������v���r�b�g���f��#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(gtools)
library(MNP)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
library(lattice)

####���ϗʐ��K�����𔭐�������֐�####
##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
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
#set.seed(8437)
##�f�[�^�̐ݒ�
hh <- 500   #�v���C���[��
pt <- rpois(hh, 10)   #�I���@�
pt <- ifelse(pt==0, 1, pt)   #�I���@���0�Ȃ�1�ɒu������
hhpt <- sum(pt)   #���T���v����
member <- 10   #�I���\�����o�[��
st <- 10   #������o�[
g <- 3   #�Z�O�����g��

##ID�̐ݒ�
id <- rep(1:hh, pt)
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}

ID <- data.frame(no=1:hhpt, id=id, t=t)   #ID�f�[�^�̌���
id_r <- matrix(1:(hhpt*(member-1)), nrow=hhpt, ncol=member-1, byrow=T)


####�����ϐ��̔���####
##���U�I�����f���̐����ϐ��̔���
#�����t���̐����ϐ��̔���
X1.cont <- matrix(rnorm(hhpt*member*2, 0, 1), nrow=hhpt, ncol=member*2)

X1.bin <- matrix(0, nrow=hhpt, ncol=member*2)
for(i in 1:(member*2)){
  X1.bin[, i]  <- rbinom(hhpt, 1, runif(1, 0.35, 0.6))
}

#������o�[�Ƃ̑��ΐ����ϐ�
X1.cont_r <- cbind(X1.cont[, 1:(member-1)] - X1.cont[, member], X1.cont[, (member+1):(2*member-1)] - X1.cont[, 2*member])
X1.bin_r <- cbind(X1.bin[, 1:(member-1)] - X1.bin[, member], X1.bin[, (member+1):(2*member-1)] - X1.bin[, 2*member])

#�����^�̐����ϐ��̔���
X2.cont <- c()
X2.bin <- c()

for(i in 1:hh){
  bin <- rbinom(1, 1, runif(1, 0.35, 0.7))
  X2.cont <- c(X2.cont, rep(rnorm(1, 0, 1), pt[i]))
  X2.bin <- c(X2.bin, rep(bin, pt[i]))
}

##�x�N�g���^�����ϐ��Ƀf�[�^�t�H�[�}�b�g��ύX
#�ؕЂ̐ݒ�
p <- c(1, rep(0, member-1))
pop <- matrix(p, nrow=hhpt*member, ncol=member-1, byrow=T)
POP <- pop[rowSums(pop) > 0, ]

#�����t�������ϐ��̐ݒ�
X1.cont_v <- cbind(as.numeric(t(X1.cont_r[, 1:(member-1)])), as.numeric(t(X1.cont_r[, member:(2*(member-1))])))
X1.bin_v <- cbind(as.numeric(t(X1.bin_r[, 1:(member-1)])), as.numeric(t(X1.bin_r[, member:(2*(member-1))])))

#�����^�����ϐ��̐ݒ�
X2.v <- matrix(0, nrow=hhpt*(member-1), ncol=(member-1)*2)
for(i in 1:hhpt){
  r <- ((i-1)*(member-1)+1):((i-1)*(member-1)+member-1)
  X2.v[r, ] <- cbind(diag(X2.cont[i], member-1), diag(X2.bin[i], member-1))
}

#�f�[�^�̌���
X <- data.frame(mu=POP, c=X1.cont_v, b=X1.bin_v, m=X2.v[, 1:(member-1)])
XM <- as.matrix(X)
round(XM, 2)


##�K�w���f���̐����ϐ��̔���
#�A���ϐ��̔���
cont <- 3
Z.cont <- matrix(rnorm(hh*cont, 0, 1), nrow=hh, ncol=cont)

#��l�ϐ��̔���
bin <- 3
Z.bin <- matrix(0, nrow=hh, ncol=bin)
for(i in 1:bin){
  Z.bin[, i] <- rbinom(hh, 1, runif(1, 0.4, 0.6))
}

#���l�ϐ��̔���
multi <- 4
p <- runif(multi, 0.25, 2)
Z.multi <- t(rmultinom(hh, 1, p))
Z.multi <- Z.multi[, -which.min(colSums(Z.multi))]

#�f�[�^�̌���
Zx <- data.frame(c=Z.cont, b=Z.bin)


##�x�N�g���^�����ϐ��Ƀf�[�^�t�H�[�}�b�g��ύX
#�ؕЂ̐ݒ�
p <- c(1, rep(0, g))
int <- matrix(p, nrow=hh*(g+1), ncol=g, byrow=T)
INT <- int[rowSums(int) > 0, -g]

#�����ϐ����x�N�g���^�ɕύX
Zi.v <- matrix(0, nrow=nrow(Zx)*g, ncol=(cont+bin)*2)

for(i in 1:hh){
  index <- ((i-1)*g+1):((i-1)*g+g)
  
  diag.x <- c()
  for(j in 1:(cont+bin)){
    diag.x <- cbind(diag.x, diag(Zx[i, j], g)[, -g])
  }
  
  #  diag.m <- matrix(0, nrow=g, ncol=(multi-1)*2)
  #  for(j in 1:(g-1)){
  #    r <- ((j-1)*g+1):((j-1)*g+g)
  #    diag.m[j, r] <- as.numeric(Zx[i, (cont+bin+1):ncol(Zx)])
  # }
  Zi.v[index, ] <- cbind(diag.x)
}

#�f�[�^�̌���
Zx.v <- cbind(INT, Zi.v)
round(Zx.v, 3)   #�f�[�^�̊m�F


####�����ϐ��̔���####
##�������W�b�g���f�����Z�O�����g�����𔭐�
#�p�����[�^�̐ݒ�
for(i in 1:1000){
  theta.z <- c(runif(g-1, -0.75, 0.75), runif(cont*(g-1), 0, 1), runif(bin*(g-1), -1.2, 1.2))

  #���W�b�g�Ɗm���̌v�Z
  logit <- matrix(Zx.v %*% theta.z, nrow=hh, ncol=g, byrow=T)   #���W�b�g�̐ݒ�
  Pr.z <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hh, ncol=g)   #�m���̌v�Z
  
  #�������z���Z�O�����g�𐶐�
  Z <- t(apply(Pr.z, 1, function(x) rmultinom(1, 1, x)))
  if(min(colMeans(Z)) > 0.25) {break}
}
Zi <- Z %*% 1:g
r_rate <- colSums(Z)/sum(Z)   #������
colSums(Z); colMeans(Z)


##ID���ƂɃZ�O�����g������
zi <- c()
z <- c()

for(i in 1:hh){
  zi <- c(zi, rep(Zi[i], pt[i]))
  z <- rbind(z, matrix(Z[i, ], nrow=pt[i], ncol=g, byrow=T))
}

ID <- data.frame(ID, z=zi)   #ID�Ɍ���

##ID���x�N�g���`���ɕύX
#�Z�O�����g���蓖�Ă�ύX
zv1 <- c(); zv2 <- c(); idv <- c(); tv <- c()

for(i in 1:(member-1)){
  zv1 <- cbind(zv1, z)
  zv2 <- cbind(zv2, ID$z)
  idv <- cbind(idv, ID$id)
  tv <- cbind(tv, ID$t)
}

#�x�N�g���`���̃f�[�^�ɕύX
z.vec <- matrix(as.numeric(t(zv1)), nrow=hhpt*(member-1), ncol=g, byrow=T)
zi.v <- as.numeric(t(zv2))
id.v <- as.numeric(t(idv))
time.v <- as.numeric(t(tv))

#ID��ύX
ID.v <- data.frame(no=1:length(id.v), id=id.v, t=time.v, z=zi.v)
cbind(ID.v, z.vec)   #�f�[�^���m�F

##�����v���r�b�g���f�����D���ȃ����o�[�𔭐�
#��A�p�����[�^�̐ݒ�
beta0.z <- rbind(c(2.3, 2.0, 1.8, 0.7, -0.3, -0.8, -0.5, 1.2, 0.4),
                 c(0.6, -0.3, 0.8, 1.9, 2.3, 1.8, 0.2, -0.3, 0.6),
                 c(0.7, -0.3, -0.4, 1.0, 0.5, 0.6, 1.7, 2.0, 2.2))
beta1.z <- matrix(runif(g*2, 0, 1.2), nrow=g, ncol=2, byrow=T)
beta2.z <- matrix(runif(g*2, -1.2, 1.2), nrow=g, ncol=2, byrow=T)
beta3.z <- matrix(runif(g*(member-1), 0, 1.2), nrow=g, ncol=member-1, byrow=T)
beta4.z <- matrix(runif(g*(member-1), -1.2, 1.3), nrow=g, ncol=member-1, byrow=T)
beta.z <- t(cbind(beta0.z, beta1.z, beta2.z, beta3.z))   #�f�[�^������
rownames(beta.z) <- 1:nrow(beta.z)

#���U�����U�s���ݒ�
Cov <- corrM(member-1, -0.55, 0.95)   #���U�����U�p�����[�^�̓Z�O�����g�ŋ���


#���p�֐��̐ݒ�
U.mean <- matrix(rowSums(XM %*% beta.z * z.vec), nrow=hhpt, ncol=member-1, byrow=T)   #���p�֐��̕��ύ\��
U <- U.mean + mvrnorm(hhpt, rep(0, member-1), Cov)

#�ΐ��ޓx�̐^�̒l
inverse.Cov <- solve(Cov)
det.Cov <- det(Cov)
LLt <- apply(cbind(U, U.mean), 1, function(x) dmv(x[1:(member-1)], x[member:(2*(member-1))], inverse.Cov, det.Cov))
LL_true <- sum(log(LLt))


#���p�ő剻�����Ɋ�Â��I�������o�[����
y <- apply(U, 1, function(x) ifelse(max(x) < 0, member, which.max(x)))

#�I�������o�[��0�A1�s��ɕύX
Y <- matrix(0, nrow=hhpt, ncol=member)
for(i in 1:hhpt){
  Y[i, y[i]] <- 1
}


table(y)   #�I�������o�[�̏W�v
round(cbind(ID, y, U, U.mean), 2)   #�I�������o�[�ƌ��p���r

####�}���R�t�A�������e�J�����@�ɂ��K�w���������v���r�b�g���f���̐���####
####MCMC����̂��߂̐��菀��####
##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##���ϗʐ��K���z�̖ޓx�֐�
dmv <- function(x, mean.vec, inv.cov, det.cov){
  LLo <- det.cov^(-1/2) * exp(-1/2*(x - mean.vec) %*% inv.cov %*% (x - mean.vec))
  return(LLo)
}

##�������W�b�g���f���̖ޓx�֐�
LL_logit <- function(theta, Z, Zx, hh, g){
  #���W�b�g�̌v�Z
  logit <- matrix(Zx %*% theta, nrow=hh, ncol=g, byrow=T)
  
  #�m���Ƒΐ��ޓx�̘a���v�Z
  P <- exp(logit)/matrix(rowSums(exp(logit)), nrow=hh, ncol=g)
  LLi <- rowSums(Z * log(P))
  LL <- sum(LLi)
  val <- list(LL=LL, P=P)
  return(val)
}

LL_opt <- function(theta, Z, Zx, hh, g){
  #���W�b�g�̌v�Z
  logit <- matrix(Zx %*% theta, nrow=hh, ncol=g, byrow=T)
  
  #�m���Ƒΐ��ޓx�̘a���v�Z
  P <- exp(logit)/matrix(rowSums(exp(logit)), nrow=hh, ncol=g)
  LLi <- rowSums(Z * log(P))
  LL <- sum(LLi)
  return(LL)
}

##���ϗʐ��K���z�̏����t�����Ғl�ƕ��U���v�Z����֐�
cdMVN <- function(mu, Cov, dependent, U){
  
  #���U�����U�s��̃u���b�N�s����`
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent]
  Cov21 <- Cov[-dependent, dependent]
  Cov22 <- Cov[-dependent, -dependent]
  
  #�����t�����U�Ə����t�����ς��v�Z
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mu[, dependent] + t(CDinv %*% t(U[, -dependent] - mu[, -dependent]))   #�����t�����ς��v�Z
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #�����t�����U���v�Z
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##MCMC�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- array(0, dim=c(R/keep))   #�ΐ��ޓx�̕ۑ��p


##�����ϐ��𑽎����z��
X.array <- array(0, dim=c(member-1, ncol(X), hhpt))
selop <- matrix(1:nrow(XM), nrow=hhpt, ncol=member-1, byrow=T)

for(i in 1:hhpt){
  X.array[, , i] <- XM[selop[i, ], ]
}

YX.array <- array(0, dim=c(member-1, ncol(X)+1, hhpt))

#����v���Z�X�̊i�[�z��
UM <- matrix(0, nrow=hhpt, ncol=member-1)
util.M <- matrix(0, nrow=hhpt, ncol=member-1) 


##���O���z�̐ݒ�
nu <- member   #�t�E�B�V���[�g���z�̎��R�x
V <- solve((1/10)*diag(member-1))    #�t�E�B�V���[�g���z�̃p�����[�^
Deltabar <- rep(0, ncol(X))  #��A�W���̕��ς̎��O���z
Adelta <- solve(100 * diag(rep(1, ncol(X)))) #��A�W���̎��O���z�̕��U
zeta <- rep(0, ncol(Zx.v))   #�K�w���f���̉�A�W���̕��ς̎��O���z
Azeta <- solve(100 * diag(rep(1, ncol(Zx.v))))    #�K�w���f���̉�A�W���̎��O���z�̕��U

##�T���v�����O���ʂ̕ۑ��p�z��
Util <- array(0, dim=c(hhpt, member-1, R/keep))
BETA1 <- matrix(0, nrow=R/keep, ncol=ncol(XM))
BETA2 <- matrix(0, nrow=R/keep, ncol=ncol(XM))
BETA3 <- matrix(0, nrow=R/keep, ncol=ncol(XM))
SIGMA <- matrix(0, nrow=R/keep, ncol=nrow(Cov)^2)
THETA <- matrix(0, nrow=R/keep, ncol=ncol(Zx.v))
Z.VEC <- matrix(0, nrow=R/keep, ncol=hh)
Prob <- matrix(0, nrow=R/keep, ncol=hh)


##�����l�̐ݒ�
##�I�����f���̏����l�̐ݒ�
#�I�����f���̉�A�W���̏����l
beta00.z <- matrix(colSums(Y[, -member])/(hhpt/10), nrow=g, ncol=member-1, byrow=T) +
  matrix(runif((member-1)*g, -1, 1), nrow=g, ncol=member-1)
beta01.z <- matrix(runif(g*2, 0, 1.4), nrow=g, ncol=2, byrow=T)
beta02.z <- matrix(runif(g*2, -1.2, 1.3), nrow=g, ncol=2, byrow=T)
beta03.z <- matrix(runif(g*(member-1), 0, 1.2), nrow=g, ncol=member-1, byrow=T)
beta04.z <- matrix(runif(g*(member-1), -1.2, 1.2), nrow=g, ncol=member-1, byrow=T)
oldbeta <- t(cbind(beta00.z, beta01.z, beta02.z, beta03.z))   #�f�[�^������

#���U�����U�s��̏����l
oldcov <- corrM(member-1, 0, 0)
invcov <- solve(oldcov)

##���ݕϐ����f���̏����l�̐ݒ�
#theta�̏����l
oldtheta <- rep(0, ncol(Zx.v))

#���ݕϐ����f���̎��O�m���̌v�Z
Pr <- rep(1/g, g)

#�Z�O�����g�𐶐�
Z1 <- t(rmultinom(hh, 1, Pr))
z1 <- Z1 %*% 1:g
ZZ1 <- Z1

##���p�̏����l
#z���x�N�g���`���ɕύX
Zi1 <- matrix(0, nrow=nrow(Y), ncol=g)
Z.vec <- matrix(0, nrow=nrow(XM), ncol=g)
pt.zeros <- c(0, pt)

for(i in 1:hh){
  r1 <- (sum(pt.zeros[1:i])*(member-1)+1):(sum(pt.zeros[1:i])*(member-1)+pt[i]*(member-1))
  r2 <- (sum(pt.zeros[1:i])+1):(sum(pt.zeros[1:i])+pt[i])
  
  Z.vec[r1, ] <- matrix(Z1[i, ], nrow=pt[i]*(member-1), ncol=g, byrow=T)
  Zi1[r2, ] <- matrix(Z1[i, ], nrow=pt[i], ncol=g, byrow=T)
}

#���p���v�Z
old.utilm <- matrix(rowSums(XM %*% oldbeta * Z.vec), nrow=hhpt, ncol=member-1, byrow=T)   #���p�̕��ύ\��
old.util <- old.utilm + mvrnorm(hhpt, rep(0, member-1), oldcov)   #�덷��������


####�}���R�t�A�������e�J�����@�ŊK�w�L�����������v���r�b�g���f���𐄒�####
for(rp in 1:R){
  
  ##�I�����ʂƐ����I�Ȑ��݌��p�𔭐�������
  #�����t�����Ғl�Ə����t�����U���v�Z
  S <- rep(0, member-1)
  
  for(j in 1:(member-1)){
    MVR <- cdMVN(mu=old.utilm, Cov=oldcov, dependent=j, U=old.util)   #�����t�����z���v�Z
    UM[, j] <- MVR$CDmu   #�����t�����Ғl�����o��
    S[j] <- sqrt(MVR$CDvar) #�����t�����U�����o��
    
    #���݌��p�𔭐�������
    #�ؒf�̈���`
    max.u <- apply(cbind(old.util[, -j], 0), 1, max)
    max.u <- ifelse(y==member, 0, max.u)
    
    #�ؒf���K���z�����ݕϐ��𔭐�
    old.util[, j] <- ifelse(y==j, rtnorm(mu=UM[, j], sigma=S[j], a=max.u, b=100), 
                            rtnorm(mu=UM[, j], sigma=S[j], a=-100, b=max.u))
    old.util[, j] <- ifelse(is.infinite(old.util[, j]), ifelse(y==j, max.u + runif(1), max.u - runif(1)), old.util[, j])
  }
  util.v <- as.numeric(t(old.util))
  
  ##beta�̕��z�ƃp�����[�^�̌v�Z
  #z.vec��X.vec
  YX.bind <- cbind(util.v, XM)
  for(i in 1:hhpt){
    YX.array[, , i] <- YX.bind[id_r[i, ], ]
  }
  
  ##��A���f���̉�A�p�����[�^���Z�O�����g�ʂɃM�u�X�T���v�����O
  #�Z�O�����g�̃C���f�b�N�X���쐬
  zi1 <- Zi1 %*% 1:g
  index.zi1 <- list()
  
  for(i in 1:g){
    index.zi1[[i]] <- subset(1:length(zi1), zi1==i)
  }
  
  #�Z�O�����g���ƂɁ@beta�̕��ςƕ��U���v�Z
  invcov <- solve(oldcov)
  B <- matrix(0, nrow=ncol(XM), ncol=g)
  er <- matrix(0, nrow=nrow(XM), ncol=g)
  util.all <- matrix(0, nrow=nrow(XM), ncol=g)
  inv_XVX <- list()
  
  for(i in 1:g){
    xvx.vec <- rowSums(apply(X.array[, , index.zi1[[i]]], 3, function(x) t(x) %*% invcov %*% x))
    XVX <- matrix(xvx.vec, nrow=ncol(XM), ncol=ncol(XM), byrow=T)
    XVY <- rowSums(apply(YX.array[, , index.zi1[[i]]], 3, function(x) t(x[, -1]) %*% invcov %*% x[, 1]))
    
    #beta�̕��z�̕��U�����U�s��̃p�����[�^
    inv_XVX[[i]] <- solve(XVX + Adelta)
    
    #beta�̕��z�̕��σp�����[�^
    B[, i] <- inv_XVX[[i]] %*% (XVY + Adelta %*% Deltabar)   #beta�̕���
    
    #���ϗʐ��K���z�����A�W�����T���v�����O
    oldbeta[, i] <- mvrnorm(1, B[, i], inv_XVX[[i]])
  }
  
  #�덷���v�Z
  util.all <- XM %*% oldbeta
  er <- matrix(util.v, nrow=length(util.v), ncol=g) - util.all
  
  ##Cov�̕��z�̃p�����[�^�̌v�Z��mcmc�T���v�����O(�S�̂ŋ���)
  #�t�E�B�V���[�g���z�̃p�����[�^���v�Z
  R.error <- matrix(rowSums((er * Z.vec)), nrow=hhpt, ncol=member-1, byrow=T)
  Sn <- nu + hhpt
  IW.R <- V + matrix(rowSums(apply(R.error, 1, function(x) x %*% t(x))), nrow=member-1, ncol=member-1, byrow=T)
  
  #�t�E�B�V���[�g���z����Cov���T���v�����O
  oldcov <- rwishart(Sn, solve(IW.R))$IW
  
  ##���ϗʐ��K���z�̖ޓx�֐����Z�O�����g�ʂɌv�Z
  det_Cov_hat <- det(oldcov)
  inv_Cov_hat <- solve(oldcov)
  LLi <- matrix(0, nrow=hhpt, ncol=g)
  
  for(i in 1:g){
    util.seg <- matrix(util.all[, i], nrow=hhpt, ncol=member-1, byrow=T)
    util.bind <- cbind(old.util, util.seg)
    
    #���ϗʐ��K���z�̖ޓx���v�Z
    LLi[, i] <- apply(util.bind, 1, function(x) dmv(x[1:(member-1)], x[member:(2*(member-1))], inv_Cov_hat, det_Cov_hat))
  }
  
  #ID�ʂɖޓx�̐ς��v�Z
  LLind <- as.matrix(data.frame(id=ID$id, L=LLi) %>%
                       dplyr::group_by(id) %>%
                       dplyr::summarise_each(funs(prod), L.1, L.2, L.3))
  
  logl <- sum(log(rowSums(LLind[, -1] * Z1)))   #�ΐ��ޓx�̘a
  
  
  ##���݃Z�O�����g�𐶐�
  #���ݕϐ�z�̊����m��
  d1 <- exp(logit) * LLind[, -1]
  Pr_Z <- d1 / matrix(rowSums(d1), nrow=hh, ncol=g)
  
  #�Z�O�����g�𑽍����z��蔭��
  Z1 <- t(apply(Pr_Z, 1, function(x) rmultinom(1, 1, x)))
  z1 <- Z1 %*% 1:g
  
  #z���x�N�g���`���ɕύX
  Zi1 <- matrix(0, nrow=nrow(Y), ncol=g)
  Z.vec <- matrix(0, nrow=nrow(XM), ncol=g)
  pt.zeros <- c(0, pt)
  
  for(i in 1:hh){
    r1 <- (sum(pt.zeros[1:i])*(member-1)+1):(sum(pt.zeros[1:i])*(member-1)+pt[i]*(member-1))
    r2 <- (sum(pt.zeros[1:i])+1):(sum(pt.zeros[1:i])+pt[i])
    
    Z.vec[r1, ] <- matrix(Z1[i, ], nrow=pt[i]*(member-1), ncol=g, byrow=T)
    Zi1[r2, ] <- matrix(Z1[i, ], nrow=pt[i], ncol=g, byrow=T)
  }
  
  
  ##z1�̃Z�O�����g�������瑽�����W�b�g���f����theta���T���v�����O
  #MH�T���v�����O�ŌŒ����beta�̃T���v�����O
  if(rp==1 | rp==100 | rp==500 | rp%%1000==0){
    print("�œK��")
    res <- optim(oldtheta, LL_opt, gr=NULL, Z=Z1, Zx=Zx.v, hh=hh, g=g, method="BFGS", hessian=TRUE, 
                 control=list(fnscale=-1))
    oldtheta <- res$par
    rw <- -diag(solve(res$hessian))
    
  } else {
    
    newtheta <- oldtheta + rnorm(length(oldtheta), 0, rw)   #�����_���E�H�[�N
    
    #�����_���E�H�[�N�T���v�����O
    #�ΐ��ޓx�Ƒΐ����O���z���v�Z
    lognew <- LL_logit(theta=newtheta, Z=Z1, Zx=Zx.v, hh, g)$LL
    logold <- LL_logit(theta=oldtheta, Z=Z1, Zx=Zx.v, hh, g)$LL
    logpnew <- lndMvn(newtheta, zeta, Azeta)
    logpold <- lndMvn(oldtheta, zeta, Azeta)
    
    #MH�T���v�����O
    alpha <- min(1, exp(lognew + logpnew - logold - logpold))
    if(alpha == "NAN") alpha <- -1
    
    #��l�����𔭐�
    u <- runif(1)
    
    #u < alpha�Ȃ�V�����Œ����beta���̑�
    if(u < alpha){
      oldtheta <- newtheta
      
      #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
    } else {
      oldtheta <- oldtheta
    }
  }
  
  ##���p�֐����ƃ��W�b�g�̍X�V
  oldcov <- cov2cor(oldcov)
    
  old.utilm <- matrix(rowSums(util.all * Z.vec), nrow=hhpt, ncol=member-1, byrow=T)
  logit <- matrix(Zx.v %*% oldtheta, nrow=hh, ncol=g, byrow=T)
  
  ##�T���v�����O���ʂ̕ۑ��ƃp�����[�^�̊m�F
  if(rp%%keep==0){
    print(rp)
    mkeep <- rp/keep
    
    #�T���v�����O���ʂ�ۑ�
    BETA1[mkeep, ] <- oldbeta[, 1]
    BETA2[mkeep, ] <- oldbeta[, 2]
    BETA3[mkeep, ] <- oldbeta[, 3]
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    THETA[mkeep, ] <- oldtheta
    Z.VEC[mkeep, ] <- z1
    Util[, , mkeep] <- old.util 
    
    #�p�����[�^���m�F
    print(round(rbind(t(oldbeta), t(beta.z)), 2))
    print(round(cbind(cov2cor(oldcov), Cov), 2))
    print(round(rbind(oldtheta, theta.z), 2))
    print(round(c(colSums(Z1)/sum(Z1), r_rate), 2))
    print(logl)
    print(alpha)
  }
}

####���茋�ʂƗv��####
##�T���v�����O���ʂ��v���b�g
#��A�W���̃v���b�g
matplot(BETA1[, 1:4], type="l", ylab="beta1�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA1[, 5:9], type="l", ylab="beta1�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA1[, 10:13], type="l", ylab="beta1�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA1[, 14:17], type="l", ylab="beta1�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA1[, 18:22], type="l", ylab="beta1�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA2[, 1:4], type="l", ylab="beta2�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA2[, 5:9], type="l", ylab="beta2�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA2[, 10:13], type="l", ylab="beta2�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA2[, 14:17], type="l", ylab="beta2�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA2[, 18:22], type="l", ylab="beta2�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA3[, 1:4], type="l", ylab="beta3�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA3[, 5:9], type="l", ylab="beta3�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA3[, 10:13], type="l", ylab="beta3�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA3[, 14:17], type="l", ylab="beta3�̉�A�W��", xlab="�T���v�����O��")
matplot(BETA3[, 18:22], type="l", ylab="beta3�̉�A�W��", xlab="�T���v�����O��")

#���U�����U�s��̃v���b�g
matplot(SIGMA[, 1:9], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 10:18], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 19:27], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 28:36], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 37:45], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 46:54], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 55:63], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 64:72], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")
matplot(SIGMA[, 73:81], type="l", ylab="sigma�̃p�����[�^", xlab="�T���v�����O��")

#�K�w���f���̉�A�W���̃p�����[�^
matplot(THETA[, 1:4], type="l", ylab="theta�̉�A�W��", xlab="�T���v�����O��")
matplot(THETA[, 5:8], type="l", ylab="theta�̉�A�W��", xlab="�T���v�����O��")
matplot(THETA[, 9:14], type="l", ylab="theta�̉�A�W��", xlab="�T���v�����O��")



