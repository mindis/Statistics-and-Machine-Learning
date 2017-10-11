#####�x�C�W�A�����ݐ��ڑ��ϗʃv���r�b�g���f��#####
library(MASS)
library(mlogit)
library(MCMCpack)
library(bayesm)
library(mvtnorm)
library(caret)
library(reshape2)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(lattice)

#set.seed(4267)

####�C�ӂ̕��U�����U�s����쐬������֐�####
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  X.Sigma <- eigen(Sigma)
  Lambda <- diag(X.Sigma$values)
  P <- X.Sigma$vector
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda)
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
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
pt <- 6
hhpt <- hh*pt
seg <- 4
choise <- 8

##ID�̐ݒ�
id <- rep(1:hh, rep(pt, hh))
time <- rep(1:pt, hh)
ID <- data.frame(no=1:hhpt, id, time)


####�����ϐ��̔���####
PRICE <- matrix(runif(hhpt*choise, 0.5, 1), nrow=hhpt, ncol=choise) - 1
DISC <- matrix(runif(hhpt*choise, 0, 0.5), nrow=hhpt, ncol=choise)

DISP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hh, 1, r)
}

CAMP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hh, 1, r)
}

income <- exp(rnorm(hh, 1.78, 0.1))
income <- scale(income)


INCOME <- rep(income, rep(pt, hh))
hist(exp(INCOME), breaks=20, col="grey", xlab="income", main="�����̕��z")


##�����ϐ����x�N�g���`���ɕϊ�
#��A�W�����S�u�����h�ŋ��ʂ̐����ϐ�
DISP.vec <- as.numeric(t(DISP))
CAMP.vec <- as.numeric(t(CAMP))

#��A�W�����u�����h�ňقȂ�����ϐ�
BP.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
PRICE.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
DISC.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
INCOME.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)

for(i in 1:hhpt){
  r <- ((i-1)*choise+1):((i-1)*choise+choise)
  BP.vec[r, ] <- diag(choise) 
  PRICE.vec[r, ] <- diag(PRICE[i, ])
  DISC.vec[r, ] <- diag(DISC[i, ])
  INCOME.vec[r, ] <- diag(INCOME[i], choise)
}

X.vec <- data.frame(bp=BP.vec, price=PRICE.vec, disc=DISC.vec, disp=DISP.vec, camp=CAMP.vec, income=INCOME.vec)


####�Z�O�����g�̐ݒ�####
#�Z�O�����g���ڊm���s��̐ݒ�
P_seg0 <- matrix(0, nrow=seg, ncol=seg)
for(i in 1:seg){
  rand <- runif(seg, 0.1, 3.5)
  P_seg0[i, ] <- rand
}

diag(P_seg0) <- runif(seg, 7.5, 25.0)   #�Ίp�s���u������
P_seg <- P_seg0 / matrix(rowSums(P_seg0), nrow=seg, ncol=seg)   #�m���ɒu������


#�Z�O�����g�𔭐�������
#����1�̃Z�O�����g�̐ݒ�
seg_m <- matrix(0, nrow=hh, ncol=pt)
seg_m[, 1] <- rep(1:seg, rep(hh/seg, seg))

#����2�`7�܂Œ����I�ɃZ�O�����g�𔭐�������
for(j in 2:pt){
  for(i in 1:hh){
    seg_m[i, j] <- t(rmultinom(1, 1, P_seg[seg_m[i, j-1], ])) %*% 1:seg
  }
}

seg_v <- as.numeric(t(seg_m))   #�Z�O�����g���x�N�g���ɕύX
table(seg_v)
seg_m <- matrix(as.numeric(table(1:hhpt, seg_v)), nrow=hhpt, ncol=seg)

r_rate0 <- do.call(rbind, tapply(seg_v, ID$time, function(x) table(x)/sum(table(x))))   #������
r_rate0 <- matrix(as.numeric(r_rate0), nrow=pt, ncol=seg)


##ID�̐ݒ�
id.v <- rep(1:hh, rep(choise*pt, hh))
pd <- rep(1:choise, hhpt)
t.vec <- rep(rep(1:pt, rep(choise, pt)), hh)
idno <- rep(1:hhpt, rep(choise, hhpt))
ID.vec <- data.frame(no=1:(hhpt*choise), idno=idno, id=id.v, t=t.vec, pd=pd)

seg_vec <- rep(seg_v, rep(choise, hhpt))


####�Z�O�����g�������牞���ϐ��𔭐�####
##�p�����[�^�̐ݒ�
#���֍s��̐ݒ�
corM <- corrM(col=choise, lower=-0.5, upper=0.8, eigen_lower=0.01, eigen_upper=0.2)   #���֍s����쐬
Sigma <- covmatrix(col=choise, corM=corM, lower=1, upper=1)   #���U�����U�s��
Cov <- Sigma$covariance

##�Ó��ȉ����ϐ�����������܂ŌJ��Ԃ�
for(i in 1:10000){
  print(i)
  
  #��A�p�����[�^�̐ݒ�
  beta0 <- matrix(c(runif(choise*seg, -2.4, 1.2)), nrow=seg, ncol=choise)
  beta1 <- matrix(c(runif(choise*seg, -2.2, -0.5)), nrow=seg, ncol=choise)
  beta2 <- matrix(c(runif(choise*seg, 0.6, 2.0)), nrow=seg, ncol=choise)
  beta3 <- matrix(c(runif(seg, 0.1, 1.4)), nrow=seg, ncol=1)
  beta4 <- matrix(c(runif(seg, 0.2, 1.5)), nrow=seg, ncol=1)
  beta5 <- matrix(c(runif(choise*seg, -0.7, 0.7)), nrow=seg, ncol=choise)
  betat <- cbind(beta0, beta1, beta2, beta3, beta4, beta5)   #��A�W��������
  
  ##���p�֐��̌v�Z�Ɖ����ϐ��̔���
  #�Z�O�����g���ƂɌ��p�֐��̌v�Z
  U.mean <- matrix(0, nrow=hhpt, ncol=choise)
  for(j in 1:seg){
    U.mean[seg_v==j, ] <- matrix(as.matrix(X.vec[seg_vec==j, ]) %*% betat[j, ], nrow=sum(seg_v==j), ncol=choise, byrow=T)
  }
  
  error  <- mvrnorm(hhpt, rep(0, choise), Cov)
  U <- U.mean + error
  
  
  #�����ϐ��̔���
  Y <- ifelse(U > 0, 1, 0)
  print(c(min(colMeans(Y)), max(colMeans(Y))))
  if(min(colMeans(Y)) > 0.2 & max(colMeans(Y)) < 0.75) break
}


####�}���R�t�A�������e�J�����@�Ő��ݐ��ڑ��ϗʃv���r�b�g���f���𐄒�####
##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##���ϗʐ��K���z�̏����t�����Ғl�Ə����t�����U���v�Z����֐�
cdMVN <- function(mean, Cov, dependent, U){
  
  #���U�����U�s��̃u���b�N�s����`
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent, drop=FALSE]
  Cov21 <- Cov[-dependent, dependent, drop=FALSE]
  Cov22 <- Cov[-dependent, -dependent]
  
  #�����t�����U�Ə����t�����ς��v�Z
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mean[, dependent] + t(CDinv %*% t(U[, -dependent] - mean[, -dependent]))   #�����t�����ς��v�Z
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #�����t�����U���v�Z
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##MCMC�̐ݒ�####
R <- 20000
sbeta <- 1.5
keep <- 4

##���O���z�̐ݒ�
nu <- choise   #�t�E�B�V���[�g���z�̎��R�x
V <- nu*diag(choise)   #�t�E�B�V���[�g���z�̃p�����[�^
Deltabar <- rep(0, ncol(X.vec))   #��A�W���̎��O���z�̕���
Adelta <- solve(100 * diag(rep(1, ncol(X.vec)))) #��A�W���̎��O���z�̕��U

#�T���v�����O���ʂ̕ۑ��p�z��
BETA <- matrix(0, nrow=R/keep, ncol=ncol(X.vec)*seg)
SIGMA <- matrix(0, nrow=R/keep, ncol=choise*choise)
Z <- matrix(0, nrow=R/keep, ncol=hhpt)


#�f�[�^�̐ݒ�
X.vec <- as.matrix(X.vec)
id_r <- matrix(1:(hhpt*choise), nrow=hhpt, ncol=choise, byrow=T)

#�v�Z�p�p�����[�^�̊i�[�p
U <- matrix(0, nrow=hhpt, ncol=choise)   #���p�֐��̊i�[�p
YX.array <- array(0, dim=c(choise, ncol(X.vec)+1, hhpt))
MVR.U <- matrix(0, nrow=hhpt, ncol=choise)
old.util <- rep(0, hhpt*choise)


##�����l�̐ݒ�
#��A�W���̏����l
#�v���r�b�g���f���̑ΐ��ޓx�̒�`
probit_LL <- function(x, Y, X){
  #�p�����[�^�̐ݒ�
  b0 <- x[1]
  b1 <- x[2:(ncol(X)+1)]
  
  #���p�֐��̒�`
  U <- b0 + as.matrix(X) %*% b1
  
  #�ΐ��ޓx���v�Z
  Pr <- pnorm(U)   #�m���̌v�Z
  LLi <- Y*log(Pr) + (1-Y)*log(1-Pr)
  LL <- sum(LLi)
  return(LL)
}

#�����ϐ����ƂɓƗ��Ƀv���r�b�g���f���𓖂Ă͂ߏ����l��ݒ�
first_beta <- c()
for(b in 1:choise){
  print(b)
  
  #�����p�����[�^�̐ݒ�
  X <- cbind(PRICE[, b], DISC[, b], DISP[ID.vec$pd==b], CAMP[ID.vec$pd==b], INCOME)
  x <- c(runif(1, -1.0, 1.0), runif(1, -1.8, -1.0), runif(1, 1.0, 1.8), runif(1, 0.5, 1.0), 
         runif(1, 0.5, 1.0), runif(1, -0.1, 0.1))
  
  #���j���[�g���@�ōő剻
  res <- optim(x, probit_LL, Y=Y[, b], X=X, method="BFGS", hessian=FALSE, 
               control=list(fnscale=-1))
  first_beta <- rbind(first_beta, res$par)
}

#�����ϐ��ɓK������悤�ɏ����l���N�����W���O
oldbeta0 <- rep(0, ncol(X.vec))
index_name <- subset(1:ncol(X.vec), colnames(X.vec) %in% c("disp", "camp"))
oldbeta0[index_name] <- colSums(first_beta[, 4:5])/choise
oldbeta0[-index_name] <- as.numeric(first_beta[, -c(4:5)])

#�Z�O�����g�ʂɏ����l�𔭐�
rw <- matrix(rnorm(length(oldbeta0)*seg, 0, 0.25), nrow=seg, ncol=length(oldbeta0))
oldbeta <- matrix(oldbeta0, nrow=seg, ncol=length(oldbeta0), byrow=T) + rw

#���֍s��̏����l
oldcov <- cor(Y)

#�Z�O�����g�����̏����l
LLind <- matrix(0, nrow=hhpt, ncol=seg)

for(j in 1:seg){
  U <- matrix(X.vec %*% oldbeta[j, ], nrow=hhpt, ncol=choise, byrow=T)
  LH <- pnorm(U)^Y * pnorm(1-U)^(1-Y)
  LLind[, j] <- rowProds(LH)
}


#���ݕϐ�z�̊����m�����v�Z
#���Ԃ��Ƃ̍������̏����l��ݒ�
r <- matrix(0, nrow=hhpt, ncol=seg)
r[ID$time==1, ] <- 1/seg
r[ID$time %in% 2:pt, ] <- LLind[ID$time %in% 1:(pt-1), ]

z_rate <- r*LLind / matrix(rowSums(r*LLind), nrow=hhpt, ncol=seg)   #�Z�O�����g�����m��

#�Z�O�����g�����𑽍����z���甭��
z <- t(apply(z_rate, 1, function(x) rmultinom(1, 1, x)))   #���ݕϐ�z�̔���

##�C���f�b�N�X��ݒ�
#�Z�O�����g�����̃C���f�b�N�X
index_z  <- list()
index_zind <- list()

for(j in 1:seg){
  index_z[[j]] <- subset(1:hhpt, z[, j]==1)
  index_zind[[j]] <- subset(1:nrow(ID.vec), ID.vec$idno %in% index_z[[j]])
}

#���Ԃ̃C���f�b�N�X
index_time <- list()
index_time[[1]] <- subset(1:hhpt, ID$time==1)
index_time[[2]] <- subset(1:hhpt, ID$time %in% 2:pt)
index_time[[3]] <- subset(1:hhpt, ID$time %in% 1:(pt-1))


##���݌��p�̏����l
old.utilm <- matrix(0, nrow=hhpt, ncol=choise)   #���p�̕��ύ\��
U <- matrix(0, nrow=hhpt, ncol=choise)   #���݌��p
er <- mvrnorm(hhpt, rep(0, choise), oldcov)   #���݌��p�̌덷

for(j in 1:seg){
  old.utilm[index_z[[j]], ] <- matrix(X.vec[index_z[[j]], ] %*% oldbeta[j, ], nrow=length(index_z[[j]]), ncol=choise, byrow=T) 
  U[index_z[[j]], ] <- old.utilm[index_z[[j]], ] + er[index_z[[j]], ]
}
old.util <- as.numeric(t(old.utilm))

#�ؒf���K���z�̐ؒf�̈���`
a <- ifelse(Y==0, -100, 0)
b <- ifelse(Y==1, 100, 0)

####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�I�����ʂƐ����I�Ȑ��ݕϐ��𔭐�������
  for(s in 1:seg){
    #�Z�O�����g�ʂɃf�[�^�𒊏o
    X.seg <- X.vec[index_zind[[s]], ]
    index <- index_z[[s]]
    
    #���p�̕��ύ\��
    u_mu <- old.utilm[index, ]
    
    #�ؒf���K���z�����݌��p�𔭐�
    for(j in 1:choise){
      MVR <- cdMVN(u_mu, oldcov, j, U[index, ])
      MVR.U <- MVR$CDmu
      MVR.S <- sqrt(MVR$CDvar)
      U[index, j] <- rtnorm(MVR.U, MVR.S, a[index, j], b[index, j])
    }
    U[is.nan(U)] <- 0
    
    ##beta�̕��z�̃p�����[�^���T���v�����O
    u_vec <- as.numeric(t(U[index, ]))   #���݌��p���x�N�g����
    
    #beta�̕��σp�����[�^���v�Z
    XX <- t(X.seg) %*% X.seg
    XY <- t(X.seg) %*% u_vec
    
    inv_XVX <- solve(XX + Adelta)
    beta <- as.numeric(inv_XVX %*% (XY + Adelta %*% Deltabar))   #��A�W���̕���
    
    #���ϗʐ��K���z����beta���T���v�����O
    oldbeta[s, ] <- mvrnorm(1, beta, inv_XVX)
    old.util[index_zind[[s]]] <- X.seg %*% oldbeta[s, ]
    
    #�����ȉ�A�W���̃p�����[�^�v�Z
    #XVX <- matrix(0, nrow=ncol(X.vec), ncol=ncol(X.vec))
    #XVY <- rep(0, ncol(X.vec))
    
    #for(i in 1:length(index)){
    #  num <- ((i-1)*choise+1):((i-1)*choise+choise)
    #  XVX <- XVX + t(X.seg[num, ]) %*% inv_cov %*% X.seg[num, ]
    #  XVY <- XVY + t(X.seg[num, ]) %*% inv_cov %*% u_vec[num]
    #}
    

  }
  
  ##���ʂ̕��U�����U�s����T���v�����O
  u <- as.numeric(t(U))   #���݌��p���x�N�g����
  old.utilm <- matrix(old.util, nrow=hhpt, ncol=choise, byrow=T)
  
  #�t�E�B�V���[�g���z�̃p�����[�^���v�Z
  R.error <- U - old.utilm
  IW.R <- solve(V) + t(R.error) %*% R.error
  
  #�t�E�B�V���[�g���z�̎��R�x���v�Z
  Sn <- nu + hhpt
  
  #�t�E�B�V���[�g���z�̎��R�x���v�Z
  Cov_hat <- rwishart(Sn, solve(IW.R))$IW
  
  
  ##���ʐ��̖���������邽�߂ɕ��U�����U�s��̑Ίp������1�ɌŒ肷��
  lambda <- diag(diag(Cov_hat)^(-1/2))
  oldcov <- cov2cor(Cov_hat)
  inv_cov <- solve(oldcov)
  old.utilm <- old.utilm %*% lambda
  U <- U %*% lambda
  
  ##���ݕϐ�����Z�O�����gz�𔭐�
  #�Z�O�����g���Ƃ̖ޓx���v�Z
  r <- matrix(0, nrow=hhpt, ncol=seg)
  r[index_time[[1]], ] <- 1
  LLind <- matrix(0, nrow=hhpt, ncol=seg)
  
  for(j in 1:seg){
    Xb <- matrix(X.vec %*% oldbeta[j, ], nrow=hhpt, ncol=choise, byrow=T) %*% lambda
    LH <- matrix(dnorm(as.numeric(U), as.numeric(Xb), 1), nrow=hhpt, ncol=choise)
    LLind[, j] <- rowProds(LH)
  }
  
  #���ݕϐ�z�̊����m�����v�Z
  #���Ԃ��Ƃ̍������̏����l��ݒ�
  #r[index_time[[2]], ] <- LLind[index_time[[3]], ]
  
  ##���O���z�����Ԃ��Ƃɋ��ʂ̏ꍇ
  r_table <- matrix(as.numeric(table(ID$time, z %*% 1:seg)), nrow=pt, ncol=seg)
  r_rate <- (r_table / matrix(hh, nrow=pt, ncol=seg))[1:(pt-1), ]
  for(i in 1:(pt-1)) {r[ID$time==i+1, ] <- matrix(r_rate[i, ], nrow=hh, ncol=seg, byrow=T)}

  #�Z�O�����g�����𑽍����z���甭��
  z_rate <- r*LLind / matrix(rowSums(r*LLind), nrow=hhpt, ncol=seg)   #�Z�O�����g�����m��
  z <- t(apply(z_rate, 1, function(x) rmultinom(1, 1, x)))   #���ݕϐ�z�̔���
  
  ##�C���f�b�N�X��ݒ�
  #�Z�O�����g�����̃C���f�b�N�X
  index_z  <- list()
  index_zind <- list()
  
  for(j in 1:seg){
    index_z[[j]] <- subset(1:hhpt, z[, j]==1)
    index_zind[[j]] <- subset(1:nrow(ID.vec), ID.vec$idno %in% index_z[[j]])
  }
  
  ##�T���v�����O���ʂ�ۑ�
  mkeep <- rp/keep
  if(rp%%keep==0){
    keep_er <- mkeep
    BETA[mkeep, ] <- as.numeric(t(oldbeta))
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    Z[mkeep, ] <- as.numeric(z %*% 1:seg)
    print(rp)
    print(round(rbind(oldbeta[, 1:20], betat[, 1:20]), 2))
    print(round(cbind(oldcov, Cov), 2))
  }
}

matplot(BETA[, 1:5], type="l")

data <- round(cbind(ID$id, z %*% 1:seg, seg_v, z_rate), 3)
