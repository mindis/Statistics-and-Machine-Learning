#####�U���^���ϗʗ��U�E�A�����f��#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(extraDistr)
library(tmvtnorm)
library(gtools)
library(MNP)
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
items <- 10   #�ϑ��u�����h��


####�����ϐ��̔���####
#�A���ϐ��̔���
cont <- 2
X.cont <- matrix(rnorm(hh*cont, 0, 1), nrow=hh, ncol=cont)

#��l�ϐ��̔���
bin <- 2
X.bin <- matrix(0, nrow=hh, ncol=bin)
for(i in 1:bin){
  X.bin[, i] <- rbinom(hh, 1, runif(1, 0.35, 0.65))
}



#���i�ƃv�����[�V���������ϐ�
Price <- matrix(0, nrow=hh, ncol=items)
Disp <- matrix(0, nrow=hh, ncol=items)
freq <- rpois(hh, 5)
freq <- ifelse(freq < 1, 1, freq)

for(i in 1:items){
  lower <- runif(1, 0.4, 0.7)
  upper <- runif(1, 0.8, 1.0)
  Price[, i] <- runif(hh, lower, upper) 
  Disp[, i] <- rbinom(hh, freq, runif(1, 0.05, 0.6))/freq
}


##�����ϐ����x�N�g����
#�ؕЂƐ����ϐ��̐ݒ�
BP.vec <- matrix(0, nrow=hh*items, ncol=items)
X1.vec <- matrix(0, nrow=hh*items, ncol=items*cont)
X2.vec <- matrix(0, nrow=hh*items, ncol=items*bin)

for(i in 1:hh){
  r <- ((i-1)*items+1):((i-1)*items+items)
  BP.vec[r, ] <- diag(items) 
  X1.vec[r, ] <- cbind(diag(X.cont[i, 1], items), diag(X.cont[i, cont], items))
  X2.vec[r, ] <- cbind(diag(X.bin[i, 1], items), diag(X.bin[i, bin], items))
}

#���i�ƃv�����[�V�����̐ݒ�
Price.v <- as.numeric(t(Price))
Disp.v <- as.numeric(t(Disp))


#�����ϐ��̌���
X.vec <- data.frame(bp=BP.vec, cont=X1.vec, bin=X2.vec, price=Price.v, disp=Disp.v)
XM.vec <- as.matrix(X.vec)
k1 <- ncol(XM.vec)

##ID�̐ݒ�
id.v <- rep(1:hh, rep(items, hh))
pd <- rep(1:items, hh)
ID.vec <- data.frame(no=1:(hh*items), id=id.v, pd=pd)


####�w���L���𑽕ϗʃv���r�b�g���f���Ŕ���������####
#���֍s��̐ݒ�
corM <- corrM(col=items, lower=-0.55, upper=0.9, eigen_lower=0.025, eigen_upper=0.35)   #���֍s����쐬
Sigma <- covmatrix(col=items, corM=corM, lower=1, upper=1)   #���U�����U�s��
Cov <- Sigma$covariance


#�Ó��ȍw�������o��܂Ńp�����[�^�̐ݒ���J��Ԃ�
for(i in 1:10000){
  print(i)
  
  #�p�����[�^�̐ݒ�
  beta0 <- runif(items, -0.7, 2.0)
  beta1 <- runif(items*cont, 0, 0.9)
  beta2 <- runif(items*bin, -0.8, 1.0)
  beta3 <- runif(1, -3.5, -2.5)
  beta4 <- runif(1, 2.0, 2.5)
  betat <- c(beta0, beta1, beta2, beta3, beta4)
  
  ##���p�֐��̌v�Z�Ɖ����ϐ��̔���
  #���p�֐��̌v�Z
  U_mean.vec <- cbind(XM.vec) %*% betat
  error.vec  <- as.numeric(t(mvrnorm(hh, rep(0, items), Cov)))
  U.vec <- U_mean.vec + error.vec
  U.mx <- matrix(U.vec, nrow=hh, ncol=items, byrow=T)
  
  #�����ϐ��̔���
  y1 <- ifelse(U.vec > 0, 1, 0)
  Y1 <- matrix(y1, nrow=hh, ncol=items, byrow=T)
  
  #�Ó��ȍw�������o��܂ŌJ��Ԃ�
  if(min(colMeans(Y1)) > 0.15 & max(colMeans(Y1)) < 0.75) break   
}



#���������������ϐ����W�v
colMeans(Y1); colSums(Y1)

####�w�����������ꍇ�A�|�A�\�����z����w�����𔭐�####
U.pois <- ifelse(U.mx < 0, 0, U.mx)   #�|�A�\�����z�̌��p���`

##�w���m���ƌ��p�̊֌W
mu <- c()
for(i in 1:items){
  mu <- c(mu, mean(U.pois[U.pois[, i] > 0, i]))
}

round(rbind(u=mu, pr=colMeans(Y1)), 3)


##�|�A�\�����z�̐����ϐ�
#�w�����͉Ƒ��l���Ǝ����Ɉˑ�
family <- rpois(hh, 2.5)
income <- rpois(hh, 5.3)

#�[����1�ɒu������
family <- ifelse(family < 1, 1, family)
income <- ifelse(income < 1, 1, income)

##�����ϐ����x�N�g���`���ɕύX
family.v <- as.numeric(t(matrix(family, nrow=hh, ncol=items)))
income.v <- as.numeric(t(matrix(income, nrow=hh, ncol=items)))
upois.v <- as.numeric(t(U.pois))

#�w�����������Ă��Ȃ��Ȃ�income��family��0�ɂ��Ă���
income.vec <- income.v * y1
family.vec <- family.v * y1

#�����ϐ�������
Z <- data.frame(u=upois.v, income=income.vec, family=family.vec)
ZM <- as.matrix(Z)
k2 <- ncol(ZM)


##�|�A�\�����z����w�����𔭐�
#�p�����[�^�̐ݒ�
theta1 <- runif(1, 0.4, 0.8)
theta2 <- runif(1, 0.05, 0.08)
theta3 <- runif(1, 0.08, 0.13)
thetat <- c(theta1, theta2, theta3)

#lambda���v�Z
lambda <- exp(ZM %*% thetat) * y1

#�|�A�\�����z�����牞���ϐ��𔭐�
y2 <- rep(0, length(lambda))

for(i in 1:length(lambda)){
  print(i)
  if(lambda[i]==0) next   #lambda��0�Ȃ玟��
  
  for(j in 1:1000){
    y2[i] <- rpois(1, lambda[i])
    if(y2[i]==0) {next} else {break}
  }
}

Y2 <- matrix(y2, nrow=hh, ncol=items, byrow=T)   #�s��`���ɕύX


##���ʂ̏W�v�Ɗm�F
summary(Y2)

Y2.mean <- c()
for(i in 1:items){
  Y2.mean <- c(Y2.mean, mean(Y2[Y2[, i] > 0, i]))
}

round(rbind(u=mu, pr=colMeans(Y1), buy=Y2.mean), 3)   #���p�A�w���m���A�w�����ʂ̊֌W


####�}���R�t�A�������e�J�����@�ő��ϗʗ��U�A�����f���𐄒�####
##�ؒf���K���z�̗����𔭐�������֐�
trunc_rnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}


##���ϗʐ��K���z�̏����t�����Ғl�Ə����t�����U���v�Z����֐�
cdMVN <- function(mu, Cov, dependent, U){
  
  #���U�����U�s��̃u���b�N�s����`
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent]
  Cov21 <- Cov[-dependent, dependent]
  Cov22 <- Cov[-dependent, -dependent]
  
  #�����t�����U�Ə����t�����ς��v�Z
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mu[, dependent] + t(CDinv %*% t(U[, -dependent] - mu[, -dependent]))   #�����t�����ς��v�Z
  CDvar <- 1 - CDinv %*% Cov21   #�����t�����U���v�Z
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}


##�|�A�\����A���f���̑ΐ��ޓx
loglike <- function(theta, y, X){
  
  #�ޓx���`����
  lambda <- exp(X %*% theta)   #���ύ\��
  LLi <- y*log(lambda)-lambda - lfactorial(y)   #�ΐ��ޓx
  LL <- sum(LLi)   #�ΐ��ޓx�̘a
  
  #���ʂ�Ԃ�
  LL.val <- list(LLi=LLi, LL=LL)
  return(LL.val)
}


####MCMC�A���S���Y���̐ݒ�####
R <- 20000
sbeta <- 1.5
keep <- 4
non_iden <- 5   #���ʂ���Ȃ��p�����[�^�̌�

#�^�̏����t�����U
cd_Cov <- c()
for(j in 1:items){
  cd_Cov <- c(cd_Cov, cdMVN(matrix(U_mean.vec, nrow=hh, ncol=items, byrow=T), Cov, j, U.mx)$CDvar)
}


##���O���z�̐ݒ�
nu <- items   #�t�E�B�V���[�g���z�̎��R�x
V <- solve(nu*diag(items))   #�t�E�B�V���[�g���z�̃p�����[�^
Deltabar <- rep(0, ncol(XM.vec))   #���ϗʃv���r�b�g���f���̉�A�W���̎��O���z�̕���
Adelta <- solve(100 * diag(rep(1, k1)))   #���ϗʃv���r�b�g���f���̉�A�W���̎��O���z�̕��U
Zetabar <- rep(0, k2)   #�|�A�\����A���f���̉�A�W���̎��O���z
Azeta <- solve(100 * diag(rep(1, k2)))   #�|�A�\����A���f���̉�A�W���̎��O���z

##�T���v�����O���ʂ̕ۑ��p�z��
Util <- matrix(0, nrow=R/keep, ncol=items)
BETA <- matrix(0, nrow=R/keep, k1)
SIGMA <- matrix(0, nrow=R/keep, ncol=items^2)
THETA <- matrix(0, nrow=R/keep, ncol=k2)


##�f�[�^�̐ݒ�
#�f�U�C���s��𑽎����z��ɕύX
X.array <- array(0, dim=c(items, k1, hh))
for(i in 1:hh){
  X.array[, , i] <- XM.vec[ID.vec$id==i, ]
}
YX.array <- array(0, dim=c(items, k1+1, hh))

#�C���f�b�N�X�̍쐬
id_r <- matrix(1:(hh*items), nrow=hh, ncol=items, byrow=T)


##�v�Z�p�p�����[�^�̊i�[�p
Z <- matrix(0, nrow=hh, ncol=items)   #���p�֐��̊i�[�p
MVR.U <- matrix(0, nrow=hh, ncol=items)
MVR.S <- rep(0, items)



##�����l�̐ݒ�
#��A�W���̏����l
##�v���r�b�g���f���̑ΐ��ޓx�̒�`
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
first_beta <- matrix(0, nrow=items, ncol=(k1-2)/items+2)

for(b in 1:items){
  print(b)
  for(i in 1:1000){
    #�����p�����[�^�̐ݒ�
    X <- cbind(X.cont, X.bin, Price[, b], Disp[, b])
    x <- c(runif(1, -0.5, 1.0), runif(cont, 0, 1), runif(bin, -0.5, 0.5), runif(1, -3.5, -2.5), runif(1, 2.0, 2.5))
    
    #���j���[�g���@�ōő剻
    res <- try(optim(x, probit_LL, Y=Y1[, b], X=X, method="BFGS", hessian=FALSE, 
                     control=list(fnscale=-1)), silent=TRUE)
    
    #�G���[����
    if(class(res) == "try-error") {
      next
    } else {
      first_beta[b, ] <- res$par
      break
    }   
  }
}

oldbeta <- c(as.numeric(first_beta[, 1:(1+cont+bin)]), mean(first_beta[, 2+cont+bin]), mean(first_beta[, 3+cont+bin]))
betaf <- oldbeta

#�|�A�\����A�̏����l
oldtheta <- c(runif(1, 0.4, 0.8), runif(2, 0, 0.1))
thetaf <- oldtheta

#���U�����U�s��̏����l
corf <- corrM(col=items, lower=0, upper=0)   #���֍s����쐬
Sigmaf <- covmatrix(col=items, corM=corf, lower=1, upper=1)   #���U�����U�s��
oldcov <- Sigmaf$covariance

#���݌��p�̏����l
old.utilm<- matrix(XM.vec %*% oldbeta, nrow=hh, ncol=items, byrow=T)   #���݌��p�̕��ύ\��
Z <- old.utilm + mvrnorm(hh, rep(0, items), oldcov)   #���ύ\��+�덷

#�ؒf���K���z�̐ؒf�̈���`
a <- ifelse(Y1==0, -100, 0)
b <- ifelse(Y1==1, 100, 0)

####�}���R�t�A�������e�J�����@�ő��ϗʗ��U�A�����f���𐄒�####
for(rp in 1:R){
  
  ##�f�[�^�g��@�Ő��݌��p�𔭐�
  #�I�����ʂƐ����I�Ȑ��݌��p�𔭐�������
  #�ؒf���K���z�����݌��p�𔭐�
  for(j in 1:items){
    MVR <- cdMVN(old.utilm, oldcov, j, Z)   #���ϗʐ��K���z�̏����t�����z���v�Z
    MVR.U[, j] <- MVR$CDmu   #�����t�����ς𒊏o
    MVR.S[j] <- sqrt(MVR$CDvar)   #�����t�����U�𒊏o
    Z[, j] <- trunc_rnorm(MVR.U[, j], MVR.S[j], a[, j], b[, j])   #�ؒf���K���z�����ݕϐ����T���v�����O
    
    #���݌��pZ�ɃG���[���������ꍇ�̏���
    if(sum(is.infinite(Z[, j])) > 0 | sum(is.nan(Z[, j]) > 0)){
      print("�G���[")
      index.error <- subset(1:nrow(Z), is.infinite(Z[, j])==TRUE | is.nan(Z[, j])==TRUE)
      Z[index.error, ] <- 0
    }
  }
  
  z.vec <- as.numeric(t(Z))   #���݌��p���x�N�g���ɕύX
  
  
  ##beta�̕��z�̃p�����[�^�̌v�Z��mcmc�T���v�����O
  #beta�̕��σp�����[�^���v�Z
  XVX <- matrix(0, nrow=ncol(XM.vec), ncol=ncol(XM.vec))
  XVY <- rep(0, ncol(XM.vec))
  invcov <- solve(oldcov)
  
  for(i in 1:hh){
    XVX <- XVX + t(X.array[, , i]) %*% invcov %*% X.array[, , i]
    XVY <- XVY + t(X.array[, , i]) %*% invcov %*% z.vec[id_r[i, ]]
  }
  XVY <- as.numeric(XVY)
  
  #beta�̕��z�̕��U�����U�s��p�����[�^
  inv_XVX <- solve(XVX + Adelta)
  
  #beta�̕��z�̕��σp�����[�^
  B <- inv_XVX %*% (XVY + Adelta %*% Deltabar)   #beta�̕���
  
  #���ϗʐ��K���z����beta���T���v�����O
  oldbeta <- mvrnorm(1, as.numeric(B), inv_XVX)
  
  #���݌��p�̕��ύ\�����X�V
  util.vec <- XM.vec %*% oldbeta
  old.util <- matrix(util.vec, nrow=hh, ncol=items, byrow=T)
  
  ##Cov�̕��z�̃p�����[�^�̌v�Z��mcmc�T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^���v�Z
  R.error <- matrix(z.vec - util.vec, nrow=hh, ncol=items, byrow=T)
  IW.R <- V + matrix(rowSums(apply(R.error, 1, function(x) x %*% t(x))), nrow=items, ncol=items, byrow=T)
  
  #�t�E�B�V���[�g���z�̎��R�x���v�Z
  Sn <- nu + hh
  
  #�t�E�B�V���[�g���z����Cov���T���v�����O
  Cov_hat <- rwishart(Sn, solve(IW.R))$IW
  
  
  ##���ʐ��̖���������邽�߂ɕ��U�����U�s��̑Ίp������1�ɌŒ肷��
  lambda <- diag(diag(Cov_hat)^(-1/2))
  oldcov <- cov2cor(Cov_hat)
  old.utilm <- old.util %*% lambda
  Z <- Z %*% lambda
  
  
  ##���p�֐�����w���ʂ��x�C�W�A���|�A�\����A���f���Ō��т���
  #�����_���E�H�[�N�T���v�����O
  if(rp %in% c(1, 100, 500, 1000, 10000)) {
    print("rw���X�V")
    res <- glm(y2 ~ -1 + ., data=data.frame(ZM), family=poisson)
    rw <- summary(res)[[12]][, 2]^2
  } 
  
  newtheta <- oldtheta + mvrnorm(1, rep(0, ncol(ZM)), diag(rw))   #�V����theta
  ZM[, 1] <- z.vec
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew <- loglike(theta=newtheta, y=y2, X=ZM)$LL
  logold <- loglike(theta=oldtheta, y=y2, X=ZM)$LL
  logpnew <- lndMvn(newtheta, Zetabar, Azeta)
  logpold <- lndMvn(oldtheta, Zetabar, Azeta)
  
  #MH�T���v�����O
  alpha <- min(1, exp(lognew + logpnew - logold - logpold))
  if(alpha == "NAN") alpha <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V�����Œ����beta���̑�
  if(u < alpha){
    oldtheta <- newtheta
    logl <- lognew
    
    #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
  } else {
    logl <- logold
  }
  
  ##�T���v�����O���ʂ�ۑ�
  mkeep <- rp/keep
  if(rp%%keep==0){
    keep_er <- mkeep
    Z_error <- Zold
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    THETA[mkeep, ] <- oldtheta
    print(rp)
    print(round(rbind(oldbeta, betat), 2))
    print(round(rbind(oldtheta, thetat), 3))
    print(round(rbind(MVR.S, sqrt(cd_Cov)), 2))
    print(round(cbind(oldcov[1:5, 1:5], Cov[1:5,1:5]), 2))
    print(round(logl, 3))
  }
}


####�T���v�����O���ʂ̊m�F�Ɨv��####
burnin <- 1000   #�o�[���C�����Ԃ�4000�T���v���܂�

##�T���v�����O���ʂ̃v���b�g
#��A�W���̃v���b�g
matplot(BETA[, 1:5], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 6:10], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 11:15], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 16:20], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 21:25], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 26:30], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 31:35], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")

#���U�����U�s��̉���
matplot(SIGMA[, 1:5], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 6:10], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 11:15], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 16:20], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 21:25], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 26:30], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 31:35], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 36:40], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 41:45], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 46:50], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 51:55], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 56:60], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 61:65], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 66:70], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 71:75], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 76:80], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 81:85], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 86:90], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA[, 91:95], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")

#�|�A�\����A�̃v���b�g
matplot(THETA[, 1:4], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")


##���茋�ʂ̗v��
#beta�̗v�񓝌v��
round(colMeans(BETA[burnin:nrow(BETA), ]), 2)   #beta�̎��㕽��
round(betat, 2)   #�^��beta
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5�����ʓ_
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95�����ʓ_
round(apply(BETA[burnin:nrow(BETA), ], 2, sd), 2)   #����W���΍�


#sigma�̗v�񓝌v��
round(colMeans(SIGMA[burnin:nrow(SIGMA), ]), 2)   #sigma�̎��㕽��
round(as.numeric(Cov), 2)   #�^��sigma
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.05)), 2)   #5�����ʓ_
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.95)), 2)   #95�����ʓ_
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, sd), 2)   #����W���΍�

#theta�̗v�񓝌v��
round(colMeans(THETA[burnin:nrow(THETA), ]), 2)   #beta�̎��㕽��
round(thetat, 2)   #�^��beta
round(apply(theta[burnin:nrow(THETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5�����ʓ_
round(apply(theta[burnin:nrow(THETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95�����ʓ_
round(apply(BETA[burnin:nrow(THETA), ], 2, sd), 2)   #����W���΍�