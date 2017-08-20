#####���݃N���X�𗘗p�����f�[�^�Z��####
library(MASS)
library(bayesm)
library(MCMCpack)
library(mlogit)
library(flexmix)
library(psych)
library(gtools)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 3000
h <- 1000   #�����v���C���Ă��郆�[�U�[��
k <- 8   #�Z�O�����g��


####���ʂ̐����ϐ��̔���####
cont <- 2; bin <- 3; multi <- 4
X.cont <- matrix(rnorm(hh*cont), nrow=hh, ncol=cont)
X.bin <- matrix(0, nrow=hh, ncol=bin)
X.multi <- matrix(0, nrow=hh, ncol=multi)

#��l�����ϐ���ݒ�
for(i in 1:bin){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(hh, 1, p)
}

#���l�����ϐ���ݒ�
p <- runif(multi)
X.multi <- t(rmultinom(hh, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))] #�璷�ȕϐ��͍폜 


##�����ϐ����x�N�g���`���ɕύX
#�ؕЂ��x�N�g���`���ɕύX
bv <- c(1, rep(0, k))
iv <- matrix(bv, nrow=hh*length(bv), ncol=k, byrow=T)
IV <- subset(iv, rowSums(iv) > 0)
IV <- IV[, -k]

#�����ϐ����x�N�g���`���ɕύX
index.z <- rep(1:hh, rep(k, hh))
Zi <- matrix(0, nrow=hh*k, ncol=(cont+bin+multi-1)*(k-1))
Xi <- cbind(X.cont, X.bin, X.multi)

for(i in 1:hh){
  x.bind <- c()
  for(j in 1:ncol(Xi)){
    x.diag <- diag(Xi[i, j], k)
    x.bind <- cbind(x.bind, x.diag[, -k])
  }
  Zi[index.z==i, ] <- x.bind
}

#�f�[�^������
Zx1 <- cbind(inter=IV, Z=Zi)

####�����p�^�[���̔���####
#�������邩�ǂ����̓����_���Ȍ���������
#�����p�^�[����3���
pattern <- 3
alpha0 <- runif(pattern-1, -0.6, 0.6)
alpha1 <- runif((pattern-1)*cont, 0, 0.7)
alpha2 <- runif((pattern-1)*bin, -0.8, 0.8)
alpha3 <- runif((pattern-1)*(multi-1), -0.9, 0.9)
alphat <- c(alpha0, alpha1, alpha2, alpha3)

index.pattern <- as.numeric(t(matrix(1:ncol(Zx1), ncol=k-1, byrow=T)[, 1:(pattern-1)]))

##���W�b�g�Ɗm���̌v�Z
logit1 <- matrix(Zx1[, index.pattern] %*% alphat, nrow=hh, ncol=pattern, byrow=T)
Pr1 <- exp(logit1)/matrix(rowSums(exp(logit1)), nrow=hh, ncol=pattern)

##�����ϐ��̔���
Z1 <- t(apply(Pr1, 1, function(x) rmultinom(1, 1, x)))
z1 <- Z1 %*% 1:pattern
z.vec <- z1
n_seg <- colSums(Z1)

colMeans(Z1); colSums(Z1)   #���ʂ��W�v


##�����p�^�[�����x�N�g���`���ɕύX
Z1.vec <- matrix(0, nrow=nrow(Zx1), ncol=(k-1)*2)

for(i in 1:hh){
  r <- ((i-1)*k+1):((i-1)*k+k)
  Z1.vec[r, ] <- cbind(diag(Z1[i, 1], k)[, -k], diag(Z1[i, 2], k)[, -k])
}

Zx2 <- cbind(Zx1, type=Z1.vec)   #�f�[�^������


####���݃N���X�𔭐�####
##�p�����[�^�̐ݒ�
theta0 <- runif(k-1, -0.55, 0.55)
theta1 <- runif((k-1)*cont, 0, 0.6)
theta2 <- runif((k-1)*bin, -0.7, 0.7)
theta3 <- runif((k-1)*(multi-1), -0.7, 0.7)
theta4 <- runif((k-1)*(pattern-1), -1.0, 1.0)
thetat <- c(theta0, theta1, theta2, theta3, theta4)

##���W�b�g�Ɗm���̌v�Z
logit2 <- matrix(Zx2 %*% thetat, nrow=hh, ncol=k, byrow=T)
Pr2 <- exp(logit2)/matrix(rowSums(exp(logit2)), nrow=hh, ncol=k)

##�����ϐ��̔���
Z2 <- t(apply(Pr2, 1, function(x) rmultinom(1, 1, x)))
z2 <- Z2 %*% 1:k
colMeans(Z2); colSums(Z2)   #���ʂ��W�v
mix_rate <- as.numeric(table(z2)/sum(table(z2))) 

#�����p�^�[���Ɛ��݃N���X�̃N���X�[
table(z1, z2)
round(table(z1, z2)/matrix(rowSums(table(z1, z2)), nrow=pattern, ncol=k), 3)


####�ϑ��ϐ��̔���####
#�����p�^�[����1�̓��u���C�u!�̂�2�̓A�C�}�X�̂�3�͗����ϑ�
k1 <- 9   #���u���C�u�f�[�^�̕ϐ�   
k2 <- 10   #�A�C�}�X�f�[�^�̕ϐ�

#�p�x�𔭐�������
ns11 <- rpois(hh, 40)
ns12 <- rpois(hh, 40)
ns21 <- rpois(hh, 30)
ns22 <- rpois(hh, 30)

##�m���x�N�g�����`���āA�����ϐ��𔭐�
P11 <- matrix(0, nrow=k, ncol=k1) 
P12 <- matrix(0, nrow=k, ncol=k2)
P21 <- rep(0, k)
P22 <- rep(0, k)

X11 <- matrix(0, nrow=hh, ncol=k1)
X12 <- matrix(0, nrow=hh, ncol=k2)
X21 <- rep(0, hh)
X22 <- rep(0, hh)

##���݃N���X���Ƃɉ����ϐ��𔭐�
for(j in 1:k){
  r <- subset(1:length(z2), z2==j)

  #�m�����v�Z
  p11 <- rgamma(k1, 0.85, 0.2)
  p12 <- rgamma(k2, 1.0, 0.15)
  P11[j, ] <- p11 / sum(p11)
  P12[j, ] <- p12 / sum(p12)
  
  P21[j] <- runif(1, 0.15, 0.7)
  P22[j] <- runif(1, 0.2, 0.6)
  
  #�����ϐ��𔭐�
  X11[r, ] <- t(apply(cbind(ns11[r], 1), 1, function(x) rmultinom(1, x[1], P11[j, ])))
  X12[r, ] <- t(apply(cbind(ns12[r], 1), 1, function(x) rmultinom(1, x[1], P12[j, ])))
  X21[r] <- apply(cbind(ns21[r], 1), 1, function(x) rbinom(1, x[1], P21[j]))
  X22[r] <- apply(cbind(ns22[r], 1), 1, function(x) rbinom(1, x[1], P22[j]))
}

##�����ϐ��̏W�v
by(X11, as.character(z2), colMeans)
by(X12, as.character(z2), colMeans)

##z1=3�ȊO�̃��[�U�[�̃f�[�^������������
X11[z1==2, ] <- 10
X12[z1==1, ] <- 10
X21[z1==2] <- 10
X22[z1==1] <- 10
ns11[z1==2] <- 10*k1
ns12[z1==1] <- 10*k2
ns21[z1==2] <- 20
ns22[z1==1] <- 20

X <- data.frame(ll1=X11, im1=X12, ll2=X21, im2=X22)
XM <- as.matrix(X)


####EM�A���S���Y���Ō����f�[�^�̃Z�O�����g�����m����\��####
####EM�A���S���Y���ɕK�v�ȑΐ��ޓx�֐����`####
##�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z���v�Z���邽�߂̊֐�
LLobz <- function(theta11, theta12, theta21, theta22, ns11, ns12, ns21, ns22, X11, X12, X21, X22, z, r, k, n_seg){
  
  #�����p�^�[�����Ƃɖޓx���v�Z
  LLind1 <- matrix(0, nrow=n_seg[1], ncol=k)
  LLind2 <- matrix(0, nrow=n_seg[2], ncol=k)
  LLind3 <- matrix(0, nrow=n_seg[3], ncol=k)
  
  #���݃N���X���Ƃɖޓx���v�Z���đ��
  for(i in 1:k){
    Li11 <- apply(cbind(ns11, X11), 1, function(x) dmultinom(x[-1], x[1], theta11[i, ]))   #���u���C�u�̑������z�̖ޓx
    Li12 <- apply(cbind(ns12, X12), 1, function(x) dmultinom(x[-1], x[1], theta12[i, ]))   #�A�C�}�X�̑������z�̖ޓx
    Li21 <- apply(cbind(ns21, X21), 1, function(x) dbinom(x[-1], x[1], theta21[i]))   #���u���C�u�̓񍀕��z�̖ޓx
    Li22 <- apply(cbind(ns22, X22), 1, function(x) dbinom(x[-1], x[1], theta22[i]))   #�A�C�}�X�̓񍀕��z�̖ޓx
    
    #�ޓx���v�Z
    Li1 <- Li11[z==1] * Li21[z==1]   #�A�C�}�X�f�[�^�������̖ޓx
    Li2 <- Li12[z==2] * Li22[z==2]   #���u���C�u�I�f�[�^�������̖ޓx
    Li3 <- Li11[z==3] * Li12[z==3] * Li21[z==3] * Li22[z==3]   #�����ϑ��f�[�^�̖ޓx
    LLind1[, i] <- Li1; LLind2[, i] <- Li2; LLind3[, i] <- Li3
  }
  
  #�ϑ��f�[�^�̖ޓx���v�Z
  LLho1 <- matrix(r[z==1, ], nrow=n_seg[1], ncol=k, byrow=T) * LLind1   #���u���C�u�݂̂̊ϑ��f�[�^�̖ޓx
  LLho2 <- matrix(r[z==2, ], nrow=n_seg[2], ncol=k, byrow=T) * LLind2   #�A�C�}�X�݂̂̊ϑ��f�[�^�̖ޓx
  LLho3 <- matrix(r[z==3, ], nrow=n_seg[3], ncol=k, byrow=T) * LLind3   #�����ϑ��̊ϑ��f�[�^�̖ޓx
  
  #���ݕϐ�z�̌v�Z
  z1 <- LLho1 / matrix(apply(LLho1, 1, sum), nrow=n_seg[1], ncol=k)  
  z2 <- LLho2 / matrix(apply(LLho2, 1, sum), nrow=n_seg[2], ncol=k) 
  z3 <- LLho3 / matrix(apply(LLho3, 1, sum), nrow=n_seg[3], ncol=k) 
  
  #�S�f�[�^�s��ɐ��ݕϐ�z����
  Z <- matrix(0, nrow=sum(n_seg), ncol=k)
  Z[z==1, ] <- z1
  Z[z==2, ] <- z2
  Z[z==3, ] <- z3
  
  #�ϑ��f�[�^�̑ΐ��ޓx�̘a���v�Z
  LLosum1 <- sum(log(apply(r[z==1, ] * LLind1, 1, sum)))   
  LLosum2 <- sum(log(apply(r[z==2, ] * LLind2, 1, sum)))  
  LLosum3 <- sum(log(apply(r[z==3, ] * LLind3, 1, sum)))  
  LLosum <- LLosum1 + LLosum2 + LLosum3

  rval <- list(LLob=LLosum, Z=Z)
  return(rval)
}

##�������W�b�g���f���̑ΐ��ޓx�֐�
loglike <- function(x, Zx, z, hh, seg){
  #�p�����[�^�̐ݒ�
  theta.z <- x
  
  #���p�֐��̐ݒ�
  U <- matrix(Zx %*% theta.z, nrow=hh, ncol=seg, byrow=T)
  
  #�ΐ��ޓx�̌v�Z
  d <- rowSums(exp(U))
  LLl <- rowSums(z * U) - log(d)
  LL <- sum(LLl)
  return(LL)
}

####EM�A���S���Y���̐ݒ�Ə����l�̐ݒ�####
##EM�A���S���Y���̐ݒ�
iter <- 0
dl <- 100
tol <- 0.5
maxit <- c(10, 20) #���j���[�g���@�̃X�e�b�v��


##�f�[�^�̏���
index.z13 <- subset(1:hh, z1 %in% c(1,3))
index.z23 <- subset(1:hh, z1 %in% c(2,3))
X11_z <- X11[index.z13, ]
X12_z <- X12[index.z23, ]
X21_z <- X21[index.z13]
X22_z <- X22[index.z23]
ns11_z <- ns11[index.z13]
ns12_z <- ns12[index.z23]
ns21_z <- ns21[index.z13]
ns22_z <- ns22[index.z23]


##�����l�̐ݒ�
#�x�X�g�ȏ����p�����[�^��I��
rp <- 20   #�J��Ԃ���
Z.list <- list()
val <- c()
x <- list()
theta.x <- matrix(0, nrow=rp, ncol=ncol(Zx2))

#�p�����[�^�̊i�[�p�z��
beta11 <- matrix(0, nrow=k, ncol=k1)
beta12 <- matrix(0, nrow=k, ncol=k2)
beta21 <- rep(0, k)
beta22 <- rep(0, k)

#�x�X�g�ȏ����p�����[�^���I�������܂Ŕ��������
for(m in 1:rp){
  for(j in 1:k){
    p11 <- abs(colSums(X11, na.rm=TRUE)/sum(X11, na.rm=TRUE) + runif(k1, -0.1, 0.15))
    p12 <- abs(colSums(X12, na.rm=TRUE)/sum(X12, na.rm=TRUE) + runif(k2, -0.1, 0.15))
    beta11[j, ] <- p11 / sum(p11)
    beta12[j, ] <- p12 / sum(p12)
    
    beta21[j] <- runif(1, 0.2, 0.6)
    beta22[j] <- runif(1, 0.2, 0.6)
  }
  beta <- list(beta11=beta11, beta12=beta12, beta21=beta21, beta22=beta22)
  x[[m]] <- beta 
  
  #���݃N���X�����̃p�����[�^�̏����l
  theta0 <- runif(k-1, -0.55, 0.55)
  theta1 <- runif((k-1)*cont, 0, 0.6)
  theta2 <- runif((k-1)*bin, -0.7, 0.7)
  theta3 <- runif((k-1)*(multi-1), -0.7, 0.7)
  theta4 <- runif((k-1)*(pattern-1), -1.0, 1.0)
  theta <- c(theta0, theta1, theta2, theta3, theta4)
  theta.x[m, ] <- theta
  
  #���O�m���̏����l
  logit <- matrix(Zx2 %*% theta, nrow=hh, ncol=k, byrow=T)
  r <- exp(logit)/matrix(rowSums(exp(logit)), nrow=hh, ncol=k)
  
  #�ϑ��f�[�^�̑ΐ��ޓx���v�Z
  oll <- LLobz(beta11, beta12, beta21, beta22, ns11, ns12, ns21, ns22, X11, X12, X21, X22, z.vec, r, k, n_seg)
  val <- c(val, oll$LLob)
  Z.list[[m]] <- oll$Z
  print(m)
}

##�x�X�g�ȃp�����[�^��I��
opt <- which.max(val)   #�x�X�g�ȑΐ��ޓx
LL1 <- val[opt]
Z <- Z.list[[opt]]
beta <- x[[opt]]
theta <- theta.x[opt, ]
beta11 <- matrix(0, nrow=k, ncol=k1)
beta12 <- matrix(0, nrow=k, ncol=k2)
beta21 <- rep(0, k)
beta22 <- rep(0, k)


####EM�A���S���Y���ɂ��f�[�^�Z��####
while(abs(dl) >= tol){   #dl��tol�ȏ�Ȃ�J��Ԃ�

  ##���S�f�[�^�ł̃p�����[�^�̐���
  #���݃N���X���Ƃɏd�ݕt���ΐ��ޓx���ő剻����
  for(j in 1:k){
    beta11[j, ] <- colSums((matrix(Z[index.z13, j], nrow=length(index.z13), ncol=k1) * X11_z) / sum(Z[index.z13, j] * ns11_z))
    beta12[j, ] <- colSums((matrix(Z[index.z23, j], nrow=length(index.z23), ncol=k2) * X12_z) / sum(Z[index.z23, j] * ns12_z))
    beta21[j] <- sum(Z[index.z13, j] * X21_z) / sum(Z[index.z13, j] * ns21_z)
    beta22[j] <- sum(Z[index.z23, j] * X22_z) / sum(Z[index.z23, j] * ns22_z)
  }
  
  ##�������W�b�g���f���Ōl���Ƃɍ������𐄒�
  #���j���[�g���@�ōŖސ���
  res <- optim(theta, loglike, gr=NULL, Zx=Zx2, z=Z, hh=hh, seg=k, method="BFGS", hessian=FALSE, 
               control=list(fnscale=-1))
  theta <- res$par   #�p�����[�^���X�V
  
  #���������X�V
  U <- matrix(Zx2 %*% theta, nrow=hh, ncol=k, byrow=T)
  r <- exp(U) / matrix(rowSums(exp(U)), nrow=hh, ncol=k)
  
  ##E�X�e�b�v�ł̑ΐ��ޓx�̊��Ғl
  #�ϑ��f�[�^�ł̑ΐ��ޓx�̊��Ғl���v�Z
  obsllz <- LLobz(beta11, beta12, beta21, beta22, ns11, ns12, ns21, ns22, X11, X12, X21, X22, z.vec, r, k, n_seg)
  LL <- obsllz$LLob
  Z <- obsllz$Z
  
  ##EM�A���S���Y���̃p�����[�^�̍X�V
  iter <- iter+1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
  print(round(rbind(mix_rate=colSums(Z)/hh, mix_rate), 3))
}

####���茋�ʂƗv��####
##�^�̃p�����[�^�őΐ��ޓx���v�Z
obsll_true <- LLobz(P11, P12, P21, P22, ns11, ns12, ns21, ns22, X11, X12, X21, X22, z.vec, Pr2, k, n_seg)
LL_true <- obsll_true$LLob

##���肳�ꂽ�p�����[�^
#beta�̐���l�Ɛ^�̒l
round(rbind(beta11, P11), 3)
round(rbind(beta12, P12), 3)
round(rbind(beta21, P21), 3)
round(rbind(beta22, P22), 3)

#theta�̐���l
round(theta, 3)
#round(tval <- theta/sqrt(-diag(solve(res$hessian))), 3)   #t�l

#���݃N���X�ƍ������̐���l
print(round(rbind(mix_rate=colSums(Z)/hh, mix_rate), 3))   #�������̔�r
round(cbind(Z, z2), 3)   #���݃N���X�̊����m���̐���l�Ɛ^�̃N���X
Z_seg <- apply(Z, 1, which.max)   #���݃N���X�̊���


##�K���x���v�Z
par <- length(beta11) + length(beta12) + length(beta21) + length(beta22)
round(c(LL, LL_true), 3)   #�ő剻���ꂽ�ΐ��ޓx
(AIC <- -2*LL + 2*par)   #AIC
(BIC <- -2*LL + log(hh)*par) #BIC


####���肳�ꂽ���݃N���X����уp�����[�^����u�[�g�X�g���b�v�@�Ō����l����####
boot <- 500   #�u�[�g�X�g���b�v��

#�����f�[�^�̃C���f�b�N�X���쐬
index.z1 <- subset(1:hh, z1==2)
index.z2 <- subset(1:hh, z1==1)

#�f�[�^�̊i�[�p�z��
X11_samp <- X11
X12_samp <- X12
X21_samp <- X21
X22_samp <- X22
X11_boot <- array(0, dim=c(hh, k1, boot))
X12_boot <- array(0, dim=c(hh, k2, boot))
X21_boot <- matrix(0, nrow=hh, ncol=boot)
X22_boot <- matrix(0, nrow=hh, ncol=boot)
seg1 <- matrix(0, nrow=length(index.z1), ncol=boot)
seg2 <- matrix(0, nrow=length(index.z2), ncol=boot)

  
##�u�[�g�X�g���b�s���O�Ō����f�[�^�����T���v�����O
for(b in 1:boot){
  print(b)
  
  #�Z�O�����g����
  seg1[, b] <- as.numeric(t(apply(Z[index.z1, ], 1, function(x) rmultinom(1, 1, x))) %*% 1:k)
  seg2[, b] <- as.numeric(t(apply(Z[index.z2, ], 1, function(x) rmultinom(1, 1, x))) %*% 1:k)
  
  #�Z�O�����g��������p�����[�^�Ɋ�Â��ă��T���v�����O
  X11_samp[index.z1, ] <- t(apply(cbind(rpois(length(ns11[index.z1]), 40), beta11[seg1[, b], ]), 1, function(x) rmultinom(1, x[1], x[-1])))
  X12_samp[index.z2, ] <- t(apply(cbind(rpois(length(ns12[index.z2]), 40), beta12[seg2[, b], ]), 1, function(x) rmultinom(1, x[1], x[-1])))
  X21_samp[index.z1] <- apply(cbind(rpois(length(ns21[index.z1]), 30), beta21[seg1[, b]]), 1, function(x) rbinom(1, x[1], x[-1]))
  X22_samp[index.z2] <- apply(cbind(rpois(length(ns22[index.z2]), 30), beta22[seg2[, b]]), 1, function(x) rbinom(1, x[1], x[-1]))
  
  #���T���v�����O���ʂ����X�g�ɕۑ�
  X11_boot[, , b] <- X11_samp
  X12_boot[, , b] <- X12_samp
  X21_boot[, b] <- X21_samp
  X22_boot[, b] <- X22_samp
}

