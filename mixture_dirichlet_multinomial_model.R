#####����-�f�B�N�����������z���f��#####
library(MASS)
library(vcd)
library(gtools)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
##���f���̐ݒ�
k <- 4   #�Z�O�����g��
hh <- 500   #�Z�O�����g�̃T���v����
n <- hh * k   #���T���v����
pt <- rpois(n, 10)   #�T���v�����Ƃ̃Z�b�V������
pt <- ifelse(pt < 2, 2, pt)
N <- sum(pt)   #���T���v����

th <- 20   #�Z�O�����g�ʂ̃p�����[�^��
freq <- rpois(N, 25)   #�l���Ƃ̂̕p�x

##ID�̐ݒ�
id <- rep(1:n, pt)   #id�̐ݒ�
t <- c()   #�Z�b�V�������̐ݒ�
for(i in 1:n){
  t <- c(t, 1:pt[i])
}  
no <- 1:N   #�����ԍ���ݒ�
ID <- data.frame(no, id, t)   #�f�[�^������

#�Z�O�����g�̐ݒ�
seg.z <- rep(1:4, rep(hh, k))

##�Z�O�����g���ƂɊm���x�N�g�����`
#�m���x�N�g���̌W���̐ݒ�
x <- matrix(rnorm(th*k, 0, 1), nrow=k, ncol=th, byrow=T)
a.vec <- matrix(0, nrow=k, ncol=th)
a <- c(0.2, 0.5, 0.9, 1.2)
for(i in 1:k){
  a.vec[i, ] <- rnorm(th, 0, a[i]) 
}

#�m���x�N�g���̌v�Z
Pr <- matrix(0, nrow=k, ncol=th)
for(i in 1:k){
  Pr[i, ] <- exp(x[i, ]*a.vec[i, ]) / sum(exp(x[i, ]*a.vec[i, ]))
}
#�����������m�����m�F
round(Pr, 3); apply(Pr, 1, summary)


##�����������m���ɂ��ƂÂ��p�x�f�[�^�𐶐�
#�p�x�f�[�^�𔭐�
Y <- matrix(0, nrow=N, ncol=th)
for(i in 1:N){
  Y[i, ] <- t(rmultinom(1, freq[i], Pr[seg.z[ID[i, 2]], ]))
}

#�o�������v�Z
YP <- matrix(0, nrow=k, ncol=th)  
for(i in 1:k){
  index.z <- subset(1:length(seg.z), seg.z==i)
  YP[i, ] <- colSums(Y[ID[, 2] %in% index.z, ]) / sum(Y[ID[, 2] %in% index.z, ])
}
round(YP, 3)
round(Pr, 3)   #�Z�O�����g�ʂ̊m��
round(colSums(Y)/sum(Y), 3)   #�S�̂ł̏o����


####EM�A���S���Y���ō�������-�f�B�N�������f���𐄒�####
##�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z���v�Z���邽�߂̊֐�
LLobz <- function(theta, r, Y, ID, n, k, freq, v){
  #�ޓx�Ƒΐ��ޓx���v�Z
  LLind <- matrix(0, nrow=N, ncol=k)
  for(i in 1:k){
    Li <- apply(cbind(Y, freq), 1, function(x) dmultinom(x[1:v], x[v+1], theta[i, ]))
    LLind[, i] <- Li 
  }
  
  #�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̌v�Z
  #������
  R <- matrix(r, nrow=n, ncol=k)
  
  #�l�ʂ̐��݊m���̌v�Z
  LLd <- matrix(0, nrow=n, ncol=k)
  LLl <- log(LLind)
  for(i in 1:n){
    LLd[i, ] <- apply(LLl[ID[, 2]==i, ], 2, sum)
  }
  
  #��������h�����߂ɑΐ��ޓx��-744�ȉ��̏ꍇ�͑ΐ��ޓx�𐓏グ����
  LL.min <- apply(LLd, 1, min)
  index.loss <- subset(1:nrow(LLd), (LL.min + 743) < 0)
  lplus <- -matrix((LL.min[index.loss] + 743), nrow=length(index.loss), ncol=k)
  LLd[index.loss, ] <- LLd[index.loss, ] + lplus
  
  #���݊m��z�̌v�Z
  LLho <- R * exp(LLd)
  z <- LLho / matrix(rowSums(LLho), nrow=n, ncol=k)
  
  #�ϑ��f�[�^�̑ΐ��ޓx���v�Z
  LLosum <- sum(log(apply(matrix(r, nrow=n, ncol=k, byrow=T) * exp(LLd), 1, sum)))
  rval <- list(LLobz=LLosum, z=z, LL=LLd)
  return(rval)
}

##EM�A���S���Y���̐ݒ�
#�X�V�X�e�[�^�X
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 1
iter <- 0

#�����l�̐ݒ�
alpha <- rep(2, th)   #���O���z�̃p�����[�^
r <- c(0.4, 0.3, 0.2, 0.2)   #�������̏����l

#�p�����[�^�̏����l
theta.f <- matrix(0, nrow=k, ncol=th)
for(i in 1:k){
  minmax <- colSums(Y)
  pf <- runif(th, min(minmax), max(minmax))
  theta.f[i, ] <- pf/sum(pf)
}

#�ΐ��ޓx�̏�����
L <- LLobz(theta=theta.f, r=r, Y=Y, ID=ID, n=n, k=k, freq=freq, v=th)
LL1 <- L$LLob
z <- L$z

##EM�A���S���Y��
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  #E�X�e�b�v�̌v�Z
  z <- L$z   #���ݕϐ�z�̏o��
  zpt <- matrix(0, nrow=N, ncol=k)
  for(i in 1:n){
    zpt[ID[, 2]==i, ] <- matrix(z[i, ], nrow=length(ID[ID[, 2]==i, 2]), ncol=k, byrow=T)
  }
  
  #M�X�e�b�v�̌v�Z�ƍœK��
  #theta�̐���
  theta <- matrix(0, nrow=k, ncol=th)
  for(j in 1:k){
    #���S�f�[�^�̑ΐ��ޓx����theta�̐���ʂ��v�Z
    theta.seg <- (colSums(zpt[, j]*Y) + r[j]*(alpha-1)) / (sum(zpt[, j]*Y) + r[j]*sum(alpha-1))
    theta[j, ] <- as.matrix(theta.seg)
  }
  #�������𐄒�
  r <- apply(z, 2, sum) / n
  
  #�ϑ��f�[�^�̑ΐ��ޓx���v�Z
  L <- LLobz(theta=theta, r=r, Y=Y, ID=ID, n=n, k=k, freq=freq, v=th)
  LL <- L$LLob   #�ϑ��f�[�^�̑ΐ��ޓx
  iter <- iter+1   
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

####���茋�ʂƗv��
##���肳�ꂽ�p�����[�^
round(theta, 3)   #���肳�ꂽ�p�����[�^
round(Pr, 3)   #�^�̃p�����[�^
round(r, 3)   #������
round(data.frame(seg=apply(z, 1, which.max), z=z), 3)   #�l�ʂ̃Z�O�����g�ւ̏����m���Ə����Z�O�����g

#���肳�ꂽ�p�����[�^�̃O���t
max.theta <- max(theta)
par(mfrow=c(2,2)) 
barplot(theta[1, ], ylim=c(0, max.theta))
barplot(theta[2, ], ylim=c(0, max.theta))
barplot(theta[3, ], ylim=c(0, max.theta))
barplot(theta[4, ], ylim=c(0, max.theta))
par(mfrow=c(1, 1))

##�K���x
round(LL, 3)   #�ő剻���ꂽ�ΐ��ޓx
round(AIC <- -2*LL + 2*(length(th*k)+k), 3)   #AIC