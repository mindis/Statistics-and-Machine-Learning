#####���݃N���X���f��#####
####�����������z���f��####
####�f�[�^�̔���####
##���f���̐ݒ�
k <- 4   #�Z�O�����g��
n <- 500   #�Z�O�����g�̃T���v����
N <- n * k   #���T���v����
ns <- 30   #�l���Ƃ̕p�x
th <- 15   #�Z�O�����g�ʂ̃p�����[�^��

##�m���x�N�g�����`
#�Z�O�����g1�̏o����
x1 <- rnorm(th, 0, 1)
a1 <- rnorm(th, 0, 0.7)
(p1 <- round(exp(a1*x1) / sum(exp(a1*x1)), 3))

#�Z�O�����g2�̏o����
x2 <- rnorm(th, 0, 1)
a2 <- rnorm(th, 0, 0.3)
(p2 <- round(exp(a2*x2) / sum(exp(a2*x2)), 3))

#�Z�O�����g3�̏o����
x3 <- rnorm(th, 0, 1)
a3 <- rnorm(th, 0, 1.0)
(p3 <- round(exp(a3*x3) / sum(exp(a3*x3)), 3))

#�Z�O�����g4�̏o����
x4 <- rnorm(th, 0, 1)
a4 <- rnorm(th, 0, 1.2)
(p4 <- round(exp(a4*x4) / sum(exp(a4*x4)), 3))

#�^�̊m��
thetaT <- rbind(p1, p2, p3, p4)

##�Z�O�����g�̊m���Ɋ�Â��p�x�f�[�^�𐶐�
#�Z�O�����g������̐l����500�l�A1�l������̕p�x��30��
#�Z�O�����g1�̃f�[�^�𔭐�
Y1 <- t(rmultinom(n, ns, p1))
dim(Y1)   #������

#�Z�O�����g2�̃f�[�^�𔭐�
Y2 <- t(rmultinom(n, ns, p2))
dim(Y2)   #������

#�Z�O�����g3�̃f�[�^�𔭐�
Y3 <- t(rmultinom(n, ns, p3))
dim(Y3)   #������

#�Z�O�����g4�̃f�[�^�𔭐�
Y4 <- t(rmultinom(n, ns, p4))
dim(Y4)   #������

#�Z�O�����g��\���x�N�g��
Seg <- rep(1:4, rep(n, 4))


##�f�[�^���������Č��ʂ��W�v
Y <- as.data.frame(rbind(Y1, Y2, Y3, Y4))   #�Z�O�����g�x�N�g���Ȃ��̃f�[�^�s��
Ys <- as.data.frame(cbind(Seg, Y))   #�Z�O�����g�x�N�g������̃f�[�^�s��
by(Ys[, 2:16], Ys[, 1], function(x) round(colMeans(x), 3))   #�Z�O�����g�ʂ̗񂲂Ƃ̔����p�x�̕���
by(Ys[, 2:16], Ys[, 1], function(x) summary(x))   #�Z�O�����g�ʂ̏W�v
by(Ys[, 2:16], Ys[, 1], function(x) round(colSums(x)/sum(x), 3))   #�Z�O�����g�ʂ̔�����
dim(Y)   #������

####EM�A���S���Y���Ő��݃N���X���f���𐄒�####
k <- 4   #�Z�O�����g��
n <- 500   #�Z�O�����g�̃T���v����
N <- n * k   #���T���v����
ns <- 30   #�l���Ƃ̕p�x
th <- 15   #�p�����[�^��

##�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z���v�Z���邽�߂̊֐�
LLobz <- function(theta, y, r, k){
  LLind <- matrix(0, nrow=nrow(y), ncol=k)
  for(i in 1:k){
    Li <- apply(y, 1, function(x) dmultinom(x, ns, theta[i, ]))   #�������z�̖ޓx���v�Z
    LLind[, i] <- Li
  }
  LLho <- matrix(r, nrow=nrow(y), ncol=k, byrow=T) * LLind   #�ϑ��f�[�^�̖ޓx
  z <- LLho / matrix(apply(LLho, 1, sum), nrow=nrow(y), ncol=k)   #���ݕϐ�z�̌v�Z
  LLosum <- sum(log(apply(matrix(r, nrow=nrow(y), ncol=k) * LLind, 1, sum)))   #�ϑ��f�[�^�̑ΐ��ޓx�̌v�Z
  rval <- list(LLob=LLosum, z=z, LL=LLind, Li=Li)
  return(rval)
}

#�����l�̐ݒ�
iter <- 0
k <- 4   #�Z�O�����g��

##theta�̏����l�̐ݒ�
#�Z�O�����g1�̏����l
minmax <- colSums(Y)
hh1 <- runif(th, min(minmax), max(minmax))
theta1 <- hh1/sum(hh1)

#�Z�O�����g2�̏����l
hh2 <- runif(th, min(minmax), max(minmax))
theta2 <- hh2/sum(hh2)

#�Z�O�����g3�̏����l
hh3 <- runif(th, min(minmax), max(minmax))
theta3 <- hh3/sum(hh3)

#�Z�O�����g4�̏����l
hh4 <- runif(th, min(minmax), max(minmax))
theta4 <- hh4/sum(hh4)

#theta�̏����l������
(theta <- rbind(theta1, theta2, theta3, theta4))

##������r�̏����l
r <- c(0.4, 0.3, 0.2, 0.2)

#�ΐ��ޓx�̏�����
L <- LLobz(theta=theta, y=Y, r=r, k=k)
LL1 <- L$LLob
z <- L$z
round(z, 3)

#�X�V�X�e�[�^�X
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 0.1  

##EM�A���S���Y��
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  #E�X�e�b�v�̌v�Z
  z <- L$z   #���ݕϐ�z�̏o��
  
  #M�X�e�b�v�̌v�Z�ƍœK��
  #theta�̐���
  theta <- matrix(0, nrow=k, ncol=th)
  for(j in 1:k){
    #���S�f�[�^�̑ΐ��ޓx����theta�̐���ʂ��v�Z
    thetaseg <- apply(matrix(z[, j], nrow=nrow(Y), ncol=th)*Y, 2, sum) / sum((z[, j])*matrix(ns, nrow=nrow(Y), ncol=1))
    theta[j, ] <- as.matrix(thetaseg)
  }
  #�������𐄒�
  r <- apply(L$z, 2, sum) / nrow(Y)
  
  #�ϑ��f�[�^�̑ΐ��ޓx���v�Z
  L <- LLobz(theta=theta, y=Y, r=r, k=k)
  LL <- L$LLob   #�ϑ��f�[�^�̑ΐ��ޓx
  iter <- iter+1   
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

####���݃N���X���f���̐��茋��####
round(theta, 3)   #theta�̐����
round(thetaT, 3)   #theta�̐^�̒l
round(r, 3)   #�������̐����
round(z, 3)   #�l�ʂ̃Z�O�����g�ւ̏����m��

L$LLob   #�ϑ��f�[�^�̑ΐ��ޓx
-2*(L$LLob) + k*nrow(theta)   #AIC
