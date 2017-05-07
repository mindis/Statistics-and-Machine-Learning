#####�������K���ϗʉ�A���f���̐���#####
library(MASS)
library(plyr)
library(lavaan)
####�f�[�^�̔���####
#set.seed(4543)
k <- 4   #������
col <- 5   #�ϐ���
n <- 4000   #�T���v����

##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
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

#�������z���Ƃ̑��֍s����쐬
corM1 <- corrM(col=5, lower=-0, upper=0)
corM2 <- corrM(col=5, lower=-0, upper=0)
corM3 <- corrM(col=5, lower=-0, upper=0)
corM4 <- corrM(col=5, lower=-0, upper=0)

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

#���U�����U�s����쐬
Sigma1 <- covmatrix(col=5, corM=corM1, lower=5, upper=12)
Sigma2 <- covmatrix(col=5, corM=corM2, lower=8, upper=14)
Sigma3 <- covmatrix(col=5, corM=corM3, lower=9, upper=16)
Sigma4 <- covmatrix(col=5, corM=corM4, lower=4, upper=18)

Sigma1[-2]
Sigma2[-2]
Sigma3[-2]
Sigma4[-2]

##���ϗʐ��K���z�̕��ς���A���f���ŕ\������
#���͂���f�[�^�̔���
#����
sex1 <- rbinom(1000, 1, 0.7)
sex2 <- rbinom(1000, 1, 0.5)
sex3 <- rbinom(1000, 1, 0.4)
sex4 <- rbinom(1000, 1, 0.3)
sex <- c(sex1, sex2, sex3, sex4)

#�N��
age1 <- t(rmultinom(1000, 1, c(0.2, 0.3, 0.2, 0.2, 0.1)))
age2 <- t(rmultinom(1000, 1, c(0.1, 0.4, 0.3, 0.1, 0.1)))
age3 <- t(rmultinom(1000, 1, c(0.1, 0.2, 0.2, 0.3, 0.2)))
age4 <- t(rmultinom(1000, 1, c(0.4, 0.3, 0.1, 0.05, 0.05)))
age <- rbind(age1, age2, age3, age4)

#�E��
job1 <- t(rmultinom(1000, 1, c(0.5, 0.2, 0.2, 0.1)))
job2 <- t(rmultinom(1000, 1, c(0.6, 0.4, 0.05, 0.05)))
job3 <- t(rmultinom(1000, 1, c(0.4, 0.4, 0.1, 0.1)))
job4 <- t(rmultinom(1000, 1, c(0.4, 0.3, 0.1, 0.2)))
job <- rbind(job1, job2, job3, job4)

#�ݐϗ��X��
cnt1 <- rpois(1000, 9)
cnt2 <- rpois(1000, 6)
cnt3 <- rpois(1000, 12)
cnt4 <- rpois(1000, 4)
cnt <- c(cnt1, cnt2, cnt3, cnt4)

#�Z�O�����g�ԍ�
segment <- rep(1:4, c(rep(1000, 4)))

#�f�[�^������
data <- as.data.frame(cbind(segment, sex, age, job, cnt))
head(data, 10); tail(data, 10)

#�ϐ��ɖ��O������
names(data)[3:12] <- c("age20", "age30", "age40", "age50", "age60u", "work", "stud", "hwife", "others", "s_cnt")
head(data, 10); tail(data, 10)
summary(data)


#��ϐ����폜����
data1 <- data[, -7]   #60��ȏ���폜
data_n <- data1[, -10]   #���̑��̐E�Ƃ��폜
data_de <- cbind(1, data_n[, -1])
names(data_de)[1] <- ("intercept")

#�Z�O�����g���Ƃ̉�A�W��������
v <- dim(data_de)[2]-1
a1 <- matrix(c(rnorm(v*col, 5.4, 2.5), runif(5, 1.0, 1.6)), nrow=v+1, ncol=col, byrow=T) 
a2 <- matrix(c(rnorm(v*col, 3.2, 2.2), runif(5, 0.8, 1.8)), nrow=v+1, ncol=col, byrow=T) 
a3 <- matrix(c(rnorm(v*col, 1.9, 1.9), runif(5, 0.2, 1.2)), nrow=v+1, ncol=col, byrow=T) 
a4 <- matrix(c(rnorm(v*col, 0.7, 1.4), runif(5, 0.4, 1.0)), nrow=v+1, ncol=col, byrow=T) 

a1; a2; a3; a4   #�ϐ����m�F

#��A�W�������X�g�\���ɂ���
a <- list(a1, a2, a3, a4)

#���ύ\���𑽕ϗʉ�A���f���Ŕ���������
y <- matrix(0, 0, 11)
for(i in 1:k){
  seg_x <- subset(data_de, data_n[, 1] == i)
  y_temp <- as.matrix(seg_x) %*% as.matrix(a[[i]])
  y_temp <- as.data.frame(y_temp)
  y <- rbind(y, y_temp)
}
round(y, 3)
names(y) <- c("y1", "y2", "y3", "y4", "y5")
y_seg <- cbind(segment, y)
by(y_seg[, 2:6], y_seg$segment, colMeans)   #�Z�O�����g���Ƃ̕���

#�Z�O�����g���Ƃ�y�𒊏o
y1 <- subset(y_seg[, 2:6], y_seg[, 1] == 1)
y2 <- subset(y_seg[, 2:6], y_seg[, 1] == 2)
y3 <- subset(y_seg[, 2:6], y_seg[, 1] == 3)
y4 <- subset(y_seg[, 2:6], y_seg[, 1] == 4)

n/col
dim(y_seg)

##���ϗʐ��K���z�����A�\���̂��鍬�����K�����𔭐�������
k = 4
n = 4000
y <- matrix(0, nrow=n, ncol=col+1)
cov <- list(Sigma1$covariance, Sigma2$covariance,Sigma3$covariance, Sigma4$covariance)
for(k in 1:4){
  yy <- subset(y_seg[, 2:6], y_seg[, 1] == k)
  cc <- as.matrix(cov[[k]])
  for(i in 1:1000){
    r <- (k-1)*1000 + i
    y[r, ] <- c(k, mvrnorm(n=1, as.matrix(yy[i, ]), cc))
  }
}
#�f�[�^������
round(y, 3)
yy <- as.data.frame(y)
by(y_seg[, 2:6], y_seg[, 1], colMeans)   #���̃f�[�^�̃Z�O�����g���Ƃ̕���
by(yy[, 2:6], yy[, 1], colMeans)   #���������������f�[�^�̃Z�O�����g���Ƃ̕���
by(yy[, 2:6], yy[, 1], cor)   #���������������f�[�^�̃Z�O�����g���Ƃ̑���

#���ʂ��v���b�g
boxplot(yy[, 2:6])   #�Z�O�����g�𖳎�����Ɓc

#�Z�O�����g�ʂɔ��Ђ��}��`��
par(mfrow=c(2, 3))
for(i in 2:6){
  boxplot(yy[, i] ~ yy[, 1], data=yy)
}
par(mfrow=c(1, 1))

plot(yy[, 2:6], col=yy[, 1])   #�U�z�}

#�A�C���X�f�[�^�̔��Ђ��}
#par(mfrow=c(2, 2))
#for (response in c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"))
#     boxplot(iris[, response] ~ Species, data=iris, ylab=response)
#par(mfrow=c(1, 1))

#####�Z�O�����g���Ƃɑ��ϗʉ�A���f���𓖂Ă͂߂�####
##��A���f���𓖂Ă͂߂鑼�̃f�[�^�̏���
#�����ϐ�
Yn <- y[, -1]   
Yn[1:10, ]; Yn[1001:1010, ]; Yn[2001:2010, ]; Yn[3001:3010, ]

iris
#�����ϐ�
head(data_de)   #�ؕЂ������f�U�C���s��
head(segment); tail(segment)   #�Z�O�����g
segment <- as.data.frame(segment)
data_seg <- cbind(segment, data_de)   #�ؕЂƃZ�O�����g�������f�U�C���s��

#���ׂĂ������s��
yx <- cbind(segment, Yn, data_de)
names(yx)[2:6] <- c("y1", "y2", "y3", "y4", "y5")
head(yx)


##�Z�O�����g���Ƃ̑��ϗʉ�A���f�������s
beta_seg <- list()   #��A�W�������郊�X�g
S_seg <- list()   #���U�����U�s������郊�X�g
for(i in 1:4){
  yx_s <- subset(yx, yx$segment==i)
  #��A�W���𐄒�
  beta_seg[[i]] <- solve(t(yx_s[, 7:16]) %*% as.matrix(yx_s[, 7:16])) %*% t(yx_s[, 7:16]) %*% as.matrix(yx_s[, 2:6])
  bs <- as.matrix(beta_seg[[i]])
  #���U�����U�s��𐄒�
  S_seg[[i]] <- t(as.matrix(yx_s[, 2:6]) - as.matrix(yx_s[, 7:16])%*%bs) %*% 
                 (as.matrix(yx_s[, 2:6]) - as.matrix(yx_s[, 7:16])%*%bs) / nrow(yx_s)
}

beta_seg   #���肳�ꂽ��A�W���s��
S_seg   #���肳�ꂽ���U�����U�s��
S_seg[[1]]
Sigma1[-2]
Sigma2[-2]
Sigma3[-2]
Sigma4[-2]

#�K���x���m�F
xs <- subset(yx, yx$segment==1)
ys <- as.matrix(xs[, 7:16]) %*% as.matrix(beta_seg[[1]])   #���肳�ꂽ�����ϐ�
yt <- xs[, 2:6]   #���̔����ϐ�
cbind(yt[, 1], ys[, 1])   
round(error <- yt - ys, 3)   #�덷

##�֐���p���Đ���
names(yx)
res2 <- by(yx, yx$segment, function(x) summary(lm(cbind(y1, y2, y3, y4, y5) ~ 
                                                        sex + age20 + age30 + age40 + age50 + 
                                                        work + stud + hwife + s_cnt, data=x)))

#�Z�O�����g���Ƃ̐��茋��
res2
res2$`1`
res2$`2`
res2$`3`
res2$`4`

####�������K���ϗʉ�A���f���̐���####
##�f�[�^�̏���
head(yx)
emyx <- yx[, -1]   #���ׂĂ̕ϐ�
round(head(emyx), 3)   #�f�[�^���m�F
emy <- yx[, 2:6]   #�����ϐ�
emx <- yx[, 7:16]   #�����ϐ�
k <- 4   #������

##�ϐ����i�[���郊�X�g
B <- matrix(0, 10, 5)
S <- matrix(0, 5, 5)

beta_em <- list(B, B, B, B)   #��A�W�����i�[���郊�X�g
S_em <- list(S, S, S, S)   #���U�����U�s����i�[���郊�X�g
r <- c(rep(0, 4))   #���������i�[���郊�X�g

#���ϗʐ��K���z�̖ޓx�֐�
dmv <- function(y, x, beta, s){
    LLo  <-  1/(sqrt((2 * pi)^nrow(s) * det(s))) * 
             exp(-t(as.matrix(y) - t(beta) %*% as.matrix(x)) %*%
             solve(as.matrix(s)) %*%
             (as.matrix(y) - t(beta) %*% as.matrix(x)) / 2)
    return(LLo)
}

##�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̒�`
LLobz <- function(yx, k, r, beta, s){
  LLind <- matrix(0, nrow=nrow(yx), ncol=k)   #�ΐ��ޓx���i�[����s��

  #���ϗʐ��`��A���f���̃Z�O�����g���Ƃ̑ΐ��ޓx
  for(kk in 1:k){
    beta_s <- beta[[kk]]
    s_s <- s[[kk]]
    Li <- apply(yx, 1, function(x) dmv(y=x[1:5], x=x[6:15], beta=beta_s, s=s_s))
    #Li <- numeric()
    #for(i in 1: 4000){
    #  Li <- c(Li, dmv(y=yx[i, 1:5], x=yx[i, 6:15], beta=beta_s, s=s_s))
    #}
    LLind[, kk] <- as.vector(Li)
  }
  LLho <- matrix(r, nrow=nrow(yx), ncol=k, byrow=T) * LLind
  z <- LLho/matrix(apply(LLho, 1, sum), nrow=nrow(yx), ncol=k)   #z�̌v�Z
  LLosum <- sum(log(apply(matrix(r, nrow=nrow(yx), ncol=k, byrow=T) * LLind, 1, sum)))   #�ϑ��f�[�^�̑ΐ��ޓx�̑��a 
  rval <- list(LLob=LLosum, z=z, LL=LLind, Li=Li)
  return(rval)
}

##EM�A���S���Y���̏����l�̐ݒ�
iter <- 0
k <- 4   #������

#�x�[�^�̏����l��ݒ�
beta_f1 <- beta_seg[[1]] + matrix(rnorm(50, 0, 1.5), 10, 5)
beta_f2 <- beta_seg[[2]] + matrix(rnorm(50, 0, 1.5), 10, 5)
beta_f3 <- beta_seg[[3]] + matrix(rnorm(50, 0, 1.5), 10, 5)
beta_f4 <- beta_seg[[4]] + matrix(rnorm(50, 0, 1.5), 10, 5)
beta <- list(beta_f1, beta_f2, beta_f3, beta_f4)


#���U�����U�s��̏����l��ݒ�
S_f1 <- S_seg[[1]] +  diag(runif(5, -7, 7))
S_f2 <- S_seg[[2]] +  diag(runif(5, -7, 7))
S_f3 <- S_seg[[3]] +  diag(runif(5, -7, 7))
S_f4 <- S_seg[[4]] +  diag(runif(5, -7, 7))
s <- list(S_f1, S_f2, S_f3, S_f4) 

#�������̏����l
r <- c(0.4, 0.2, 0.3, 0.1)

#�ΐ��ޓx�̏�����
L <- LLobz(yx=emyx, k=k, r=r, beta=beta, s=s)
L1 <- L$LLob

#�X�V�X�e�[�^�X
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 1

##EM�A���S���Y���ɂ�鐄��
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  #E�X�e�b�v�̌v�Z
  z <- L$z   #���ݕϐ�z�̏o��
  
  #M�X�e�b�v�̌v�Z
  #��A�W���̐���
  BETA <- list()
  for(i in 1:k){
    R <- diag(z[, i])
    betan <- solve(t(emyx[, 6:15]) %*% R %*% as.matrix(emyx[, 6:15])) %*% 
             t(emyx[, 6:15]) %*% R %*% as.matrix(emyx[, 1:5])
    BETA[[i]] <- betan
  }

  #���U�����U�s��̐���
  SIGMA <- list()
  for(j in 1:k){
  sigman <- (t(z[, j]*as.matrix(emyx[, 1:5]) - z[, j]*as.matrix(emyx[, 6:15]) %*% BETA[[j]]) %*%
              (z[, j]*as.matrix(emyx[, 1:5]) - z[, j]*as.matrix(emyx[, 6:15]) %*% BETA[[j]])) / sum(z[, j])
  SIGMA[[j]] <- sigman 
  }
  
  #�������̐���
  r <- apply(L$z, 2, sum) / nrow(emyx) 
  
  L <- LLobz(yx=emyx, k=k, r=r, beta=BETA, s=SIGMA)   #�ϑ��f�[�^�̑ΐ��ޓx���v�Z

  LL <- L$LLob   #�ϑ��f�[�^�̑ΐ��ޓx
  iter <- iter+1
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

BETA[[1]]   #�Z�O�����g���Ƃ̐��肳�ꂽ��A�W��
beta_seg[[1]]   #�Z�O�����g���Ƃ̍ŏ����@�Ő��肳�ꂽ��A�W��
a1   #�Z�O�����g���Ƃ̐^�̉�A�W��

SIGMA[[3]]   #�Z�O�����g���Ƃ̐��肳�ꂽ���U�����U�s��
S_seg[[3]]   #�Z�O�����g���Ƃ̍ŏ����@�Ő��肳�ꂽ���U�����U�s��
Sigma3   #�Z�O�����g���Ƃ̐^�̕��U�����U�s��

r   #���肳�ꂽ������

round(L$z[1:20, ], 3)   #�l���Ƃ̃Z�O�����g�ւ̏����m��(�^�̃Z�O�����g1)
round(L$z[1001:1020, ], 3)   #�l���Ƃ̃Z�O�����g�ւ̏����m��(�^�̃Z�O�����g2)
round(L$z[2001:2020, ], 3)   #�l���Ƃ̃Z�O�����g�ւ̏����m��(�^�̃Z�O�����g3)
round(L$z[3001:3020, ], 3)   #�l���Ƃ̃Z�O�����g�ւ̏����m��(�^�̃Z�O�����g4)
 
BETA   #���肳�ꂽ��A�W��
beta_seg   #�ŏ����@�Ő��肳�ꂽ��A�W��
a1; a2; a3; a4   #�^�̉�A�W��

L$LLob   #�ϑ��f�[�^�̑ΐ��ޓx
-2*(L$LLob) + 2*k*nrow(BETA[[1]])*ncol(BETA[[1]])   #AIC

