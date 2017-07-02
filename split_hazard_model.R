#####split hazard model#####
library(MASS)
library(survival)
library(reshape2)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
zk <- 2   #���ݕϐ���
n <- 2000   #�T���v����
pt <- rpois(n, 10)   #�T���v�����Ƃ̊ϑ���
pt <- ifelse(pt < 1, 1, pt)
nmax <- sum(pt)   #�ő�T���v����

##���ݕϐ��̔���
#�f���O���t�B�b�N�ϐ��̔���
#�A���ϐ��̔���
cont_z <- matrix(rnorm(n*2), nrow=n, ncol=2)

#��l�ϐ��̔���
bin_z <- matrix(0, nrow=n, ncol=3)
pbin <- c(0.4, 0.3, 0.6)
for(i in 1:length(pbin)){
  bin_z[, i] <- rbinom(n, 1, pbin[i])
}

#���l�ϐ��̔���
pmulti <- c(0.2, 0.1, 0.3, 0.3, 0.1)
multi_z <- t(rmultinom(n, 1, pmulti))
multi_z <- multi_z[, -5]

X_z <- cbind(cont_z, bin_z, multi_z)

#�p�����[�^�̐ݒ�
bz1 <- runif(ncol(cont_z), -0.7, 1.0)
bz2 <- runif(ncol(bin_z), -0.9, 1.2)
bz3 <- runif(ncol(multi_z), -1.1, 1.5)
bz0 <- c(0.5)
bz_t <- c(bz1, bz2, bz3)

#���p�֐��̒�`
U <- bz0[1] + cbind(cont_z, bin_z, multi_z) %*% bz_t

#�m���̌v�Z�Ɛ��ݕϐ�z�̔���
Pr_z <- exp(U) / (1+exp(U))
z.vec <- apply(cbind(Pr_z, 1-Pr_z), 1, function(x) rbinom(1, 1, x[1]))

#���ݕϐ��̗v��
sum(z.vec)
round(mean(z.vec), 3)

##ID�Ɛ��ݕϐ��̐ݒ�
#���ݕϐ��̊���
z.index <- c()
for(i in 1:n){
  z.index <- c(z.index, rep(z.vec[i], pt[i]))
}

#id�̊���
id <- c()
for(i in 1:n){
  id <- c(id, rep(i, pt[i]))
}

#time�̊���
time <- c()
for(i in 1:n){
  time <- c(time, 1:pt[i])
}

##�n�U�[�h���f���̐����ϐ��̔���
page_cnt <- 10

##�y�[�W�{���񐔂Ɖ{�������̔���
#�y�[�W�{���񐔂̔���
lam_lower <- 5
lam_upper <- 9
p_cnt.zero <- rpois(nmax, runif(nmax, lam_lower, lam_upper))
p_cnt <- ifelse(p_cnt.zero==0, 1, p_cnt.zero)
hist(p_cnt, breaks=15, col="grey", xlab="�y�[�W�{����", main="�y�[�W�{�����̕��z")

#�y�[�W�{�������̔���
p_rate <- runif(page_cnt)
p_hist <- matrix(0, nmax, page_cnt)

for(i in 1:nmax){
  p_hist[i, ] <- t(rmultinom(1, p_cnt[i], p_rate))
}
p_hist
p_hist.r <- p_hist / rowSums(p_hist)

#�Ō�Ɍ��Ă����y�[�W�̔���
p_last <- t(rmultinom(nmax, 1, p_rate))

#�O��̃y�[�W�{���ł����Ƃ����Ă����y�[�W
index.m <- subset(1:nrow(p_hist), apply(p_hist, 1, max) > 1)

pm <- t(apply((p_hist-2), 1, function(x) x-abs(max(x))))
p_most1 <- ifelse(pm==0, 1, 0)

#1�Ԗڂ̃A�N�Z�X�͂��ׂ�0�ɂȂ�
p_most2 <- rbind(0, p_most1)
p_most2[time==1, ] <- rep(0, page_cnt)
p_most <- p_most2[1:nrow(p_hist), ]


##�O�񂩂�̃A�N�Z�X�o�ߎ���(�P��=��)�̔���
shape <- 1.65
rate <- 0.6
t <- round(rgamma(nmax, shape, rate), 0)
index.t <- subset(1:length(t), t == 0)
t[index.t] <- round(runif(length(index.t), 1, 5), 0)


##�璷�ȕϐ����폜���ăf�[�^������
index.h <- which.min(colSums(p_hist))
ph_hist <- p_hist[, -index.h]
ph_hist.r <- p_hist.r[, -index.h]
ph_last <- p_last[, -index.h]
p_most <- p_most[, -index.h]

X <- data.frame(time=time, page=ph_hist.r, page_l=ph_last, page_m=p_most, cnt_p=p_cnt, t=t)   #�f�[�^�̌���
round(X, 3)

##�p�����[�^�̐ݒ�
beta1 <- rep(0, ncol(X)+1)
beta2 <- c(-1.55, -0.08, runif(page_cnt-1, -0.6, 0.8), runif(2*(page_cnt-1), -0.6, 0.8), 0.07, -0.14)
betat <- beta2

##�����ϐ��̔���
#���p�֐��̒�`
M_prov1 <- as.matrix(cbind(1, X)) %*% betat
U_prov1 <- M_prov1 + rnorm(nmax, 0, 1)
Y_prov1 <- ifelse(U_prov1 > 0, 1, 0)

#ID���Ƃ̃Z�O�����g�ɉ����ϐ���Ή�����y=1�����������ꍇ�ł��؂�
#�����ϐ����m��
Y_prov2 <- ifelse(z.index==0, 0, Y_prov1)   #�Z�O�����g�ŉ����ϐ���Ή�������
U_prov2 <- ifelse(z.index==0, 0, U_prov1)   #�Z�O�����g�Ō��p�֐���Ή�������
M_prov2 <- ifelse(z.index==0, 0, M_prov1)   #�Z�O�����g�Ō��p�֐��̕��ς�Ή�������

round(YZ_prov <- data.frame(id=id, tz=time, z=z.index, Y=Y_prov2, U=U_prov2, M=M_prov2), 3)   #�f�[�^�̌����Ɗm�F

#�ł��؂�ϐ���ݒ�
censored <- matrix(3, nrow(YZ_prov), ncol=1)

for(i in 1:n){
  index.ind <- subset(1:nrow(YZ_prov), YZ_prov$id==i)
  index.y <- subset(1:nrow(YZ_prov), YZ_prov$id==i & YZ_prov$Y==1)[1]
  
  if(is.na(index.y)==FALSE & index.y==index.ind[1]) {
    print("cv����")
    censored[index.y, ] <- 1
    censored[(index.y+1):index.ind[length(index.ind)], ] <- 2
    
  } else if(is.na(index.y)==FALSE & index.y!=index.ind[1]) { 
    print("cv����")
    censored[index.y, ] <- 1 
    censored[index.ind[1]:(index.y-1), ] <- 0
    censored[(index.y+1):index.ind[length(index.ind)], ] <- 2
  }
  else {
    print("cv�Ȃ�")
  }
}

##�ł��؂肵���f�[�^�͎�菜��
YZ <- data.frame(YZ_prov, censored)
index.censor <- subset(1:nrow(YZ), censored==2)
YX <- cbind(YZ, X)[-index.censor, ]
YX$censored[YX$censored==3] <- 0
Y <- YX$Y; X <- X[-index.censor, ]; id <- YX$id

#�[���ޓx�̂��߂�y�̐ݒ�
index.id <- subset(id, Y==1)
Yn <- ifelse(id %in% index.id, 1, 0)

table(z.vec)
table(Y)


####EM�A���S���Y����Split_Hazard_model�𐄒�####
##�������W�b�g���f���̑ΐ��ޓx�֐�
mlogit_LL <- function(x, Z, X_z, r){
  #�p�����[�^�̐ݒ�
  b10 <- x[r[1]]
  b11 <- x[r[2]:r[3]]
  b20 <- x[r[4]]
  b21 <- x[r[5]:r[6]]
  
  #���p�֐��̒�`
  U1 <- b10 + X_z %*% b11
  U2 <- b20 + X_z %*% b21
  U3 <- 0
  U <- cbind(U1, U2, U3)   #���p�֐�������
  
  #�ΐ��ޓx���v�Z
  LLi <- rowSums(Z * U) - log(rowSums(exp(U)))
  LL <- sum(LLi)
  return(LL)
}

##�v���r�b�g���f���̑ΐ��ޓx�̒�`
probit_LL <- function(x, Y, X){
  #�p�����[�^�̐ݒ�
  b0 <- x[1]
  b1 <- x[2:(ncol(X)+1)]
  
  #���p�֐��̒�`
  U <- b0 + as.matrix(X) %*% b1
  
  #�ΐ��ޓx���v�Z
  Pr <- pnorm(U)   #�m��
  LLi <- Y*log(Pr) + (1-Y)*log(1-Pr)
  LL <- sum(LLi)
  return(LL)
}

##���S�f�[�^�ł�Split_hazard_mode�̑ΐ��ޓx
#�p�����[�^�̐ݒ�
cll <- function(b, Y, Yn, X, zpt, zk){
  beta0 <- b[1]
  beta1 <- b[2:(ncol(X)+1)] 
  
  #���p�֐����`
  U <- beta0 + as.matrix(X) %*% beta1
  
  #�v���r�b�g���f���̊m���Ƒΐ��ޓx�̌v�Z
  Pr <- pnorm(U)
  LLc <- Y*log(Pr) + (1-Y)*log(1-Pr)
  
  #�[���ޓx�̌v�Z
  LLzero <- dbinom((1-Yn), 1, 1)   #�[���ޓx�̌v�Z
  
  LL <- sum(zpt[, 1]*LLzero + zpt[, 2]*LLc)
  return(LL)
}


##�ϑ��f�[�^�ł̖ޓx�Ɛ��ݕϐ�z�̌v�Z
ollz <- function(b, Y, Yn, X, r, id, n, zk){
  #�p�����[�^�̐ݒ�
  beta0 <- b[1]
  beta1 <- b[2:(ncol(X)+1)] 
  
  #�Z�O�����g���Ƃ̌��p�֐����`
  U <- beta0 + as.matrix(X) %*% beta1
  
  #�v���r�b�g���f���̊m���Ƒΐ��ޓx�̌v�Z
  Pr <- pnorm(U)   #�m���̌v�Z
  LLi <- Pr^Y * (1-Pr)^(1-Y)   #�ޓx�̌v�Z
  LLzero <- dbinom((1-Yn), 1, 1)   #�[���ޓx�̌v�Z
  LCo <- cbind(LLzero, LLi)   #�ޓx�̌���
  
  #ID�ʂɖޓx�̐ς����
  LLho <- matrix(0, nrow=n, ncol=zk)
  for(i in 1:n){
    if(sum(id==i)==1){
      LLho[i, ] <- LCo[YX$id==i, ]
    } else {
      LLho[i, ] <- apply(LCo[id==i, ], 2, prod) 
    }
  }
  
  #�ϑ��f�[�^�ł̑ΐ��ޓx
  LLo <- sum(log(apply(r * LLho, 1, sum)))
  
  #���ݕϐ�z�̌v�Z
  z0 <- r * LLho   #���ݕϐ�z�̕��q
  z1 <- z0 / rowSums(z0)   #���ݕϐ�z�̌v�Z
  
  rval <- list(LLo=LLo, z1=z1)
  return(rval)
}


##EM�A���S���Y���̏����l�̐ݒ�
#�f�[�^�̐ݒ�
index.cv <- subset(YX$id, YX$Y==1)
Yf <- YX[YX$id %in% index.cv, "Y"]
Xf <- X[YX$id %in% index.cv, ]
zf <- as.matrix(data.frame(id=YX$id, Y) %>%
                  dplyr::group_by(id) %>%
                  dplyr::summarise(z=sum(Y)))


#�v���r�b�g���f���Ő������f���̃p�����[�^�̏����l��ݒ�
for(i in 1:1000){
  #�����p�����[�^�̐ݒ�
  print(i)
  x <- c(-0.5, runif(ncol(X)-2, -1, 1), runif(2, -0.2, 0.2)) 
  
  #���j���[�g���@�ōő剻
  res.b <- try(optim(x, probit_LL, Y=Yf, X=Xf, method="BFGS", hessian=FALSE, 
                     control=list(fnscale=-1)), silent=TRUE)
  if(class(res.b) == "try-error") {next} else {break}   #�G���[����
}
beta <- res.b$par   #�������f���̏����l


#�v���r�b�g���f���ō�����r�̃p�����[�^�̏����l��ݒ�
for(i in 1:1000){
  #�����p�����[�^�̐ݒ�
  print(i)
  x <- c(runif(ncol(X_z)+1)) 
  
  #���j���[�g���@�ōő剻
  res.z <- try(optim(x, probit_LL, Y=zf[, 2] , X=X_z, method="BFGS", hessian=FALSE, 
                     control=list(fnscale=-1)), silent=TRUE)
  if(class(res.z) == "try-error") {next} else {break}   #�G���[����
}
bz <- res.z$par + c(runif(1, 0, 1.5), runif(ncol(X_z), -0.2, 0.5))
Pz <- pnorm(cbind(1, X_z) %*% bz)    #�m���̌v�Z
r <- cbind(1-Pz, Pz)   #�������̏����l


#���ݕϐ�z�Ɗϑ��f�[�^�̑ΐ��ޓx�̏����l�̐ݒ�
oll <- ollz(b=beta, Y=Y, Yn=Yn, X=X, r=r, id=id, n=n, zk=zk)
z <- oll$z1
LL1 <- oll$LLo

#EM�A���S���Y���̐ݒ�
iter <- 0
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l��ݒ�
tol <- 1   
zpt <- matrix(0, nrow=nrow(X), ncol=zk)


##EM�A���S���Y���ɂ��Split_Hazard_model�̐���
while(abs(dl) >= tol){   #dl��tol�ȏ�Ȃ�J��Ԃ�
  #���ݕϐ�z���p�l���`���ɕύX
  for(i in 1:n){
    zpt[id==i,] <- matrix(z[i, ], nrow=length(id[id==i]), ncol=zk, byrow=T)
  }
  
  #���S�f�[�^�ł̐������f���̍Ŗސ���(M�X�e�b�v)
  res <- optim(beta, cll, Y=Y, Yn=Yn, X=X, zpt=zpt, zk=zk, method="BFGS", hessian=FALSE, control=list(fnscale=-1))
  beta <- res$par   #�p�����[�^�̍X�V
  
  #������r�̍X�V
  res.z <- optim(bz, probit_LL, Y=z[, 2], X=X_z, method="BFGS", hessian=FALSE, control=list(fnscale=-1))
  bz <- res.z$par   #���ݕϐ��̃p�����[�^
  Pz <- pnorm(cbind(1, X_z) %*% bz)    #�m���̌v�Z
  r <- cbind(1-Pz, Pz)   #�������̍X�V
    
  #E�X�e�b�v�ł̑ΐ��ޓx�̊��Ғl�̌v�Z
  oll <- ollz(b=beta, Y=Y, Yn=Yn, X=X, r=r, id=id, n=n, zk=zk)
  z <- oll$z1
  LL <- oll$LLo
  
  #EM�A���S���Y���̃p�����[�^�̍X�V
  iter <- iter+1
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

####���茋�ʂƗv��####
##���肳�ꂽ�p�����[�^�Ɛ^�̃p�����[�^�̔�r
#�������f���̉�A�p�����[�^
round(beta, 2)
round(betat, 2)

#���ݕϐ�z�̉�A�p�����[�^
round(bz, 2)
round(c(bz0, bz_t), 2)

#���ݕϐ�z�Ɛ^�̐��ݕϐ��̔�r
round(data.frame(z=z, p1=1-Pz, p2=Pz, zt=z.vec, y=zf[, 2]), 3)

##AIC��BIC�̌v�Z
round(LL, 3)   #�ő剻���ꂽ�ϑ��f�[�^�̑ΐ��ޓx
round(AIC <- -2*LL + 2*(length(beta)), 3)   #AIC
round(BIC <- -2*LL + log(nrow(X))*length(beta), 3) #BIC