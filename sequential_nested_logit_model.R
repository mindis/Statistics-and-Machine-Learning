#####�����^�l�X�e�b�h���W�b�g���f��#####
library(MASS)
library(mlogit)
library(nnet)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

#####�f�[�^�̔���#####
#set.seed(89361)
##�f�[�^�̐ݒ�
hh <- 1000   #�v���C���[��
pt <- 120   #�ϑ�����
hhpt <- hh*pt   #�S�T���v����

#ID�̐ݒ�
ID <- matrix(0, hhpt, 3)
ID[, 1] <- 1:nrow(ID)
ID[, 2] <- rep(1:hh, rep(pt, hh))
ID[, 3] <- rep(1:pt, hh)


##�p�����[�^�̐ݒ�
#���O�C���L���̃p�����[�^�̐ݒ�
login.para <- 7   #���O�C�����f���̃p�����[�^��
alpha.s0 <- runif(1, -1.6, -1.3)   #�ؕ�
alpha.s1 <- runif(1, 0.9, 1.5)   #�ߋ�60���̃��O�C����
alpha.s2 <- runif(1, -0.4, -0.3)   #�O�񂩂�̃��O�C���o�ߎ���(�P��=��)
alpha.s3 <- runif(1, 0.5, 1.2)   #�O���̃��O�C���L��
alpha.s4 <- runif(1, 0.13, 0.24)   #1��������̃v���C�����ς̑ΐ�
alpha.s5 <- runif(1, 0.5, 1.0)   #�C�x���g�L��
alpha.s6 <- runif(1, 0.3, 0.6)   #�v���C���[�̃X�L�����x��
alpha.s7 <- runif(1, 0.4, 0.9)   #���O�T���ϐ��̃p�����[�^
alpha.s <- c(alpha.s1, alpha.s2, alpha.s3, alpha.s4, alpha.s4, alpha.s6, alpha.s7)

#�K�`���L���̃p�����[�^�̐ݒ�
buy.para <- 14   #�ۋ��L���̃p�����[�^��
beta.s0 <- runif(1, -4.2, -3.6)   #�ؕ�
beta.s1 <- runif(1, 0.7, 1.2)   #�ߋ�60���̃K�`����
beta.s2 <- runif(1, -0.7, -0.5)   #�ߋ�7���̃K�`���L��
beta.s3 <- runif(1, 0.15, 0.3)   #1�񂠂���̃K�`����
beta.s4 <- runif(1, 0.3, 0.6)   #�O��̃K�`������̌o�ߓ����̑ΐ�
beta.s5 <- runif(1, 0.1, 0.22)   #1��������̃v���C�����ς̑ΐ�
beta.s6 <- runif(9, 0.8, 1.6)   #UR�̊m���A�b�v�L����
beta.s <- c(beta.s1, beta.s2, beta.s3, beta.s4, beta.s5, beta.s6)

##�����l�̐ݒ�
X.login <- matrix(0, nrow=hhpt, ncol=login.para)
X.buy <- matrix(0, nrow=hhpt, ncol=buy.para)

#1���ڂ����o��
index.f <- subset(1:nrow(X.login), ID[, 3]==1)

##�����f�[�^�̔���
##���O�C���f�[�^�̏����l
#�ߋ�60���̃��O�C������
login.hist <- matrix(0, nrow=hh, ncol=60)
for(i in 1:hh){
 p <- runif(1, 0.2, 1.0) 
 login.hist[i, ] <- rbinom(60, 1, p)
}
login.ratef <- rowMeans(login.hist)
X.login[index.f, 1] <- login.ratef   #�ߋ�60���̃��O�C����

#�O�񂩂�̃��O�C���o�ߎ���
login.last <- apply(login.hist, 1, function(x) tail(subset(1:length(x), x==1), 1))
X.login[index.f, 2] <- 61-login.last   #�O�񂩂�̃��O�C���o�ߎ���   

#�O�����O�C���̗L��
login.pre <- ifelse(61-login.last==1, 1, 0)
X.login[index.f, 3] <- login.pre

#1��������̃v���C������
play.M <- matrix(0, nrow=hh, ncol=60) 

for(i in 1:hh){
  #���σv���C�񐔂̋L�^
  pois <- runif(1, 5, 18)
  play <- rpois(sum(login.hist[i, ]), pois)
  X.login[index.f[i], 4] <- log(sum(play)/sum(login.hist[i, ]))
  
  #�v���C�����̋L�^
  index.play <- subset(1:length(login.hist[i, ]), login.hist[i, ]==1)
  play.M[i, index.play] <- play
}  
hist(X.login[index.f, 4], breaks=20, col="grey", xlab="1��������v���C������",main="1��������v���C�����ς̕��z")

#�C�x���g�̗L��
X.login[, 5] <- rep(c(rep(0, 5), rep(1, 10)), hh)


#�v���C���[�X�L��
for(i in 1:hh){
  X.login[ID[, 2]==i, 6] <- rnorm(1, 0, 1)   
}

##�K�`���f�[�^�̏����l
#�ߋ�60���̃K�`����
buy.hist <- matrix(0, nrow=hh, ncol=60)

for(i in 1:hh){
  logit.bpast <- beta.s0 + runif(1, 0.5, 3.0)
  Pr.gpast <- exp(logit.bpast)/(1+exp(logit.bpast))
  buy.hist[i, ] <- rbinom(60, 1, Pr.gpast)
}
X.buy[index.f, 1]  <- rowSums(buy.hist)/60

#�ߋ�7���̃K�`���L��
X.buy[index.f, 2] <- ifelse(rowSums(buy.hist[, 54:60]) > 0, 1, 0)   


#1�񂠂���̃K�`����
cnt.zeros <- rpois(hh, 2.5)
g.cnt <- ifelse(cnt.zeros==0, 1, cnt.zeros)
g.cnt[index.pp==0] <- 0
X.buy[, 3] <- g.cnt

#1��������̃v���C������
buy.M <- matrix(0, nrow=hh, ncol=60) 

for(i in 1:hh){
  #���σK�`���񐔂̋L�^
  
  play <- rpois(sum(login.hist[i, ]), pois)
  X.login[index.f[i], 4] <- log(sum(play)/sum(login.hist[i, ]))
  
  #�v���C�����̋L�^
  index.play <- subset(1:length(login.hist[i, ]), login.hist[i, ]==1)
  play.M[i, index.play] <- play
}  
hist(X.login[index.f, 4], breaks=20, col="grey", xlab="1��������v���C������",main="1��������v���C�����ς̕��z")


#�O��̃K�`������̌o�ߓ���
buy.progress <- rep(0, hh)
index.buy <- subset(1:nrow(buy.hist), rowSums(buy.hist) > 0)
buy.last <- apply(buy.hist[index.buy, ], 1, function(x) tail(subset(1:length(x), x==1), 1))
buy.progress[index.buy] <- 61-buy.last
X.buy[index.f, 4] <- ifelse(buy.progress > 0, log(buy.progress), 0)   #�O�񂩂�̃��O�C���o�ߎ���  


day.r <- c()
for(i in 1:hh){
  day.r <- c(day.r, (pt+1) - max(sample(1:pt, index.pp[i])))
}
day.rl <- log(day.r)
day.rl[is.infinite(day.rl)==TRUE] <- 0
X.buy[index.f, 5] <- day.rl


#1��������̃v���C���̑ΐ�
X.buy[index.f, 6] <- X.login[index.f, 4]

#UR�̊m���A�b�v�L����
#�m���A�b�v���Ƃ����ł͂Ȃ��������
cycle <- pt/15
r.no <- matrix(0, nrow=cycle, 5)
r.new1 <- matrix(0, nrow=cycle, 5)
r.new2 <- matrix(0, nrow=cycle, 5)

for(c in 1:cycle){
  r.no[c, ] <- (c-1) * 10 + ((c-1)*5+1):((c-1)*5+5)
  r.new1[c, ] <- (c-1) * 10 + ((c-1)*5+6):((c-1)*5+10)
  r.new2[c, ] <- (c-1) * 10 + ((c-1)*5+11):((c-1)*5+15)
}
r.cycle <- list(r.no, r.new1, r.new2)


#UR�m���A�b�v�L�����𔭐�������
#UR�m���A�b�v�L�����̂��炩���ߔ��������Ă���
UR <- matrix(0, nrow=cycle*(length(r.cycle)-1), 9)

ur <- rep(0, length(beta.s6))
index.ur <- sample(1:length(beta.s6), 2)
ur[index.ur] <- 1
UR[1, ] <- ur

#�O���UR�L�����Ƃ��Ԃ�Ȃ��悤�ɂ��Ă���
for(c in 1:(cycle*(length(r.cycle)-1))){
  for(i in 1:500){
    ur <- rep(0, length(beta.s6))
    index.ur <- sample(1:length(beta.s6), 2)
    ur[index.ur] <- 1
    UR[c, ] <- ur
    if(max(UR[c, ]+UR[c-1, ]) == 1) {break}  
  }
}


#�f�[�^�Z�b�g��UR�L�������L�^
for(t in 1:cycle){
  for(c in 1:length(r.cycle)){
    if(c==1) {next} else {
      tf <- r.cycle[[c]][t, ]
    }  
    #UR�L�������L�^
    X.buy[ID[, 3] %in% tf, 6:ncol(X.buy)] <- matrix(UR[2*t-2+(c-1), ], nrow=sum(ID[, 3] %in% tf), ncol=9, byrow=T)
  }
}
cbind(ID, X.buy[, 6:ncol(X.buy)])[1:120, ]   #�f�[�^���m�F 


##�f�[�^�̏����l���m�F
round(cbind(ID[index.f, ], X.login[index.f, ]), 3)   #���O�C���f�[�^
round(cbind(ID[index.f, ], X.buy[index.f, ]), 3)   #�K�`���f�[�^


####�p�l�����Ƃɒ����I�Ƀf�[�^�𔭐�������####
##�f�[�^�X�V�p�̕ۑ��z��
login.hist <- cbind(login.hist, matrix(0, nrow=hh, ncol=pt))
play.M <- cbind(play.M, matrix(0, nrow=hh, ncol=pt))
buy.hist <- cbind(buy.hist, matrix(0, nrow=hh, ncol=pt))


##���W�b�g�̒�`
#�K�`���̃��W�b�g
logit.g <- beta.s0 + X.buy[ID[, 2]==1 & ID[, 3]==1, ] %*% beta.s   

#���O�T���ϐ��̒�`
logsum <- log(1+exp(logit.g))
X.login[ID[, 2]==1 & ID[, 3]==1, 7] <- logsum

#���O�C���̃��W�b�g
logit.l <- alpha.s0 + X.login[ID[, 2]==1 & ID[, 3]==1, ] %*% alpha.s   #���O�C���̃��W�b�g

##�m�����v�Z���ĉ����ϐ��𔭐�
Y <- matrix(0, nrow=hhpt, ncol=2)

#���O�C���̔���
Pr.l <- exp(logit.l)/(1+exp(logit.l))   #�m�����v�Z
Y[ID[, 2]==1 & ID[, 3]==1, 1] <- rbinom(1, 1, Pr.l)

#�K�`���̔���
#���O�C�����������ꍇ�ɂ̂ݔ���
Pr.g <- exp(logit.g)/(1+exp(logit.g))   #�m�����v�Z
if(Y[ID[, 2]==1 & ID[, 3]==1, 1]==1) {
  Y[ID[, 2]==1 & ID[, 3]==1, 2] <- rbinom(1, 1, Pr.g)} else {
    Y[ID[, 2]==1 & ID[, 3]==1, 2] <- 0
}



##�f�[�^�̍X�V
#���O�C�������̍X�V

