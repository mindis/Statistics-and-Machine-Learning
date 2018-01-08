#####Switching�}���`���[�_��LDA���f��#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(5723)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
s <- 2   #�f�[�^��
k1 <- 4   #����1�̓Ɨ������g�s�b�N
k2 <- 4   #����2�̓Ɨ������g�s�b�N
k3 <- 5   #���ʂ̃g�s�b�N
k <- k1 + k2 + k3   #���g�s�b�N��
d <- 2000   #������
v1 <- 100   #�f�[�^1�Ɋ֌W�̂���g�s�b�N�̌�b��
v2 <- 100   #�f�[�^2�Ɋ֌W�̂���g�s�b�N�̌�b��
v3 <- 100   #���ʂ̃g�s�b�N�Ɋ֌W�̂����b��
v4 <- 100   #�g�s�b�N�Ɋ֌W�̂Ȃ���b��
v <-   v1 + v2 + v3 + v4   #����b��
w1 <- rpois(d, rgamma(d, 45, 0.50))   #1����������̒P�ꐔ
w2 <- rpois(d, rgamma(d, 50, 0.50))
f1 <- sum(w1)
f2 <- sum(w2)

##�p�����[�^�̐ݒ�
#�f�B���N�����z�̃p�����[�^��ݒ�
alpha01 <- rep(0.4, k)
alpha11 <- c(rep(0.5, v1), rep(0.001, v2+v3+v4))
alpha12 <- rep(0.001, v)
alpha13 <- c(rep(0.001, v1+v2), rep(0.5, v3), rep(0.001, v4))
alpha14 <- rep(0.001, v)
alpha15 <- c(rep(0.001, v1), rep(0.5, v2), rep(0.001, v3+v4))
alpha16 <- c(rep(0.001, v1+v2), rep(0.5, v3), rep(0.001, v4))
alpha17 <- c(rep(0.01, v1+v2+v3), rep(30, v4))


#�p�����[�^�̐���
thetat <- theta <- extraDistr::rdirichlet(d, alpha01)
phit <- phi <- rbind(extraDistr::rdirichlet(k1, alpha11), extraDistr::rdirichlet(k2, alpha12), 
                     extraDistr::rdirichlet(k3, alpha13))
omega <- omegat <- rbind(extraDistr::rdirichlet(k1, alpha14), extraDistr::rdirichlet(k2, alpha15),
                         extraDistr::rdirichlet(k3, alpha16))
gammat <- gamma <- extraDistr::rdirichlet(1, alpha17)
betat <- beta <- rbeta(d, 25, 15)


##�������z����f�[�^�𐶐�
WX1 <- matrix(0, nrow=d, ncol=v)
WX2 <- matrix(0, nrow=d, ncol=v)
Z1_list <- list()
Z2_list <- list()
y1_list <- list()
y2_list <- list()

for(i in 1:d){
  ##����1�̒P��𐶐�
  #����1�̃g�s�b�N�𐶐�
  z1 <- rmnom(w1[i], 1, theta[i, ])   #�����̃g�s�b�N���z�𔭐�
  z1_vec <- as.numeric(z1 %*% c(1:k))   #�g�s�b�N�������x�N�g����
 
  #��ʌꂩ�ǂ����𐶐�
  y01 <- rbinom(w1[i], 1, beta[i])
  index_y01 <- which(y01==1)

  #�g�s�b�N����P��𐶐�
  wn1 <- rmnom(sum(y01), 1, phi[z1_vec[index_y01], ])   #����1�̃g�s�b�N����P��𐶐�
  wn2 <- rmnom(1, sum(1-y01), gamma)   #��ʌ�𐶐�
  wdn <- colSums(wn1) + colSums(wn2)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  
  WX1[i, ] <- wdn
  Z1_list[[i]] <- z1
  y1_list[[i]] <- y01
  
  
  ##����2�̒P��𐶐�
  #����2�̃g�s�b�N�𐶐�
  z2 <- rmnom(w2[i], 1, theta[i, ])   #�����̃g�s�b�N���z�𔭐�
  z2_vec <- as.numeric(z2 %*% c(1:k))   #�g�s�b�N�������x�N�g����
  
  #��ʌꂩ�ǂ����𐶐�
  y02 <- rbinom(w2[i], 1, beta[i])
  index_y02 <- which(y02==1)

  #�g�s�b�N����P��𐶐�
  wn1 <- rmnom(sum(y02), 1, omega[z2_vec[index_y02], ])   #����1�̃g�s�b�N����P��𐶐�
  wn2 <- rmnom(1, sum(1-y02), gamma)   #��ʌ�𐶐�
  wdn <- colSums(wn1) + colSums(wn2)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  WX2[i, ] <- wdn
  Z2_list[[i]] <- z2
  y2_list[[i]] <- y02
}

#���X�g�`����ϊ�
Z1 <- do.call(rbind, Z1_list)
y1_vec <- unlist(y1_list)
Z2 <- do.call(rbind, Z2_list)
y2_vec <- unlist(y2_list)


####�g�s�b�N���f������̂��߂̃f�[�^�Ɗ֐��̏���####
##�f�[�^����pID���쐬
ID1_list <- list()
wd1_list <- list()
ID2_list <- list()
wd2_list <- list()

#���l���Ƃɋ��lID����ђP��ID���쐬
for(i in 1:d){
  print(i)
  
  #����1�̒P���ID�x�N�g�����쐬
  ID1_list[[i]] <- rep(i, w1[i])
  num1 <- (WX1[i, ] > 0) * (1:v)
  num2 <- which(num1 > 0)
  W1 <- WX1[i, (WX1[i, ] > 0)]
  number <- rep(num2, W1)
  wd1_list[[i]] <- number
  
  #����2�̒P���ID�x�N�g�����쐬
  ID2_list[[i]] <- rep(i, w2[i])
  num1 <- (WX2[i, ] > 0) * (1:v)
  num2 <- which(num1 > 0)
  W1 <- WX2[i, (WX2[i, ] > 0)]
  number <- rep(num2, W1)
  wd2_list[[i]] <- number
}

#���X�g���x�N�g���ɕϊ�
ID1_d <- unlist(ID1_list)
wd1 <- unlist(wd1_list)
ID2_d <- unlist(ID2_list)
wd2 <- unlist(wd2_list)
wd <- c(wd1, wd2)

##�C���f�b�N�X���쐬
doc1_list <- list()
word1_list <- list()
doc2_list <- list()
word2_list <- list()
word_list <- list()

for(i in 1:length(unique(ID1_d))) {doc1_list[[i]] <- which(ID1_d==i)}
for(i in 1:v) {word1_list[[i]] <- which(wd1==i)}
for(i in 1:length(unique(ID2_d))) {doc2_list[[i]] <- which(ID2_d==i)}
for(i in 1:v) {word2_list[[i]] <- which(wd2==i)}
for(i in 1:v) {word_list[[i]] <- which(wd==i)}
gc(); gc()


####�}���R�t�A�������e�J�����@�őΉ��g�s�b�N���f���𐄒�####
##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #���S�W���̊i�[�p
  for(j in 1:k){
    #���S�W�����v�Z
    Bi <- rep(theta[, j], w) * phi[j, wd]   #�ޓx
    Bur[, j] <- Bi   
  }
  
  Br <- Bur / rowSums(Bur)   #���S���̌v�Z
  r <- colSums(Br) / sum(Br)   #�������̌v�Z
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##�A���S���Y���̐ݒ�
R <- 10000   #�T���v�����O��
keep <- 2   #2���1��̊����ŃT���v�����O���ʂ��i�[
iter <- 0
disp <- 20
burnin <- 1000/keep

##���O���z�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z
alpha01 <- 1.0
beta01 <- 0.25
beta02 <- 0.25
beta03 <- c(f1/20, f1/20)
beta04 <- c(f2/20, f2/20)


##�p�����[�^�̏����l
#tfidf�ŏ����l��ݒ�
idf11 <- log(nrow(WX1)/colSums(rbind(WX1, 1) > 0))
idf12 <- log(nrow(WX1)/colSums(rbind(WX1, 1)==0))
idf21 <- log(nrow(WX2)/colSums(rbind(WX2, 1) > 0))
idf22 <- log(nrow(WX2)/colSums(rbind(WX2, 1)==0))


theta <- extraDistr::rdirichlet(d, rep(1, k))   #�����g�s�b�N�̃p�����[�^�̏����l
phi <- extraDistr::rdirichlet(k, idf11*10)   #����1�̒P��g�s�b�N�̃p�����[�^�̏����l
omega <- extraDistr::rdirichlet(k, idf21*10)   # ����2�̒P��g�s�b�N�̃p�����[�^�̏����l
gamma <- extraDistr::rdirichlet(1, 1/(idf11+idf21)*10)   #��ʌ�̃p�����[�^�̏����l
r1<- 0.5; r2 <- 0.5   #�������̏����l


##�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, R/keep))
OMEGA <- array(0, dim=c(k, v, R/keep))
GAMMA <- matrix(0, nrow=R/keep, v)
SEG1 <- matrix(0, nrow=f1, ncol=k)
SEG2 <- matrix(0, nrow=f2, ncol=k)
Y1 <- rep(0, f1)
Y2 <- rep(0, f2)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"
gc(); gc()


##MCMC����p�z��
wsum0 <- matrix(0, nrow=d, ncol=k)
vf0 <- matrix(0, nrow=k, ncol=v)
wf0 <- matrix(0, nrow=k, ncol=v) 
df0 <- rep(0, v)


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##����1�̒P��g�s�b�N���T���v�����O
  #�g�s�b�N�̏o�����Ɩޓx�𐄒�
  out1 <- burden_fr(theta, phi, wd1, w1, k)
  LH1 <- out1$Bur
  word_rate1 <- out1$Br
  
  ##����2�̒P��g�s�b�N���T���v�����O
  #�g�s�b�N�̏o�����Ɩޓx�𐄒�
  out2 <- burden_fr(theta, omega, wd2, w2, k)
  LH2 <- out2$Bur
  word_rate2 <- out2$Br
  
  #���֌W�̃g�s�b�N�̏o�����Ɩޓx�𐄒�
  LH01 <- gamma[wd1]
  LH02 <- gamma[wd2]

  ##��ʌꂩ�ǂ������T���v�����O
  #����1�̃X�C�b�`���O�ϐ����T���v�����O
  Bur11 <- r1 * rowSums(LH1)
  Bur12 <- (1-r1) * LH01
  switch_rate1 <- Bur11 / (Bur11 + Bur12)
  y1 <- rbinom(f1, 1, switch_rate1)
  index_y1 <- which(y1==1)
  
  
  #����2�̃X�C�b�`���O�ϐ����T���v�����O
  Bur21 <- r2 * rowSums(LH2)
  Bur22 <- (1-r2) * LH02
  switch_rate2 <- Bur21 / (Bur21 + Bur22)
  y2 <- rbinom(f2, 1, switch_rate2)
  index_y2 <- which(y2==1)
  
  #�x�[�^���z���獬�������X�V
  par1 <- sum(y1); par2 <- sum(y2)
  r1 <- rbeta(1, par1+beta03[1], f1-par1+beta03[2])
  r2 <- rbeta(1, par2+beta04[1], f2-par2+beta04[2])
  
  ##�������z����P��g�s�b�N���T���v�����O
  #����1�̃g�s�b�N���T���v�����O
  Zi1 <- rmnom(f1, 1, word_rate1)   
  Zi1[-index_y1, ] <- 0
  z1_vec <- as.numeric(Zi1 %*% 1:k)
  
  #����2�̃g�s�b�N���T���v�����O
  Zi2 <- rmnom(f2, 1, word_rate2)   
  Zi2[-index_y2, ] <- 0
  z12_vec <- as.numeric(Zi2 %*% 1:k)
  
  ##�p�����[�^���T���v�����O
  #�g�s�b�N���z���f�B�N�������z����T���v�����O
  wsum01 <- matrix(0, nrow=d, ncol=k)
  wsum02 <- matrix(0, nrow=d, ncol=k)
  for(i in 1:d){
    wsum01[i, ] <- colSums(Zi1[doc1_list[[i]], ])
    wsum02[i, ] <- colSums(Zi2[doc2_list[[i]], ])
  }
  wsum <- wsum01  + alpha01
  theta <- extraDistr::rdirichlet(d, wsum)
  
  #�g�s�b�N��̕��z���f�B�N�������z����T���v�����O
  vf0 <- matrix(0, nrow=k, ncol=v)
  wf0 <- matrix(0, nrow=k, ncol=v)
  for(j in 1:v){
    vf0[, j] <- colSums(Zi1[word1_list[[j]], , drop=FALSE])
    wf0[, j] <- colSums(Zi2[word2_list[[j]], , drop=FALSE])
  }
  
  vf <- vf0 + beta01
  wf <- wf0 + beta01
  phi <- extraDistr::rdirichlet(k, vf)
  omega <- extraDistr::rdirichlet(k, wf)
  
  #��ʌ�̕��z���f�B�N�������z����T���v�����O
  y1_zeros <- 1-y1
  y2_zeros <- 1-y2
  df01 <- rep(0, v)
  df02 <- rep(0, v)
 
  for(j in 1:v){
    df01[j] <- sum(y1_zeros[word1_list[[j]]])
    df02[j] <- sum(y2_zeros[word2_list[[j]]])
  }
  df1 <- df01 + df02 + beta02
  gamma <- extraDistr::rdirichlet(1, df1)

  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    OMEGA[, , mkeep] <- omega
    GAMMA[mkeep, ] <- gamma
    
    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(rp%%keep==0 & rp >= burnin){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
      Y1 <- Y1 + y1
      Y2 <- Y2 + y2
    }
    
    if(rp%%disp==0){
      #�T���v�����O���ʂ��m�F
      print(rp)
      print(c(mean(y1), mean(y1_vec)))
      print(c(mean(y2), mean(y2_vec)))
      #print(round(cbind(theta[1:10, ], thetat[1:10, ]), 3))
      print(round(cbind(phi[, 296:305], phit[, 296:305]), 3))
      print(round(rbind(gamma[296:305], gammat[296:305]), 3))
    }
  }
}


####�T���v�����O���ʂ̉����Ɨv��####
burnin <- 1000/keep   #�o�[���C������
RS <- R/keep

##�T���v�����O���ʂ̉���
#�����̃g�s�b�N���z�̃T���v�����O����
matplot(t(THETA[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA[100, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA[1000, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA[2000, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

#�P��̏o���m���̃T���v�����O����
matplot(t(PHI[1, 296:305, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N1�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[2, 296:305, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N2�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[3, 296:305, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N3�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[4, 296:305, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N4�̒P��̏o�����̃T���v�����O����")
matplot(t(OMEGA[1, 296:305, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N1�̒P��̏o�����̃T���v�����O����")
matplot(t(OMEGA[2, 296:305, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N2�̒P��̏o�����̃T���v�����O����")
matplot(t(OMEGA[3, 296:305, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N3�̒P��̏o�����̃T���v�����O����")
matplot(t(OMEGA[4, 296:305, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N4�̒P��̏o�����̃T���v�����O����")

#��ʌ�̏o���m���̃T���v�����O����
matplot(GAMMA[, 286:295], type="l", ylab="�p�����[�^", main="�P��̏o�����̃T���v�����O����")
matplot(GAMMA[, 296:305], type="l", ylab="�p�����[�^", main="�P��̏o�����̃T���v�����O����")
matplot(GAMMA[, 306:315], type="l", ylab="�p�����[�^", main="�P��̏o�����̃T���v�����O����")


##�T���v�����O���ʂ̗v�񐄒��
#�g�s�b�N���z�̎��㐄���
topic_mu <- apply(THETA[, , burnin:(R/keep)], c(1, 2), mean)   #�g�s�b�N���z�̎��㕽��
round(cbind(topic_mu, thetat), 3)
round(topic_sd <- apply(THETA[, , burnin:(R/keep)], c(1, 2), sd), 3)   #�g�s�b�N���z�̎���W���΍�

#�P��o���m���̎��㐄���
word_mu1 <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #�P��̏o�����̎��㕽��
word1 <- round(t(rbind(word_mu1, phit)), 3)
word_mu2 <- apply(OMEGA[, , burnin:(R/keep)], c(1, 2), mean)   #�P��̏o�����̎��㕽��
word2 <- round(t(rbind(word_mu2, omegat)), 3)
word <- round(t(rbind(word_mu1, word_mu2, phit, omegat)), 3)
colnames(word) <- 1:ncol(word)

word_mu3 <- apply(GAMMA[burnin:(R/keep), ], 2, mean)   #�P��̏o�����̎��㕽��
round(rbind(word_mu3, gamma=gammat), 3)


##�g�s�b�N�̎��㕪�z�̗v��
round(seg1_mu <- SEG1 / rowSums(SEG1), 3)
round(seg2_mu <- SEG2 / rowSums(SEG2), 3)




