#####����LDA���f��#####
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

#set.seed(2578)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
k1 <- 10
k2 <- 4
d <- 2000   #������
v1 <- 250   #���e�Ɋ֌W�̂����b��
v2 <- 100   #���e�Ɋ֌W�̂Ȃ���b��
v <- v1 + v2   #��b�� 
s <- rpois(d, 14)   #���͐�
s[s < 5] <- ceiling(runif(sum(s < 5), 5, 10))
a <- sum(s)   #�����͐�
w <- rpois(a, 11.5)   #���͂�����̒P�ꐔ
w[w < 5] <- ceiling(runif(sum(w < 5), 5, 10))
f <- sum(w)   #���P�ꐔ

#����ID�̐ݒ�
u_id <- rep(1:d, s)
t_id <- c()
for(i in 1:d){t_id <- c(t_id, 1:s[i])}
words <- as.numeric(tapply(w, u_id, sum))


##�p�����[�^��ݒ�
#�f�B���N�����z�̃p�����[�^
alpha0 <- rep(0.15, k1)
alpha1 <- c(rep(0.5, v1), rep(0.005, v2))
alpha2 <- rep(1, k2)
alpha3 <- c(rep(0.005, v1), rep(0.5, v2))

#�f�B���N�����z���p�����[�^�𐶐�
thetat <- theta <- extraDistr::rdirichlet(d, alpha0)
phit <- phi <- extraDistr::rdirichlet(k1, alpha1)
omegat <- omega <- matrix(1/k2, nrow=d, ncol=k2)
gammat <- gamma <- extraDistr::rdirichlet(k2, alpha3)
delta <- deltat <- rbeta(d, 20, 15)


##���͂��ƂɒP��𐶐�����
WX <- matrix(0, nrow=a, ncol=v)
y_list <- list()
Z1_list <- list()
Z2_list <- list()
index_v1 <- 1:v1
index_v2 <- (v1+1):v 

for(i in 1:a){
  
  ##���͂��ƂɃg�s�b�N�𐶐�
  id <- u_id[i]
  z1 <- rmnom(w[i], 1, theta[id, ])
  z1_vec <- as.numeric(z1 %*% 1:k1)
  
  ##��ʌ�̔���ƈ�ʌ�g�s�b�N�̐���
  #��ʌꂩ�ǂ����̐���
  y0 <- rbinom(w[i], 1, delta[id])
  index <- which(y0==1)
  
  #���͂̃g�s�b�N�𐶐�
  z2 <- rmnom(1, 1, omega[id, ])
  z2_vec <- as.numeric(z2 %*% 1:k2)
  
  #�g�s�b�N���z�Ɋ�Â��P��𐶐�
  n1 <- length(index)
  n2 <- w[i]-length(index)
  w2 <- w1 <- matrix(0, nrow=1, ncol=v)
  
  if(n1 > 0){
    w1 <- rmnom(n1, 1, phi[z1_vec[index], ])   #�g�s�b�N��𐶐�
  }
  if(n2 > 0){
    w2 <- rmnom(n2, 1, gamma[z2_vec, ])   #��ʌ�𐶐�
  }
  
  #bag of words�s����쐬
  wdn <- colSums(w1) + colSums(w2)
  WX[i, ] <- wdn
  
  #�p�����[�^���i�[
  Z1_list[[i]] <- z1
  Z2_list[[i]] <- z2
  y_list[[i]] <- y0 
}


#���X�g�`����ϊ�
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)
y_vec <- unlist(y_list)


####�g�s�b�N���f������̂��߂̃f�[�^�Ɗ֐��̏���####
##���ꂼ��̕������̒P��̏o������ѕ⏕���̏o�����x�N�g���ɕ��ׂ�
##�f�[�^����pID���쐬
ID_list <- list()
td1_list <- list()
td2_list <- list()
wd_list <- list()

#����id���Ƃɕ���id����ђP��id���쐬
for(i in 1:a){
  
  #����ID���L�^
  ID_list[[i]] <- rep(u_id[i], w[i])
  td1_list[[i]] <- rep(i, w[i])
  td2_list[[i]] <- rep(t_id[i], w[i])
  
  #�P��ID���L�^
  num1 <- WX[i, ] * 1:v
  num2 <- which(num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
}

#���X�g���x�N�g���ɕϊ�
ID_d <- unlist(ID_list)
td1_d <- unlist(td1_list)
td2_d <- unlist(td2_list)
wd <- unlist(wd_list)

##�C���f�b�N�X���쐬
doc_list <- list()
id_list <- list()
sent_list <- list()
word_list <- list()
for(i in 1:length(unique(ID_d))){doc_list[[i]] <- which(ID_d==i)}
for(i in 1:d){id_list[[i]] <- which(u_id==i)}
for(i in 1:length(unique(td1_d))){sent_list[[i]] <- which(td1_d==i)}
for(i in 1:length(unique(wd))){word_list[[i]] <- which(wd==i)}


####�}���R�t�A�������e�J�����@�ō���LDA���f���𐄒�####
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
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep

##���O���z�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z
alpha01 <- 1  
beta01 <- 0.5
gamma01 <- 0.5 

##�p�����[�^�̏����l
#tfidf�ŏ����l��ݒ�
tf <- WX/rowSums(WX)
idf1 <- log(nrow(WX)/colSums(WX > 0))
idf2 <- log(nrow(WX)/colSums(WX==0))

#�P��g�s�b�N�P�ʂ̃p�����[�^�̏����l
theta <- extraDistr::rdirichlet(d, rep(1, k1))   #���[�U�[�g�s�b�N�̏����l
phi <- extraDistr::rdirichlet(k1, idf1*10)   #�]���Ώی�̏o���m���̏����l

#��ʌ�g�s�b�N�P�ʂ̃p�����[�^�̏����l
omega <- 1/k2
gamma <- extraDistr::rdirichlet(k2, idf2*100)   #��ʌ�̏o���m���̏����l
y <- rbinom(f, 1, 0.5)
r <- rep(0.5, f)


##�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(d, k1, R/keep))
PHI <- array(0, dim=c(k1, v, R/keep))
OMEGA <- array(0, dim=c(d, k2, R/keep))
GAMMA <- array(0, dim=c(k2, v, R/keep))
SEG1 <- matrix(0, nrow=f, ncol=k1)
SEG2 <- matrix(0, nrow=a, ncol=k2)
Y <- matrix(0, nrow=f, ncol=2)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"

##MCMC����p�z��
wsum0 <- matrix(0, nrow=d, ncol=k1)
vf0 <- matrix(0, nrow=k1, ncol=v)
dsum0 <- matrix(0, nrow=d, ncol=k2)
sf0 <- matrix(0, nrow=k2, ncol=v)


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�g�s�b�N�֘A�ꂲ�ƂɃg�s�b�N���z�̃p�����[�^�𐄒�
  #�g�s�b�N�̏o�����Ɩޓx���v�Z
  word_par1 <- burden_fr(theta, phi, wd, words, k1)
  word_rate1 <- word_par1$Br   #�g�s�b�N�����m��
  LH1 <- rowSums(word_par1$Bur)   #�ޓx
  
  ##��ʌꂩ�ǂ������T���v�����O
  #��ʌ�̃g�s�b�N�̖ޓx�Ɗ����m�����v�Z
  word_par2 <- matrix(0, nrow=f, ncol=k2)
  for(j in 1:k2){
    word_par2[, j] <- omega * gamma[j, wd]
  }
  LH2 <- rowSums(word_par2)   #��ʌ�̖ޓx
  
  
  #�񍀕��z����ʌꂩ�ǂ������T���v�����O
  y_rate <- r*LH1 / (r*LH1 + (1-r)*LH2)   #���ݕϐ��̊����m��
  y <- rbinom(f, 1, y_rate)   #�񍀕��z�����ݕϐ����T���v�����O
  index <- which(y==1)

  #�������̍X�V
  r <- 0.5
  
  ##���������X�C�b�`���O�ϐ��Ɋ�Â��g�s�b�N�𐶐�
  #�P��P�ʂ̃g�s�b�N���T���v�����O
  Zi1 <- matrix(0, nrow=f, ncol=k1)
  zi1_vec <- rep(0, f)
  Zi1[index, ] <- rmnom(length(index), 1, word_rate1[index, ])
  zi1_vec[index] <- as.numeric(Zi1[index, ] %*% 1:k1)
  
  
  #���͒P�ʂ̃g�s�b�N���T���v�����O
  LH2 <- word_par2 * 10^10   #�ޓx�����������Ȃ��悤�ɒ萔��������
  LH2[index, ] <- 1
  LL <- matrix(0, nrow=a, ncol=k2)
  for(i in 1:a){
    LL[i, ] <- colProds(LH2[sent_list[[i]], ])
  }
  index_ones <- which(LL[, 1]==1)
  
  
  #�������z��蕶�͒P�ʂ̃g�s�b�N���T���v�����O
  sentence_rate <- LL / rowSums(LL)
  Zi2 <- rmnom(a, 1, sentence_rate)
  Zi2[index_ones, ] <- 0
  zi2_vec <- as.numeric(Zi2 %*% 1:k2)
  
  
  #�g�s�b�N�̍�����omega���T���v�����O
  #for(i in 1:d){
  #  dsum0[i, ] <- colSums(Zi2[sent_list[[i]], ])
  #}
  #dsum <- dsum0 + alpha01
  #omega <- extraDistr::rdirichlet(d, dsum)
  
  #�P�ꕪ�z���T���v�����O
  zi2_word <- rep(zi2_vec, w)[-index]
  wd0 <- wd[-index]
  
  for(j in 1:k2){
    index_seg <- which(zi2_word==j)
    sf0[j, ] <- plyr::count(c(wd0[index_seg], 1:v))$freq
  }
  gamma <- extraDistr::rdirichlet(k2, sf0)
  
  
  ##�P��P�ʂ̃p�����[�^���T���v�����O
  #�g�s�b�N���ztheta���T���v�����O
  for(i in 1:d){
    wsum0[i, ] <- colSums(Zi1[doc_list[[i]], ])
  }
  wsum <- wsum0 + alpha01
  theta <- extraDistr::rdirichlet(d, wsum)
  
  
  #�P�ꕪ�zphi���T���v�����O
  for(j in 1:v){
    vf0[, j] <- colSums(Zi1[word_list[[j]], ] )
  }
  vf <- vf0 + beta01
  phi <- extraDistr::rdirichlet(k1, vf)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , mkeep] <- phi
    #OMEGA[, , mkeep] <- omega
    GAMMA[, , mkeep] <- gamma
    
    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(mkeep >= burnin & rp%%keep==0){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
      Y <- Y + y
    }
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(c(mean(y), mean(y_vec)))
    print(round(cbind(theta[1:7, ], thetat[1:7, ]), 3))
    print(round(cbind(phi[, 246:255], phit[, 246:255]), 3))
    print(round(cbind(gamma[, 246:255], gammat[, 246:255]), 3))
  }
}
