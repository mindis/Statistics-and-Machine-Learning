#####�g�s�b�N�ǐՃ��f��####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(gtools)
library(bayesm)
library(extraDistr)
library(monomvn)
library(glmnet)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(54876)

####�f�[�^�̔���####
#set.seed(423943)
##�f�[�^�̐ݒ�
hh <- 300   #�ϑ��l��
pt <- 12   #�ϑ�����
d <- hh*pt   #��������
k <- 10   #�g�s�b�N��
v <- 200   #��b��
w <- rpois(d, rgamma(d, 100, 0.8))   #1����������̒P�ꐔ

##ID�̐ݒ�
id <- rep(1:hh, rep(pt, hh))
time <- rep(1:pt, hh)
ID <- data.frame(no=1:d, id, time)

##�p�����[�^�̐ݒ�
#�f�B�N�������O���z��ݒ�
alpha01 <- rep(0.4, k)   #1���ڂ̕����̃f�B�N�������O���z�̃p�����[�^
alpha02 <- rep(1, v)   #1���ڂ̒P��̃f�B�N�������O���z�̃p�����[�^
alpha11 <- rep(50, k)   #2���ڈȍ~�̕����̃f�B�N�������O���z�̃p�����[�^
alpha12 <- rep(1000, v)   #2���ڈȍ~�̒P��̃f�B�N�������O���z�̃p�����[�^


##���Ԃ��ƂɃf�B�N�������z���當���ƒP��̃g�s�b�N���z�𐶐�
#1���ڂ̃g�s�b�N���z�𐶐�
thetat <- theta <- matrix(0, nrow=d, ncol=k)
phit <- phi <- array(0, dim=c(k, v, pt))

index_time <- which(ID$time==1)
thetat[index_time, ] <- theta[index_time, ] <- extraDistr::rdirichlet(hh, alpha01)
phit[, , 1] <- phi[, , 1] <- extraDistr::rdirichlet(k, alpha02)

#2���ڈȍ~�̃g�s�b�N���z�𒀎�����
for(j in 2:pt){
  index_time <- which(ID$time==j)
  
  #�����g�s�b�N�𐶐�
  thetat[index_time, ] <- theta[index_time, ] <- 
    extraDistr::rdirichlet(hh, matrix(alpha11, nrow=hh, ncol=k, byrow=T) * thetat[index_time-1, ])
  thetat[index_time, ] <- theta[index_time, ]<- 
    (thetat[index_time, ] < 0.0001)*0.0001 + (thetat[index_time, ] >= 0.0001)*thetat[index_time, ]
  
  #�P��g�s�b�N�𐶐�
  phit[, , j] <- phi[, , j] <- 
    extraDistr::rdirichlet(k, matrix(alpha12, nrow=k, ncol=v, byrow=T) * phit[, , j-1])
  phit[, , j] <- phit[, , j]<- 
    (phit[, , j] < 0.00001)*0.00001 + (phit[, , j] >= 0.00001)*phit[, , j]
}

##�������z���當���f�[�^�𐶐�
WX <- matrix(0, nrow=d, ncol=v)
Z1 <- list()

for(i in 1:d){
  print(i)
  
  #�����̃g�s�b�N���z�𔭐�
  z1 <- extraDistr::rmnom(w[i], 1, theta[i, ])   #�����̃g�s�b�N���z�𔭐�
  
  #�����̃g�s�b�N���z����P��𔭐�
  zn <- z1 %*% c(1:k)   #0,1�𐔒l�ɒu��������
  zdn <- cbind(zn, z1)   #apply�֐��Ŏg����悤�ɍs��ɂ��Ă���
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi[x[1], , ID$time[i]])))   #�����̃g�s�b�N����P��𐶐�
  wdn <- colSums(wn)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  WX[i, ] <- wdn  
  
  #�����������g�s�b�N���i�[
  Z1[[i]] <- zdn[, 1]
}
colSums(WX)

#�f�[�^�s��𐮐��^�s��ɕύX
storage.mode(WX) <- "integer"


####�g�s�b�N���f���̂��߂̃f�[�^�Ɗ֐��̏���####
##���ꂼ��̕������̒P��̏o������ѕ⏕���̏o�����x�N�g���ɕ��ׂ�
##�f�[�^����pID���쐬
ID_list <- list()
hh_list <- list()
pt_list <- list()
wd_list <- list()

#���l���Ƃɋ��lID����ђP��ID���쐬
for(i in 1:nrow(WX)){
  print(i)
  
  #�P���ID�x�N�g�����쐬
  ID_list[[i]] <- rep(i, w[i])
  hh_list[[i]] <- rep(ID$id[i], w[i])
  pt_list[[i]] <- rep(ID$time[i], w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
  num2 <- subset(num1, num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
}

#���X�g���x�N�g���ɕϊ�
ID_d <- unlist(ID_list)
hh_d <- unlist(hh_list)
pt_d <- unlist(pt_list)
wd <- unlist(wd_list)

##�C���f�b�N�X���쐬
doc1_list <- list()
doc2_list <- list()
doc3_list <- list()
word_list <- list()

for(i in 1:length(unique(ID_d))) {doc1_list[[i]] <- which(ID_d==i)}
for(i in 1:length(unique(hh_d))) {doc2_list[[i]] <- which(hh_d==i)}
for(i in 1:length(unique(pt_d))) {doc3_list[[i]] <- which(pt_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- which(wd==i)}

#���Ԃ��ƂɃC���f�b�N�X���쐬
ptdoc_list <- list()
for(j in 1:pt){
  pd <- list()
  for(i in 1:hh){
    pd[[i]] <- subset(doc3_list[[j]], doc3_list[[j]] %in% doc2_list[[i]])
  }
  ptdoc_list[[j]] <- pd
}

ptwd_list <- list()
for(j in 1:pt){
  pd <- list()
  for(i in 1:length(unique(wd))){
    pd[[i]] <- subset(doc3_list[[j]], doc3_list[[j]] %in% word_list[[i]])
  }
  ptwd_list[[j]] <- pd
}
gc(); gc()


####�}���R�t�A�������e�J�����@�őΉ��g�s�b�N���f���𐄒�####
##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #���S�W���̊i�[�p
  for(kk in 1:k){
    #���S�W�����v�Z
    Bi <- rep(theta[, kk], w) * phi[kk, c(wd)]   #�ޓx
    Bur[, kk] <- Bi   
  }
  
  Br <- Bur / rowSums(Bur)   #���S���̌v�Z
  r <- colSums(Br) / sum(Br)   #�������̌v�Z
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}


##�A���S���Y���̐ݒ�
R <- 5000   #�T���v�����O��
keep <- 2   #2���1��̊����ŃT���v�����O���ʂ��i�[
iter <- 0
burnin <- 1000/keep

##���O���z�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z
alpha01 <- matrix(rep(1, k), nrow=hh, ncol=k, byrow=T)   
alpha02 <- matrix(0.5, nrow=k, ncol=v, byrow=T) 
alpha11 <- array(25, dim=c(hh, k, pt))
alpha12 <- array(100, dim=c(k, v, pt))

##�p�����[�^�̏����l
theta.ini <- extraDistr::rdirichlet(hh, rep(1, k))
phi.ini <- extraDistr::rdirichlet(k, rep(1, v))
theta <- matrix(0, nrow=d, ncol=k)
phi <- array(0, dim=c(k, v, pt))
for(j in 1:pt){
  theta[ID$time==j, ] <- theta.ini   #�����g�s�b�N�̃p�����[�^�̏����l
  phi[, , j] <- phi.ini   #�P��g�s�b�N�̃p�����[�^�̏����l
}

##�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(k, v, pt, R/keep))
W_SEG <- matrix(0, nrow=sum(w), ncol=k)
storage.mode(W_SEG) <- "integer"
gc(); gc()

#�C���f�b�N�X���쐬
index_id <- list()
index_pt <- list()
for(i in 1:hh){index_id[[i]] <- which(ID$id==i)}
for(i in 1:pt){index_pt[[i]] <- which(ID$time==i)}
vec <- 1/1:k

####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##1���ڂ̒P��g�s�b�N���T���v�����O
  #�P�ꂲ�ƂɃg�s�b�N�̏o�������v�Z
  word_rate <- burden_fr(theta[index_pt[[1]], ], phi[, , 1], wd[doc3_list[[1]]], w[index_pt[[1]]], k)$Br
  
  #�������z����P��g�s�b�N���T���v�����O
  Zi <- matrix(0, nrow=length(wd), ncol=k)
  word_cumsums <- rowCumsums(word_rate)
  rand <- matrix(runif(nrow(word_rate)), nrow=nrow(word_rate), ncol=k)   #��l����
  Zi0 <- ((k+1) - (word_cumsums > rand) %*% rep(1, k)) %*% vec   #�g�s�b�N���T���v�����O
  Zi0[Zi0!=1] <- 0
  Zi[doc3_list[[1]], ] <- Zi0
  
  ##�P��g�s�b�N�̃p�����[�^���X�V
  #�f�B�N�������z����theta���T���v�����O
  wsum0 <- matrix(0, nrow=hh, ncol=k)
  for(i in 1:hh){
    wsum0[i, ] <- colSums(Zi[ptdoc_list[[1]][[i]], ])
  }
  wsum <- wsum0 + alpha01 + 1   #�f�B�N�������z�̃p�����[�^
  theta[index_pt[[1]], ] <- extraDistr::rdirichlet(hh, wsum)   #�f�B�N�������z����g�s�b�N�������T���v�����O
  
  #�f�B�N�������z����phi���T���v�����O
  vf01 <- matrix(0, nrow=k, ncol=v)
  for(j in 1:v){
    vf01[, j] <- colSums(Zi[ptwd_list[[1]][[j]], ])
  }
  vf <- vf01 + alpha02 + 1   #�f�B�N�������z�̃p�����[�^
  phi[, , 1] <- extraDistr::rdirichlet(k, vf)   #�f�B�N�������z����phi���T���v�����O
  
  
  ##2���ڈȍ~�̒P��g�s�b�N�𒀎��T���v�����O
  #�P�ꂲ�ƂɃg�s�b�N�̏o�������v�Z
  wsum0 <- array(0, dim=c(hh, k, pt))
  vf01 <- array(0, dim=c(k, v, pt))
  
  for(p in 2:pt){

    #�P�ꂲ�ƂɃg�s�b�N�̏o�������v�Z
    word_rate <- burden_fr(theta[index_pt[[p]], ], phi[, , p], wd[doc3_list[[p]]], w[index_pt[[p]]], k)$Br
    
    #�������z����P��g�s�b�N���T���v�����O
    word_cumsums <- rowCumsums(word_rate)
    rand <- matrix(runif(nrow(word_rate)), nrow=nrow(word_rate), ncol=k)   #��l����
    Zi0 <- ((k+1) - (word_cumsums > rand) %*% rep(1, k)) %*% vec   #�g�s�b�N���T���v�����O
    Zi0[Zi0!=1] <- 0
    Zi[doc3_list[[p]], ] <- Zi0
    
    ##�P��g�s�b�N�̃p�����[�^���X�V
    #�f�B�N�������z����theta���T���v�����O
    index_now1 <- index_pt[[p]]
    index_obs1 <- index_pt[[p-1]]
    for(i in 1:hh){
      wsum0[i, , p] <- colSums(Zi[ptdoc_list[[p]][[i]], ])
    }
    wsum <- wsum0[, , p] + alpha11[, , p]*theta[index_obs1, ] + 1   #�f�B�N�������z�̃p�����[�^
    theta[index_now1, ] <- extraDistr::rdirichlet(hh, wsum)   #�f�B�N�������z����g�s�b�N�������T���v�����O
    
    #�f�B�N�������z����phi���T���v�����O
    for(j in 1:v){
      vf01[, j, p] <- colSums(Zi[ptwd_list[[p]][[j]], ])
    }
    vf <- vf01[, , p] + alpha12[, , p]*phi[, , p-1] + 1   #�f�B�N�������z�̃p�����[�^
    phi[, , p] <- extraDistr::rdirichlet(k, vf)   #�f�B�N�������z����phi���T���v�����O
    
    ##�n�C�p�[�p�����[�^���X�V
    ##�����g�s�b�N���z�̃n�C�p�[�p�����[�^alpha���X�V
    par <- alpha11[, , p]*theta[index_obs1, ]
    alpha_m <- rowSums(theta[index_obs1, ] * (digamma(wsum0[, , p] + par) - digamma(par)))
    alpha_d <- digamma(w[index_now1] + alpha11[, 1, p]) - digamma(alpha11[, 1, p])
    alpha11[, , p] <- matrix(alpha11[, 1, p] * alpha_m / alpha_d, nrow=hh, ncol=k)
    
    #�P��g�s�b�N���z�̃n�C�p�[�p�����[�^beta���X�V
    #par <- alpha12[, , p]*phi[, , p-1]
    #beta_m <- rowSums(phi[, , p-1] * (digamma(vf01[, , p] + par) - digamma(par)))
    #beta_d <- digamma(rowSums(vf01[, , p]) + alpha12[, 1, p]) - digamma(alpha12[, 1, p])
    #alpha12[, , p] <- matrix(alpha12[, 1, p] * beta_m / beta_d, nrow=k, ncol=v)
  }
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI[, , , mkeep] <- phi
    if(rp>=burnin){
      W_SEG <- W_SEG + Zi
    }
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(round(cbind(theta[13:24, ], thetat[13:24, ]), 2)) 
    print(round(alpha11[1:10, 1, 2], 3))
    print(round(alpha12[, 1, 10], 3))
    #print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
  }
}

digamma(20)


matplot(t(THETA[10, 4:5, ]), type="l")