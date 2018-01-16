#####�x�C�W�A���L�������}���R�t���ڃ��f��#####
library(MASS)
library(Matrix)
library(flexmix)
library(mclust)
library(matrixStats)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(90345)

####�f�[�^�̔���####
hh <- 5000   #���[�U�[��
S <- 20   #�y�[�W��
k <- 10   #������
seg <- as.numeric(rmnom(hh, 1, rep(1, k)) %*% 1:k)
rt <- as.numeric(table(seg)/hh)

##�p�����[�^�̐ݒ�
#�}���R�t���ڍs��̐ݒ�
Pr <- array(0, dim=c(S-1, S, k))
for(i in 1:k){
  for(j in 1:(S-1)){
    if(j==1){
      Pr[j, -S, i] <- extraDistr::rdirichlet(1, rep(1, S-1))   
    } else {
      Pr[j, , i] <- extraDistr::rdirichlet(1, rep(1, S))
    }
  }
}

##���[�U�[���ƂɃR���o�[�W��������܂Ńf�[�^�𒀎�����
Data_list <- list()
id_list <- list()

for(i in 1:hh){
  data <- matrix(0, nrow=1000, ncol=S)
  data[1, ] <- rmnom(1, 1, Pr[1, , seg[i]])   #1�A�N�Z�X�ڂ̃��O�𐶐�
  
  for(j in 2:1000){
    data[j, ] <- rmnom(1, 1, Pr[which(data[j-1, ]==1), , seg[i]])   #2�A�N�Z�X�ȍ~�̃��O�𐶐�
    if(data[j, S]==1) break   #�R���o�[�W�������Ă�����break
  }
  Data_list[[i]] <- data[rowSums(data) > 0, ]
  id_list[[i]] <- rep(i, sum(rowSums(data) > 0))
}
Data <- do.call(rbind, Data_list)
id <- unlist(id_list)
as.matrix(data.frame(id, Data) %>%
            dplyr::group_by(id) %>%
            dplyr::summarize_all(funs(sum)))
r <- rep(1/k, k)


####�}���R�t�A�������e�J�����@�ŗL�������}���R�t���ڃ��f���𐄒�####
##���ڃx�N�g�����쐬
index_list <- list()
for(i in 1:hh){
  data <- Data[id==i, ]
  index <- rep(0, nrow(data))
  index[1] <- 1
  index[-1] <- data[1:(nrow(data)-1), ] %*% 1:S
  index_list[[i]] <- index
}
index_trans <- unlist(index_list)

##���ݕϐ�z���v�Z���邽�߂̊֐�
LLobz <- function(theta, r, Data, id, index, hh, k){
  index <- index_trans
  #���ݕϐ����Ƃ̖ޓx���v�Z
  LLind0 <- matrix(0, nrow=hh, ncol=k)
  log_theta <- log(theta)   #�p�����[�^��ΐ��ϊ�
  for(j in 1:k){
    Li <- rowSums(Data * log_theta[index, , j])
    LLind0[, j] <- tapply(Li, id, sum)
  }
  LLind <- exp(LLind0 - apply(LLind0, 1, max))
  
  #���ݕϐ��̊����m�����v�Z
  LLho <- matrix(r, nrow=hh, ncol=k, byrow=T) * LLind   #�ϑ��f�[�^�̖ޓx
  z <- LLho / matrix(rowSums(LLho), nrow=hh, ncol=k)   #���ݕϐ�z�̊����m��
  rval <- list(z=z)
  return(rval)
}

##�A���S���Y���̐ݒ�
R <- 2000
keep <- 2  
iter <- 0
burnin <- 200/keep0
disp <- 10

##���O���z�̐ݒ�
alpha <- 1
beta <- 1

##�����l�̐ݒ�
#�p�����[�^�̏����l
theta <- array(0, dim=c(S-1, S, k))
for(i in 1:k){
  for(j in 1:(S-1)){
    if(j==1){
      theta[j, -S, i] <- extraDistr::rdirichlet(1, rep(5, S-1))   
    } else {
      theta[j, , i] <- extraDistr::rdirichlet(1, rep(5, S))
    }
  }
}
theta <- theta + 0.0001

#�������̏����l
r <- rep(1/k, k)

##�T���v�����O���ʂ̕ۑ��p�z��
THETA <- array(0, dim=c(S-1, S, k, R/keep))
RATE <- matrix(0, nrow=R/keep, ncol=k)
SEG <- matrix(0, nrow=hh, ncol=k)
storage.mode(Z) <- "integer"


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){

  ##���[�U�[���ƂɃZ�O�����g���T���v�����O
  z <- LLobz(theta, r, Data, id, index_trans, hh, k)$z   #�Z�O�����g�����m��
  Zi <- rmnom(hh, 1, z)   #�������z����Z�O�����g�𐶐�

  ##�p�����[�^���T���v�����O
  #���������T���v�����O
  wsum <- colSums(Zi) + alpha
  r <- as.numeric(extraDistr::rdirichlet(1, wsum))

  #�}���R�t���ڍs����T���v�����O
  theta <- array(0, dim=c(S-1, S, k))
  for(j in 1:k){
    index <- which(Zi[id, j]==1)
    data <- Data[index, ]
    theta0 <- as.matrix(data.frame(id=index_trans[index], data=data) %>%
                          dplyr::group_by(id) %>%
                          dplyr::summarise_all(funs(sum)))[, 2:(S+1)]
    dsum <- theta0 + beta
    theta[, , j] <- extraDistr::rdirichlet(S-1, dsum)
  }
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA[, , , mkeep] <- theta
    RATE[mkeep, ] <- r

    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(mkeep >= burnin & rp%%keep==0){
      SEG <- SEG + Zi
    }
    
    #�T���v�����O���ʂ��m�F
    if(rp%%disp==0){
      print(rp)
      print(round(rbind(r, rt), 3))
    }
  }
}

####���茋�ʂ̊m�F####
burnin <- 200/keep0
RS <- R/keep

##�T���v�����O���ʂ̉���
matplot(t(THETA[1, , 1, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA[2, , 2, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA[3, , 3, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA[4, , 4, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA[5, , 5, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(RATE, type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

#�p�����[�^�̐���l
theta_mu <- array(0, dim=c(S-1, S, k))
for(j in 1:k){
  theta_mu[, , j] <- apply(THETA[, , j, burnin:RS], c(1, 2), mean)
}

round(theta_mu, 3)   #�}���R�t���ڊm���̐���l
round(Pr, 3)   #�}���R�t���ڊm���̐^�l
round(rbind(r1=colMeans(RATE[burnin:RS, ]), rt), 3)   #������


#���ݕϐ��̊���
round(Z <- SEG / rowSums(SEG))   #���ݕϐ�z�̊����m��
round(cbind(seg0=seg, seg1=apply(Z, 1, which.max), Z), 3)   #���肳�ꂽ���ݕϐ��Ɛ^�̐��ݕϐ�

##�K���x
#���j�O�����̑ΐ��ޓx
LLst1 <- sum(Data %*% log(colMeans(Data)))

#�}���R�t���f���̑ΐ��ޓx
par <- matrix(0, nrow=S-1, ncol=S)
for(j in 1:max(index_trans)){
  index <- which(index_trans==j)
  par[j, ] <- colMeans(Data[index, ])
}
log_par <- log(par)
LLi <- Data * log_par[index_trans, ]
LLi[is.nan(LLi)] <- 0
LLst2 <- sum(LLi)

#�����}���R�t���f���̑ΐ��ޓx
LLi <- matrix(0, nrow=nrow(Data), ncol=k)
for(j in 1:k){
  log_theta <- log(theta_mu[, , j])
  LLi[, j] <- rowSums(Z[id, j] * Data * log_theta[index_trans, ])
}
LL <- sum(LLi)

#�K���x�̔�r
round(LLc <- c(LLst1, LLst2, LL), 3)   #�ΐ��ޓx
round(exp(-LLc / nrow(Data)), 3)   #Perplexity