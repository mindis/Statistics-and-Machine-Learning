#####Switching Multinomial LDA model#####
options(warn=0)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
r <- 5   #�]���X�R�A��
s <- 3   #�ɐ��l��
a <- 3   #����
k11 <- 10   #���[�U�[�̃e�L�X�g�̃g�s�b�N��
k12 <- 15   #�A�C�e���̃e�L�X�g�̃g�s�b�N��
hh <- 1000   #���r���A�[��
item <- 200   #�A�C�e����
v1 <- 300   #�]���X�R�A�̌�b��
v2 <- 350   #���[�U�[�g�s�b�N�̌�b��
v3 <- 350   #�A�C�e���g�s�b�N�̌�b��
v <- v1 + v2 + v3   #����b��
spl <- matrix(1:v1, nrow=s, ncol=v1/s, byrow=T)
v1_index <- 1:v1
v2_index <- (v1+1):v2
v3_index <- (v2+1):v

##ID�ƌ����x�N�g���̍쐬
#ID�����ݒ�
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)

#�����x�N�g�����쐬
for(rp in 1:100){
  m_vec <- rep(0, hh*item)
  for(i in 1:item){
    prob <- runif(1, 0.025, 0.16)
    m_vec[item_id0==i] <- rbinom(hh, 1, prob)
  }
  m_index <- which(m_vec==1)
  
  #���S��ID��ݒ�
  user_id <- user_id0[m_index]
  item_id <- item_id0[m_index]
  d <- length(user_id)   #�����r���[��
  
  #���ׂẴp�^�[��������������break
  if(length(unique(user_id))==hh & length(unique(item_id))==item) break
}

#�P�ꐔ��ݒ�
w <- rpois(d, rgamma(d, 25, 0.5))   #����������̒P�ꐔ
f <- sum(w)   #���P�ꐔ
n_user <- plyr::count(user_id)$freq
n_item <- plyr::count(item_id)$freq

#�P��ID��ݒ�
u_id <- rep(user_id, w)
i_id <- rep(item_id, w)
d_id <- rep(1:d, w)

#�C���f�b�N�X��ݒ�
user_index <- list()
item_index <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
}

##�p�����[�^�̐ݒ�
#�f�B���N�����z�̎��O���z�̐ݒ�
alpha11 <- rep(0.15, k11)
alpha12 <- rep(0.15, k12)
alpha21 <- c(rep(0.5, v1/s), rep(0.025, v1/s), rep(0.0025, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha22 <- c(rep(0.3, v1/s), rep(0.1, v1/s), rep(0.025, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha23 <- c(rep(0.2, v1/s), rep(1.0, v1/s), rep(0.2, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha24 <- c(rep(0.025, v1/s), rep(0.1, v1/s), rep(0.3, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha25 <- c(rep(0.0025, v1/s), rep(0.025, v1/s), rep(0.5, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha2 <- rbind(alpha21, alpha22, alpha23, alpha24, alpha25)
alpha31 <- c(rep(0.001, v1/s), rep(0.001, v1/s), rep(0.001, v1/s), rep(0.1, v2), rep(0.002, v3))
alpha32 <- c(rep(0.001, v1/s), rep(0.001, v1/s), rep(0.001, v1/s), rep(0.002, v2), rep(0.1, v3))
beta1 <- c(1.6, 4.8, 5.6)

##���ׂĂ̒P�ꂪ�o������܂Ńf�[�^�̐����𑱂���
for(rp in 1:1000){
  print(rp)
  
  #���O���z����p�����[�^�𐶐�
  theta11 <- thetat11 <- extraDistr::rdirichlet(hh, alpha11)
  theta12 <- thetat12 <- extraDistr::rdirichlet(item, alpha12)
  omega <- omegat <- extraDistr::rdirichlet(r, alpha2)
  phi <- phit <- extraDistr::rdirichlet(k11, alpha31)
  gamma <- gammat <- extraDistr::rdirichlet(k12, alpha32)
  lambda <- lambdat <- extraDistr::rdirichlet(hh, beta1)
  
  ##���f���Ɋ�Â��f�[�^�𐶐�
  WX <- matrix(0, nrow=d, ncol=v)
  y <- rep(0, d)
  Z1_list <- Z21_list <- Z22_list <- wd_list <- list()
  
  for(i in 1:d){
    #���[�U�[�ƃA�C�e���𒊏o
    u_index <- user_id[i]
    i_index <- item_id[i]
    
    #�]���X�R�A�𐶐�
    y[i] <- as.numeric(rmnom(1, 1, c(0.1, 0.225, 0.3, 0.25, 0.125)) %*% 1:r)

    #�������z����X�C�b�`���O�ϐ��𐶐�
    z1 <- rmnom(w[i], 1, lambda[u_index, ])
    z1_vec <- as.numeric(z1 %*% 1:a)
    index_z11 <- which(z1[, 1]==1)
    
    #���[�U�[�g�s�b�N�𐶐�
    z21 <- matrix(0, nrow=w[i], ncol=k11)
    index_z21 <- which(z1[, 2]==1)
    if(sum(z1[, 2]) > 0){
      z21[index_z21, ] <- rmnom(sum(z1[, 2]), 1, theta11[u_index, ])
    }
    z21_vec <- as.numeric(z21 %*% 1:k11)
    
    #�A�C�e���g�s�b�N�𐶐�
    z22 <- matrix(0, nrow=w[i], ncol=k12)
    index_z22 <- which(z1[, 3]==1)
    if(sum(z1[, 3]) > 0){
      z22[index_z22, ] <- rmnom(sum(z1[, 3]), 1, theta12[i_index, ])
    }
    z22_vec <- as.numeric(z22 %*% 1:k12)
    
    #�g�s�b�N����P��𐶐�
    words <- matrix(0, nrow=w[i], ncol=v)
    if(sum(z1[, 1]) > 0){
      words[index_z11, ] <- rmnom(sum(z1[, 1]), 1, omega[y[i], ])
    }
    if(sum(z1[, 2]) > 0){
      words[index_z21, ] <- rmnom(sum(z1[, 2]), 1, phi[z21_vec[index_z21], ])
    }
    if(sum(z1[, 3]) > 0){
      words[index_z22, ] <- rmnom(sum(z1[, 3]), 1, gamma[z22_vec[index_z22], ])
    }
    word_vec <- as.numeric(words %*% 1:v)
    WX[i, ] <- colSums(words)
    
    #�f�[�^���i�[
    wd_list[[i]] <- word_vec
    Z1_list[[i]] <- z1
    Z21_list[[i]] <- z21
    Z22_list[[i]] <- z22
  }
  if(min(colSums(WX)) > 0) break
}

#���X�g��ϊ�
wd <- unlist(wd_list)
Z1 <- do.call(rbind, Z1_list)
Z21 <- do.call(rbind, Z21_list)
Z22 <- do.call(rbind, Z22_list)
storage.mode(Z1) <- "integer"
storage.mode(Z21) <- "integer"
storage.mode(Z22) <- "integer"
storage.mode(WX) <- "integer"


####�}���R�t�A�������e�J�����@��Switching Binomial LDA�𐄒�####
##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, w, k){
  #���S�W�����v�Z
  Bur <- theta[w, ] * t(phi)[wd, ]   #�ޓx
  Br <- Bur / rowSums(Bur)   #���S��
  r <- colSums(Br) / sum(Br)   #������
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##�A���S���Y���̐ݒ�
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##���O���z�̐ݒ�
alpha11 <- 0.1; alpha12 <- 0.1
alpha21 <- 0.1; alpha22 <- 0.1; alpha23 <- 0.1
beta <- 0.5

##�p�����[�^�̐^�l
theta11 <- thetat11
theta12 <- thetat12
phi <- phit
gamma <- gammat
omega <- omegat
lambda <- lambdat

##�p�����[�^�̏����l��ݒ�
#�g�s�b�N���z�̏����l
theta11 <- extraDistr::rdirichlet(hh, rep(1.0, k11))
theta12 <- extraDistr::rdirichlet(item, rep(1.0, k12))

#�P�ꕪ�z�̏����l
phi <- extraDistr::rdirichlet(k11, rep(1.0, v))   #���[�U�[�̒P�ꕪ�z�̏����l
gamma <- extraDistr::rdirichlet(k12, rep(1.0, v))   #�A�C�e���̒P�ꕪ�z�̏����l
omega <- extraDistr::rdirichlet(r, rep(5.0, v))   #�]���X�R�A�̒P�ꕪ�z�̏����l

#�X�C�b�`���O�ϐ��̏����l
lambda <- matrix(1/s, nrow=hh, ncol=s)


##�p�����[�^�̊i�[�p�z��
THETA11 <- array(0, dim=c(hh, k11, R/keep))
THETA12 <- array(0, dim=c(item, k12, R/keep))
PHI <- array(0, dim=c(k11, v, R/keep))
GAMMA <- array(0, dim=c(k12, v, R/keep))
OMEGA <- array(0, dim=c(r, v, R/keep))
LAMBDA <- array(0, dim=c(hh, s, R/keep))
SEG1 <- matrix(0, nrow=f, ncol=a)
SEG21 <- matrix(0, nrow=f, ncol=k11)
SEG22 <- matrix(0, nrow=f, ncol=k12)
storage.mode(SEG21) <- "integer"
storage.mode(SEG22) <- "integer"

##�f�[�^�ƃC���f�b�N�X�̐ݒ�
#�C���f�b�N�X�̐ݒ�
user_list <- user_vec <- list()
item_list <- item_vec <- list()
wd_list <- wd_vec <- list()
user_n <- rep(0, hh)
item_n <- rep(0, item)
for(i in 1:hh){
  user_list[[i]] <- which(u_id==i)
  user_vec[[i]] <- rep(1, length(user_list[[i]]))
  user_n[i] <- sum(user_vec[[i]])
}
for(i in 1:item){
  item_list[[i]] <- which(i_id==i)
  item_vec[[i]] <- rep(1, length(item_list[[i]]))
  item_n[i] <- sum(item_vec[[i]])
}
for(j in 1:v){
  wd_list[[j]] <- which(wd==j)
  wd_vec[[j]] <- rep(1, length(wd_list[[j]]))
}

#�f�[�^�̐ݒ�
y_vec <- y[d_id]
y_data <- matrix(as.numeric(table(1:f, y_vec)), nrow=f, ncol=r)
storage.mode(y_data) <- "integer"
r_vec <- rep(1, r)
a_vec <- rep(1, a)
vec11 <- rep(1, k11)
vec12 <- rep(1, k12)

##�ΐ��ޓx�̊�l
par <- colSums(WX) / f
LLst <- sum(WX %*% log(par))


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){

  ##�������z���X�C�b�`���O�ϐ����T���v�����O
  #�]���X�R�A�A���[�U�[����уA�C�e���̊��Җޓx��ݒ�
  Li_score <- as.numeric((t(omega)[wd, ] * y_data) %*% r_vec)   #�X�R�A�ޓx
  Li_user <- theta11[u_id, ] * t(phi)[wd, ]   #���[�U�[�ޓx
  par_user <- as.numeric(Li_user %*% vec11)   #���[�U�[�̊��Җޓx
  Li_item <- theta12[i_id, ] * t(gamma)[wd, ]   #�A�C�e���ޓx
  par_item <- as.numeric(Li_item %*% vec12)   #�A�C�e���̊��Җޓx
  par <- cbind(Li_score, par_user, par_item)
  
  #���݊m������X�C�b�`���O�ϐ��𐶐�
  lambda_r <- lambda[u_id, ]   #�X�C�b�`���O�ϐ��̎��O���z
  par_r <- lambda_r * par
  s_prob <- par_r / as.numeric(par_r %*% a_vec)   #�X�C�b�`���O�ϐ��̊����m��
  Zi1 <- rmnom(f, 1, s_prob)   #�������z����X�C�b�`���O�ϐ��𐶐�
  Zi1_T <- t(Zi1)
  index_z21 <- which(Zi1[, 2]==1)
  index_z22 <- which(Zi1[, 3]==1)
  
  #�f�B���N�����z���獬�������T���v�����O
  rsum0 <- matrix(0, nrow=hh, ncol=a)
  for(i in 1:hh){
    rsum0[i, ] <- Zi1_T[, user_list[[i]]] %*% user_vec[[i]]
  }
  rsum <- rsum0 + beta   #�f�B���N�����z�̃p�����[�^
  lambda <- extraDistr::rdirichlet(hh, rsum)   #�f�B���N�����z����lambda���T���v�����O
  
  
  ##���[�U�[����уA�C�e���̃g�s�b�N���T���v�����O
  #�g�s�b�N�̊����m���𐄒�
  z_rate1 <- Li_user[index_z21, ] / par_user[index_z21]   #���[�U�[�̃g�s�b�N�����m��
  z_rate2 <- Li_item[index_z22, ] / par_item[index_z22]   #�A�C�e���̃g�s�b�N�����m��
  
  #�������z����g�s�b�N�𐶐�
  Zi21 <- matrix(0, nrow=f, ncol=k11)
  Zi22 <- matrix(0, nrow=f, ncol=k12)
  Zi21[index_z21, ] <- rmnom(nrow(z_rate1), 1, z_rate1)
  Zi22[index_z22, ] <- rmnom(nrow(z_rate2), 1, z_rate2)
  Zi21_T <- t(Zi21)
  Zi22_T <- t(Zi22)
  
  
  ##�g�s�b�N���f���̃p�����[�^���T���v�����O
  #���[�U�[�̃g�s�b�N���z���T���v�����O
  wusum0 <- matrix(0, nrow=hh, ncol=k11)
  for(i in 1:hh){
    wusum0[i, ] <- Zi21_T[, user_list[[i]]] %*% user_vec[[i]]
  }
  wusum <- wusum0 + alpha21   #�f�B���N�����z�̃p�����[�^
  theta11 <- extraDistr::rdirichlet(hh, wusum)   #�f�B���N�����z����theta21���T���v�����O
  
  #�A�C�e���̃g�s�b�N���z���T���v�����O
  wisum0 <- matrix(0, nrow=item, ncol=k12)
  for(i in 1:item){
    wisum0[i, ] <- Zi22_T[, item_list[[i]]] %*% item_vec[[i]]
  }
  wisum <- wisum0 + alpha22   #�f�B���N�����z�̃p�����[�^
  theta12 <- extraDistr::rdirichlet(item, wisum)   #�f�B���N�����z����theta22���T���v�����O
  
  
  ##�]���X�R�A�A���[�U�[����уA�C�e���̒P�ꕪ�z���T���v�����O
  y_data_t <- t(y_data * Zi1[, 1])
  vssum0 <- matrix(0, nrow=r, ncol=v)
  vusum0 <- matrix(0, nrow=k11, ncol=v)
  visum0 <- matrix(0, nrow=k12, ncol=v)
  for(j in 1:v){
    vssum0[, j] <- y_data_t[, wd_list[[j]], drop=FALSE] %*% wd_vec[[j]]
    vusum0[, j] <- Zi21_T[, wd_list[[j]], drop=FALSE] %*% wd_vec[[j]]
    visum0[, j] <- Zi22_T[, wd_list[[j]], drop=FALSE] %*% wd_vec[[j]]
  }
  vssum <- vssum0 + alpha21; vusum <- vusum0 + alpha22; visum <- visum0 + alpha23
  omega <- extraDistr::rdirichlet(r, vssum)
  phi <- extraDistr::rdirichlet(k11, vusum)
  gamma <- extraDistr::rdirichlet(k12, visum)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA11[, , mkeep] <- theta11
    THETA12[, , mkeep] <- theta12
    PHI[, , mkeep] <- phi
    GAMMA[, , mkeep] <- gamma
    OMEGA[, , mkeep] <- omega
    LAMBDA[, , mkeep] <- lambda
  }  
  
  #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
  if(rp%%keep==0 & rp >= burnin){
    SEG1 <- SEG1 + Zi1
    SEG21 <- SEG21 + Zi21
    SEG22 <- SEG22 + Zi22
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    LL <- sum(log(rowSums(Zi1 * par)))
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(c(LL, LLst))
    print(round(c(colMeans(Zi1), colMeans(Z1)), 3))
    print(round(cbind(phi[, (v1-5):(v1+4)], phit[, (v1-5):(v1+4)]), 3))
  }
}

####�T���v�����O���ʂ̉����Ɨv��####
burnin <- 1000/keep
RS <- R/keep

##�T���v�����O���ʂ̃v���b�g
#�X�C�b�`���O�ϐ��̍������̉���
matplot(t(LAMBDA[, 1, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(LAMBDA[, 2, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(LAMBDA[, 3, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

#�g�s�b�N���z�̃T���v�����O���ʂ��v���b�g
matplot(t(THETA11[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA11[100, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA11[250, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA11[500, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA11[1000, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA12[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA12[50, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA12[100, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA12[150, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA12[200, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

#�P�ꕪ�z�̃T���v�����O���ʂ̉���
matplot(t(PHI[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(PHI[3, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(PHI[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(PHI[7, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(PHI[9, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(GAMMA[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(GAMMA[4, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(GAMMA[8, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(GAMMA[12, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(GAMMA[15, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(OMEGA[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(OMEGA[2, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(OMEGA[3, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(OMEGA[4, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(OMEGA[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")


##�T���v�����O���ʂ̗v��
#�T���v�����O���ʂ̎��㕽��
round(cbind(apply(LAMBDA[, , burnin:RS], c(1, 2), mean), lambdat), 3)   #�X�C�b�`���O�ϐ��̍������̎��㕽��
round(cbind(apply(THETA11[, , burnin:RS], c(1, 2), mean), thetat11), 3)   #���[�U�[�̃g�s�b�N�����̎��㕽��
round(cbind(apply(THETA12[, , burnin:RS], c(1, 2), mean), thetat12), 3)   #�A�C�e���̃g�s�b�N�����̎��㕽��
round(cbind(t(apply(PHI[, , burnin:RS], c(1, 2), mean)), t(phit)), 3)   #���[�U�[�̒P�ꕪ�z�̎��㕽��
round(cbind(t(apply(GAMMA[, , burnin:RS], c(1, 2), mean)), t(gammat)), 3)   #�A�C�e���̒P�ꕪ�z�̎��㕽��
round(cbind(t(apply(OMEGA[, , burnin:RS], c(1, 2), mean)), t(omegat)), 3)   #�]���X�R�A�̒P�ꕪ�z�̎��㕽��


#�T���v�����O���ʂ̎���M�p���
round(apply(LAMBDA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.025)), 3)
round(apply(LAMBDA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.975)), 3)
round(apply(THETA1[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.025)), 3)
round(apply(THETA2[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.975)), 3)
round(t(apply(PHI[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.025))), 3)
round(t(apply(PHI[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.975))), 3)
round(t(apply(GAMMA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.025))), 3)
round(t(apply(GAMMA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.975))), 3)
round(t(apply(OMEGA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.025))), 3)
round(t(apply(OMEGA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.975))), 3)


##�T���v�����O���ꂽ���ݕϐ��̗v��
n <- max(SEG1)
round(cbind(SEG1/n, Z1), 3)
round(cbind(rowSums(SEG21), SEG21/max(rowSums(SEG21)), Z21 %*% 1:k11), 3)
round(cbind(rowSums(SEG22), SEG22/max(rowSums(SEG22)), Z22 %*% 1:k12), 3)
