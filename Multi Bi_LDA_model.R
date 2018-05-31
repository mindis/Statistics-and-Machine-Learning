#####Multi Bi-LDA model#####
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
k11 <- 5   #���[�U�[�̕]���X�R�A�̃g�s�b�N��
k12 <- 5   #�A�C�e���̕]���X�R�A�̃g�s�b�N��
K1 <- matrix(1:(k11*k12), nrow=k11, ncol=k12, byrow=T)   #�g�s�b�N�̔z��
k21 <- 10   #���[�U�[�̃e�L�X�g�̃g�s�b�N��
k22 <- 15   #�A�C�e���̃e�L�X�g�̃g�s�b�N��
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
user_index <- user_ones <- list()
item_index <- item_ones <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
  user_ones[[i]] <- rep(1, length(user_index[[i]]))
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
  item_ones[[j]] <- rep(1, length(item_index[[j]]))
}

##�p�����[�^�̐ݒ�
#�f�B���N�����z�̎��O���z�̐ݒ�
alpha11 <- rep(0.2, k11)
alpha12 <- rep(0.2, k12)
alpha21 <- rep(0.15, k21)
alpha22 <- rep(0.15, k22)
alpha3 <- c(0.1, 0.225, 0.3, 0.25, 0.125) * r
alpha41 <- c(rep(0.5, v1/s), rep(0.025, v1/s), rep(0.0025, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha42 <- c(rep(0.3, v1/s), rep(0.1, v1/s), rep(0.025, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha43 <- c(rep(0.2, v1/s), rep(1.0, v1/s), rep(0.2, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha44 <- c(rep(0.025, v1/s), rep(0.1, v1/s), rep(0.3, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha45 <- c(rep(0.0025, v1/s), rep(0.025, v1/s), rep(0.5, v1/s), rep(0.001, v2), rep(0.001, v3))
alpha4 <- rbind(alpha41, alpha42, alpha43, alpha44, alpha45)
alpha51 <- c(rep(0.001, v1/s), rep(0.001, v1/s), rep(0.001, v1/s), rep(0.1, v2), rep(0.002, v3))
alpha52 <- c(rep(0.001, v1/s), rep(0.001, v1/s), rep(0.001, v1/s), rep(0.002, v2), rep(0.1, v3))
beta1 <- c(1.6, 4.8, 5.6)

##���ׂĂ̒P�ꂪ�o������܂Ńf�[�^�̐����𑱂���
for(rp in 1:1000){
  print(rp)
  
  #���O���z����p�����[�^�𐶐�
  theta11 <- thetat11 <- extraDistr::rdirichlet(hh, alpha11)
  theta12 <- thetat12 <- extraDistr::rdirichlet(item, alpha12)
  theta21 <- thetat21 <- extraDistr::rdirichlet(hh, alpha21)
  theta22 <- thetat22 <- extraDistr::rdirichlet(item, alpha22)
  eta <- etat <- extraDistr::rdirichlet(k11*k12, alpha3)
  omega <- omegat <- extraDistr::rdirichlet(r, alpha4)
  phi <- phit <- extraDistr::rdirichlet(k21, alpha51)
  gamma <- gammat <- extraDistr::rdirichlet(k22, alpha52)
  lambda <- lambdat <- extraDistr::rdirichlet(hh, beta1)
  

  ##���f���Ɋ�Â��f�[�^�𐶐�
  WX <- matrix(0, nrow=d, ncol=v)
  y <- rep(0, d)
  U1 <- matrix(0, nrow=d, ncol=k11)
  U2 <- matrix(0, nrow=d, ncol=k12)
  Z1_list <- Z21_list <- Z22_list <- wd_list <- list()
  
  for(i in 1:d){
    #���[�U�[�ƃA�C�e���𒊏o
    u_index <- user_id[i]
    i_index <- item_id[i]
    
    #�]���X�R�A�̃g�s�b�N�𐶐�
    u1 <- as.numeric(rmnom(1, 1, theta11[u_index, ]))
    u2 <- as.numeric(rmnom(1, 1, theta12[i_index, ]))
    
    #�]���X�R�A�̃g�s�b�N����X�R�A�𐶐�
    y[i] <- as.numeric(rmnom(1, 1, eta[K1[which.max(u1), which.max(u2)], ]) %*% 1:r)
    K1
    #�������z����X�C�b�`���O�ϐ��𐶐�
    z1 <- rmnom(w[i], 1, lambda[u_index, ])
    z1_vec <- as.numeric(z1 %*% 1:a)
    index_z11 <- which(z1[, 1]==1)
    
    #���[�U�[�g�s�b�N�𐶐�
    z21 <- matrix(0, nrow=w[i], ncol=k21)
    index_z21 <- which(z1[, 2]==1)
    if(sum(z1[, 2]) > 0){
      z21[index_z21, ] <- rmnom(sum(z1[, 2]), 1, theta21[u_index, ])
    }
    z21_vec <- as.numeric(z21 %*% 1:k21)
    
    #�A�C�e���g�s�b�N�𐶐�
    z22 <- matrix(0, nrow=w[i], ncol=k22)
    index_z22 <- which(z1[, 3]==1)
    if(sum(z1[, 3]) > 0){
      z22[index_z22, ] <- rmnom(sum(z1[, 3]), 1, theta22[i_index, ])
    }
    z22_vec <- as.numeric(z22 %*% 1:k22)
    
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
    U1[i, ] <- u1
    U2[i, ] <- u2
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
alpha21 <- 0.1; alpha22 <- 0.1
alpha31 <- 0.1; alpha32 <- 0.1; alpha33 <- 0.1
beta <- 0.5

##�p�����[�^�̐^�l
theta11 <- thetat11
theta12 <- thetat12
theta21 <- thetat21
theta22 <- thetat22
eta <- etat
phi <- phit
gamma <- gammat
omega <- omegat
lambda <- lambdat

##�p�����[�^�̏����l��ݒ�
#�g�s�b�N���z�̏����l
theta11 <- extraDistr::rdirichlet(hh, rep(1.0, k11))
theta12 <- extraDistr::rdirichlet(item, rep(1.0, k12))
theta21 <- extraDistr::rdirichlet(hh, rep(1.0, k21))
theta22 <- extraDistr::rdirichlet(item, rep(1.0, k22))

#�P�ꕪ�z�̏����l
eta <- extraDistr::rdirichlet(k11*k12, rep(1.0, r))   #�]���X�R�A���z�̏����l
phi <- extraDistr::rdirichlet(k21, rep(1.0, v))   #���[�U�[�̒P�ꕪ�z�̏����l
gamma <- extraDistr::rdirichlet(k22, rep(1.0, v))   #�A�C�e���̒P�ꕪ�z�̏����l
omega <- extraDistr::rdirichlet(r, rep(5.0, v))   #�]���X�R�A�̒P�ꕪ�z�̏����l

#�X�C�b�`���O�ϐ��̏����l
lambda <- matrix(1/s, nrow=hh, ncol=s)


##�p�����[�^�̊i�[�p�z��
THETA11 <- array(0, dim=c(hh, k11, R/keep))
THETA12 <- array(0, dim=c(item, k12, R/keep))
THETA21 <- array(0, dim=c(hh, k21, R/keep))
THETA22 <- array(0, dim=c(item, k22, R/keep))
ETA <- array(0, dim=c(k11*k12, r, R/keep))
PHI <- array(0, dim=c(k21, v, R/keep))
GAMMA <- array(0, dim=c(k22, v, R/keep))
OMEGA <- array(0, dim=c(r, v, R/keep))
LAMBDA <- array(0, dim=c(hh, s, R/keep))
U_SEG1 <- matrix(0, nrow=d, ncol=k11)
U_SEG2 <- matrix(0, nrow=d, ncol=k12)
SEG1 <- matrix(0, nrow=f, ncol=a)
SEG21 <- matrix(0, nrow=f, ncol=k21)
SEG22 <- matrix(0, nrow=f, ncol=k22)
storage.mode(U_SEG1) <- "integer"
storage.mode(U_SEG2) <- "integer"
storage.mode(SEG21) <- "integer"
storage.mode(SEG22) <- "integer"

##�f�[�^�ƃC���f�b�N�X�̐ݒ�
#�C���f�b�N�X�̐ݒ�
user_list <- user_vec <- list()
item_list <- item_vec <- list()
wd_list <- wd_vec <- list()
y_list <- y_ones <- list()
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
for(j in 1:r){
  y_list[[j]] <- which(y==j)
  y_ones[[j]] <- rep(1, length(y_list[[j]]))
}
index_k11 <- rep(1:k12, k11)
index_k12 <- rep(1:k11, rep(k12, k11))

#�f�[�^�̐ݒ�
y_vec <- y[d_id]
y_data <- matrix(as.numeric(table(1:f, y_vec)), nrow=f, ncol=r)
storage.mode(y_data) <- "integer"
r_vec <- rep(1, r)
a_vec <- rep(1, a)
vec11 <- rep(1, k11)
vec12 <- rep(1, k12)
vec21 <- rep(1, k21)
vec22 <- rep(1, k22)
K11 <- matrix(1:(k11*k12), nrow=k11, ncol=k12, byrow=T)
K12 <- matrix(1:(k11*k12), nrow=k12, ncol=k11, byrow=T)


##�ΐ��ޓx�̊�l
par <- colSums(WX) / f
LLst <- sum(WX %*% log(par))


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�]�_�X�R�A�̃��[�U�[�g�s�b�N���T���v�����O
  #���[�U�[�g�s�b�N�̏����t���m��
  eta11 <- t(eta)[y, ] * theta12[item_id, index_k11]
  par_u0 <- matrix(0, nrow=d, ncol=k11)
  for(j in 1:k11){
    par_u0[, j] <- eta11[, K1[j, ]] %*% vec11
  }
  par_u1 <- theta11[user_id, ] * par_u0   #���[�U�[�g�s�b�N�̊��Җޓx
  
  #���ݕϐ��̊����m������g�s�b�N���T���v�����O
  u1_rate <- par_u1 / as.numeric(par_u1 %*% vec11)
  Ui1 <- rmnom(d, 1, u1_rate)
  Ui1_T <- t(Ui1)
  
  
  ##�]�_�X�R�A�̃A�C�e���g�s�b�N���T���v�����O
  #�A�C�e���g�s�b�N�̏����t���m��
  eta12 <- t(eta)[y, ] * theta11[item_id, index_k12]
  par_u0 <- matrix(0, nrow=d, ncol=k12)
  for(j in 1:k12){
    par_u0[, j] <- eta12[, K1[, j]] %*% vec12
  }
  par_u2 <- theta12[item_id, ] * par_u0   #�A�C�e���g�s�b�N�̊��Җޓx
  
  #���ݕϐ��̊����m������g�s�b�N���T���v�����O
  u2_rate <- par_u2 / as.numeric(par_u2 %*% vec12)
  Ui2 <- rmnom(d, 1, u2_rate)
  Ui2_T <- t(Ui2)

  #���[�U�[�ƃA�C�e���g�s�b�N�𓝍�
  Ui <- Ui2[, index_k11] * Ui1[, index_k12]
  Ui_T <- t(Ui)
  
  
  ##�]�_�X�R�A�̃��[�U�[����уA�C�e���̃g�s�b�N���T���v�����O
  #���[�U�[�̃g�s�b�N���z���T���v�����O
  wusum0 <- matrix(0, nrow=d, ncol=k11)
  for(i in 1:hh){
    wusum0[i, ] <- Ui1_T[, user_index[[i]]] %*% user_ones[[i]]
  }
  wusum <- wusum0 + alpha11   #�f�B���N�����z�̃p�����[�^
  theta11 <- extraDistr::rdirichlet(hh, wusum)   #�f�B���N�����z����theta11���T���v�����O
  
  #�A�C�e���̃g�s�b�N���z���T���v�����O
  wisum0 <- matrix(0, nrow=item, ncol=k12)
  for(i in 1:item){
    wisum0[i, ] <- Ui2_T[, item_index[[i]]] %*% item_ones[[i]]
  }
  wisum <- wisum0 + alpha12   #�f�B���N�����z�̃p�����[�^
  theta12 <- extraDistr::rdirichlet(item, wisum)   #�f�B���N�����z����theta11���T���v�����O
  
  
  ##�]���X�R�A���z���T���v�����O
  vsum0 <- matrix(0, nrow=k11*k12, ncol=r)
  for(j in 1:r){
    vsum0[, j] <- Ui_T[, y_list[[j]]] %*% y_ones[[j]]
  }
  vsum <- vsum0 + alpha31   #�f�B���N�����z�̃p�����[�^
  eta <- extraDistr::rdirichlet(k11*k12, vsum)   #�f�B���N�����z����eta���T���v�����O
  
  
  ##�������z���X�C�b�`���O�ϐ����T���v�����O
  #�]���X�R�A�A���[�U�[����уA�C�e���̊��Җޓx��ݒ�
  Li_score <- as.numeric((t(omega)[wd, ] * y_data) %*% r_vec)   #�X�R�A�ޓx
  Li_user <- theta21[u_id, ] * t(phi)[wd, ]   #���[�U�[�ޓx
  par_user <- as.numeric(Li_user %*% vec21)   #���[�U�[�̊��Җޓx
  Li_item <- theta22[i_id, ] * t(gamma)[wd, ]   #�A�C�e���ޓx
  par_item <- as.numeric(Li_item %*% vec22)   #�A�C�e���̊��Җޓx
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
  Zi21 <- matrix(0, nrow=f, ncol=k21)
  Zi22 <- matrix(0, nrow=f, ncol=k22)
  Zi21[index_z21, ] <- rmnom(nrow(z_rate1), 1, z_rate1)
  Zi22[index_z22, ] <- rmnom(nrow(z_rate2), 1, z_rate2)
  Zi21_T <- t(Zi21)
  Zi22_T <- t(Zi22)
  
  
  ##�g�s�b�N���f���̃p�����[�^���T���v�����O
  #���[�U�[�̃g�s�b�N���z���T���v�����O
  wusum0 <- matrix(0, nrow=hh, ncol=k21)
  for(i in 1:hh){
    wusum0[i, ] <- Zi21_T[, user_list[[i]]] %*% user_vec[[i]]
  }
  wusum <- wusum0 + alpha21   #�f�B���N�����z�̃p�����[�^
  theta21 <- extraDistr::rdirichlet(hh, wusum)   #�f�B���N�����z����theta21���T���v�����O
  
  #�A�C�e���̃g�s�b�N���z���T���v�����O
  wisum0 <- matrix(0, nrow=item, ncol=k22)
  for(i in 1:item){
    wisum0[i, ] <- Zi22_T[, item_list[[i]]] %*% item_vec[[i]]
  }
  wisum <- wisum0 + alpha22   #�f�B���N�����z�̃p�����[�^
  theta22 <- extraDistr::rdirichlet(item, wisum)   #�f�B���N�����z����theta22���T���v�����O
  
  
  ##�]���X�R�A�A���[�U�[����уA�C�e���̒P�ꕪ�z���T���v�����O
  y_data_t <- t(y_data * Zi1[, 1])
  vssum0 <- matrix(0, nrow=r, ncol=v)
  vusum0 <- matrix(0, nrow=k21, ncol=v)
  visum0 <- matrix(0, nrow=k22, ncol=v)
  for(j in 1:v){
    vssum0[, j] <- y_data_t[, wd_list[[j]], drop=FALSE] %*% wd_vec[[j]]
    vusum0[, j] <- Zi21_T[, wd_list[[j]], drop=FALSE] %*% wd_vec[[j]]
    visum0[, j] <- Zi22_T[, wd_list[[j]], drop=FALSE] %*% wd_vec[[j]]
  }
  vssum <- vssum0 + alpha31; vusum <- vusum0 + alpha32; visum <- visum0 + alpha33
  omega <- extraDistr::rdirichlet(r, vssum)
  phi <- extraDistr::rdirichlet(k21, vusum)
  gamma <- extraDistr::rdirichlet(k22, visum)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA11[, , mkeep] <- theta11
    THETA12[, , mkeep] <- theta12
    THETA21[, , mkeep] <- theta21
    THETA22[, , mkeep] <- theta22
    ETA[, , mkeep] <- eta
    PHI[, , mkeep] <- phi
    GAMMA[, , mkeep] <- gamma
    OMEGA[, , mkeep] <- omega
    LAMBDA[, , mkeep] <- lambda
  }  
  
  #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
  if(rp%%keep==0 & rp >= burnin){
    U_SEG1 <- U_SEG1 + Ui1
    U_SEG2 <- U_SEG2 + Ui2
    SEG1 <- SEG1 + Zi1
    SEG21 <- SEG21 + Zi21
    SEG22 <- SEG22 + Zi22
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    LL1 <- sum(log(rowSums(par_u1)))
    LL2 <- sum(log(rowSums(Zi1 * par)))
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(c(LL1, LL2, LLst))
    print(round(c(colMeans(Zi1), colMeans(Z1)), 3))
    print(round(cbind(phi[, (v1-5):(v1+4)], phit[, (v1-5):(v1+4)]), 3))
  }
}

