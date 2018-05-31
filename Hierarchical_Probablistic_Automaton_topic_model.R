#####Hierarchical Probablistic Automaton topic model#####
options(warn=0)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(data.table)
library(ggplot2)

#set.seed(5723)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
k1 <- 6   #��ʊK�w
k21 <- 10   #�g�s�b�N��
k22 <- 4   #�@�\���BIES�I�[�g�}�g��
k3 <- 4   #�g�s�b�N���BIES�I�[�g�}�g��
d <- 3000   #������
v1 <- 700   #��ʌ�
v2 <- 400   #�Ɨ���
v3 <- 100   #�@�\��
v <- v1 + v2 + v3   #����b��
w <- rpois(d, rgamma(d, 70, 0.4))   #��b��
f <- sum(w)   #����b��
B <- 1; I <- 2; E <- 3; S <- 4   #BIES�̃C���f�b�N�X
topic_index <- matrix(1:v1, nrow=k21, ncol=v1/k21, byrow=T)
functional_index <- matrix((v1+1):(v1+v2), nrow=k1-1, ncol=v2/(k1-1), byrow=T)

#ID�̐ݒ�
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}


##�p�����[�^�̎��O���z�̐ݒ�
#�������z�̎��O���z
pi01 <- c(rep(30, 2), rep(5, k1-2))

#�}���R�t���ڍs��̎��O���z
alpha011 <- c(60.0, rep(15.0, k1-1))   #�g�s�b�N�����Ă̏����t���ŏ�ʊK�w�̎��O���z
alpha012 <- c(10.0, rep(1.0, k1-1))   #�@�\�ꊄ���̍ŏ�ʊK�w�̎��O���z
alpha021 <- rep(0.1, k21)   #�g�s�b�N���z�̎��O���z
alpha022 <- c(15, 30)   #BIE��stop word�̎��O���z
alpha023 <- c(10, 20)   #�@�\���BIES�I�[�g�}�g���̎��O���z
alpha031 <- c(15, 30)   #BIE��stop word�̎��O���z
alpha032 <- c(10, 20)   #�g�s�b�N���BIES�I�[�g�}�g���̎��O���z


##�S�P�ꂪ���������܂Ńf�[�^�̐������p��
for(rp in 1:1000){
  print(rp)
  
  ##�p�����[�^�̐���
  #�}���R�t���ڍs��̃p�����[�^�𐶐�
  eta <- etat <- as.numeric(extraDistr::rdirichlet(1, pi01))
  theta1 <- thetat1 <- rbind(as.numeric(extraDistr::rdirichlet(1, alpha011)), extraDistr::rdirichlet(k1-1, alpha012))
  theta21 <- thetat21 <- extraDistr::rdirichlet(d, alpha021)
  beta21 <- betat21 <- rbeta(k1-1, alpha022[1], alpha022[2])
  beta22 <- betat22 <- rbeta(k1-1, alpha023[1], alpha023[2])
  beta31 <- betat31 <- rbeta(k21, alpha031[1], alpha031[2])
  beta32 <- betat32 <- rbeta(k21, alpha032[1], alpha032[2])
  
  ##�g�s�b�N�̒P�ꕪ�z�𐶐�
  phi <- array(0, dim=c(k21, v, k22))
  
  for(j in 1:k21){
    #���O���z��ݒ�
    alpha11 <- c(rep(0.15, v1), rep(0.001, v2), rep(0.001, v3))   #�g�s�b�N���B�̎��O���z
    alpha12 <- c(rep(0.25, v1), rep(0.001, v2), rep(0.2, v3))   #�g�s�b�N���I�̎��O���z
    alpha13 <- c(rep(0.15, v1), rep(0.001, v2), rep(0.01, v3))   #�g�s�b�N���E�̎��O���z
    alpha14 <- c(rep(0.15, v1), rep(0.001, v2), rep(0.001, v3))   #�g�s�b�N���S�̎��O���z
    
    #�g�s�b�N�ɉ����Ď��O���z��ݒ�
    alpha11[as.numeric(t(topic_index[-j, ]))] <- 0.015
    alpha12[as.numeric(t(topic_index[-j, ]))] <- 0.025
    alpha13[as.numeric(t(topic_index[-j, ]))] <- 0.015
    alpha14[as.numeric(t(topic_index[-j, ]))] <- 0.015
    
    #�p�����[�^�𐶐�
    phi[j, , 1] <- extraDistr::rdirichlet(1, alpha11)
    phi[j, , 2] <- extraDistr::rdirichlet(1, alpha12)
    phi[j, , 3] <- extraDistr::rdirichlet(1, alpha13)
    phi[j, , 4] <- extraDistr::rdirichlet(1, alpha14)
  }
  phit <- phi
  
  ##�@�\��̒P�ꕪ�z�𐶐�
  omega <- array(0, dim=c(k1-1, v, k3))
  
  for(j in 1:(k1-1)){
    #���O���z��ݒ�
    alpha21 <- c(rep(0.001, v1), rep(1.0, v2), rep(0.001, v3))   #�Ɨ����B�̎��O���z
    alpha22 <- c(rep(0.001, v1), rep(0.75, v2), rep(0.5, v3))   #�Ɨ����I�̎��O���z
    alpha23 <- c(rep(0.001, v1), rep(1.0, v2), rep(0.1, v3))   #�Ɨ����E�̎��O���z
    alpha24 <- c(rep(0.001, v1), rep(1.0, v2), rep(0.001, v3))   #�Ɨ����S�̎��O���z
    
    #�g�s�b�N�ɉ����Ď��O���z��ݒ�
    alpha21[as.numeric(t(functional_index[-j, ]))] <- 0.1
    alpha22[as.numeric(t(functional_index[-j, ]))] <- 0.15
    alpha23[as.numeric(t(functional_index[-j, ]))] <- 0.1
    alpha24[as.numeric(t(functional_index[-j, ]))] <- 0.1
   
    #�p�����[�^�𐶐�
    omega[j, , 1] <- extraDistr::rdirichlet(1, alpha21)
    omega[j, , 2] <- extraDistr::rdirichlet(1, alpha22)
    omega[j, , 3] <- extraDistr::rdirichlet(1, alpha23)
    omega[j, , 4] <- extraDistr::rdirichlet(1, alpha24)
  }
  omegat <- omega
  
  #��؂蕶���̃p�����[�^�𐶐�
  delta <- rbeta(d, 15, 60)
  
  
  ##���f���Ɋ�Â��f�[�^�𐶐�
  z1_list <- z2_list <- z3_list <- topic_list <- r2_list <- r3_list <- delimiter_list <- wd_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  
  for(i in 1:d){
    if(i%%100==0){
      print(i)
    }
    #�f�[�^�̊i�[�p�z��
    z1 <- matrix(0, nrow=w[i], ncol=k1)
    z1_vec <- rep(0, w[i])
    delimiter <- rep(0, w[i])
    topic <- matrix(0, nrow=w[i], ncol=k21)
    topic_vec <- rep(0, w[i])
    r3 <- r2 <- rep(0, w[i])
    z2 <- matrix(0, nrow=w[i], ncol=k22)
    z3 <- matrix(0, nrow=w[i], ncol=k3)
    word <- matrix(0, nrow=w[i], ncol=v)
    
    #1�P�ꂲ�Ƃɐ��ݕϐ��ƒP��𐶐�
    for(j in 1:w[i]){
      
      #�ŏ�ʊK�w�𐶐�
      if(j==1 | sum(delimiter[j-1])==1){
        z1[j, ] <- as.numeric(rmnom(1, 1, eta))
        z1_vec[j] <- which.max(z1[j, ])
      }
      if(j > 1 & sum(delimiter[j-1])==0){
        if(sum(z2[j-1, c(E, S)]) > 0 | sum(z3[j-1, c(E, S)]) > 0){
          z1[j, ] <- rmnom(1, 1, theta1[z1_vec[j-1], ])
          z1_vec[j] <- which.max(z1[j, ])
        }
        if(sum(z2[j-1, c(E, S)])==0 & sum(z3[j-1, c(E, S)])==0){
          z1[j, ] <- z1[j-1, ]
          z1_vec[j] <- which.max(z1[j-1, ])      
        }
      }
      
      #�@�\��̃X�C�b�`���O�ϐ��ƃg�s�b�N�𐶐�
      if(z1_vec[j] > 1){
        if(j==1 | sum(z2[j-1, c(E, S)]) > 0 | sum(z3[j-1, c(E, S)]) > 0){
          #�@�\��̃X�C�b�`���O�ϐ��𐶐�
          r2[j] <- rbinom(1, 1, beta21[z1_vec[j]-1])
          z2[j, S] <- 1-r2[j]
        }
      }
      
      if(z1_vec[j]==1){
        if(j==1 | sum(z2[j-1, c(E, S)]) > 0 | sum(z3[j-1, c(E, S)]) > 0){
          #�g�s�b�N�𐶐�
          topic[j, ] <- rmnom(1, 1, theta21[i, ])
          topic_vec[j] <- which.max(topic[j, ])
            
          #�g�s�b�N�̃X�C�b�`���O�ϐ��𐶐�
          r3[j] <- rbinom(1, 1, beta31[topic_vec[j]])
          z3[j, S] <- 1-r3[j]
        }
        if(j > 1 & sum(z2[j-1, c(E, S)])==0 & sum(z3[j-1, c(E, S)])==0){
          topic[j, ] <- topic[j-1, ]
          topic_vec[j] <- topic_vec[j-1]
        }
      }
      
      #�@�\��ƃg�s�b�N���BIE�𐶐�
      #�����lB�𐶐�
      if(z1_vec[j] > 0){
        if(z1_vec[j] > 1 & r2[j]==1){
          z2[j, B] <- 1
        }
        if(z1_vec[j]==1 & r3[j]==1){
          z3[j, B] <- 1
        }
      }
      
      #�p���lI�𐶐�
      if(j > 1 & sum(z2[j-1, B])==1){
        z2[j, I] <- 1
        z1[j, ] <- z1[j-1, ]
        z1_vec[j] <- z1_vec[j-1]
      }
      
      if(j > 1 & sum(z3[j-1, B])==1){
        z3[j, I] <- 1
        z1[j, ] <- z1[j-1, ]
        z1_vec[j] <- z1_vec[j-1]
        topic[j, ] <- topic[j-1, ]
        topic_vec[j] <- topic_vec[j-1]
      }
    
      #BIE���I��邩�ǂ����𐶐�
      if(sum(z2[j-1, I])==1){
        x <- rbinom(1, 1, beta22[z1_vec[j-1]-1])
        z2[j, I] <- x
        z2[j, E] <- 1-x
      }
      if(sum(z3[j-1, I])==1){
        x <- rbinom(1, 1, beta32[topic_vec[j-1]])
        z3[j, I] <- x
        z3[j, E] <- 1-x
      }
      
      #��؂蕶���𐶐�
      if(sum(z2[j, c(E, S)]) > 0 | sum(z3[j, c(E, S)]) > 0){
        delimiter[j] <- rbinom(1, 1, delta[i])
      }
      
      #�P��𐶐�
      if(z1_vec[j] > 1){
        word[j, ] <- rmnom(1, 1, omega[z1_vec[j]-1, , which.max(z2[j, ])])
      }
      if(z1_vec[j]==1){
        word[j, ] <- rmnom(1, 1, phi[topic_vec[j], , which.max(z3[j, ])])
      }
    }
    
    #�f�[�^���i�[
    z1_list[[i]] <- z1
    z2_list[[i]] <- z2
    z3_list[[i]] <- z3
    topic_list[[i]] <- topic
    r2_list[[i]] <- r2
    r3_list[[i]] <- r3
    delimiter_list[[i]] <- delimiter
    wd_list[[i]] <- as.numeric(word %*% 1:v)
    WX[i, ] <- colSums(word)
  }
  if(min(colSums(WX)) > 0) break
}

#���X�g��ϊ�
Z1 <- do.call(rbind, z1_list)
Z2 <- do.call(rbind, z2_list)
Z3 <- do.call(rbind, z3_list)
topic <- do.call(rbind, topic_list)
r2 <- unlist(r2_list)
r3 <- unlist(r3_list)
delimiter <- unlist(delimiter_list)
wd <- unlist(wd_list)

#�X�p�[�X�s��ɕϊ�
sparse_data <- sparseMatrix(1:f, wd, dims=c(f, v))
sparse_data_T <- t(sparse_data)


####�}���R�t�A�������e�J�����@��Hierarchical Probablistic Automaton topic model�𐄒�####
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
burnin <- 1000/keep
disp <- 10

##���O���z�̐ݒ�
alpha1 <- 0.1
alpha2 <- 1/k21
alpha3 <- c(1, 1) 
alpha4 <- 0.01

#�p�����[�^�̐^�l
eta <- etat
theta1 <- thetat1
theta21 <- thetat21
beta21 <- betat21
beta22 <- betat22
beta31 <- betat31
beta32 <- betat32
phi <- phit
omega <- omegat


##BIES�̒P�ꂲ�Ƃ̖ޓx�֐���ݒ�
#�g�s�b�N���BIES�ޓx
Lit_B <- theta21_d * t(phi[, , B])[wd, ]; Lit_I <- theta21_d * t(phi[, , I])[wd, ]
Lit_E <- theta21_d * t(phi[, , E])[wd, ]; Lit_S <- theta21_d * t(phi[, , S])[wd, ]

#�@�\����BIES�ޓx
Lif_B <- t(omega[, , B])[wd, ]; Lif_I <- t(omega[, , I])[wd, ]
Lif_E <- t(omega[, , E])[wd, ]; Lif_S <- t(omega[, , S])[wd, ]


eta

vec_topic <- rep(1, k21)
vec_functional <- rep(1, k1-1)

Lit_B %*% vec_topic
Lit_S %*% vec_topic

##��ʊK�w�̊��Җޓx�𐄒�
#�g�s�b�N���z�̊K�w�I���Җޓx
beta31_d <- matrix(beta31, nrow=f, ncol=k21, byrow=T)
Lit_switch1 <- beta31_d * Lit_B; Lit_switch2 <- (1-beta31_d) * Lit_S 
Lit_switch <- eta[1] * ((Lit_switch1 + Lit_switch2) %*% vec_topic)

Lit_I
theta


#�@�\�ꕪ�z�̊K�w�I���Җޓx
beta21_d <- matrix(beta21 , nrow=f, ncol=k1-1, byrow=T)
Lif_switch1 <- beta21_d * Lif_B; Lif_switch2 <- (1-beta21_d) * Lif_S
Lif_switch <- (matrix(eta[-1], nrow=f, ncol=k1-1, byrow=T) * (Lif_switch1 + Lif_switch2)) %*% vec_functional

Li_switch1 <- cbind(Lit_switch, Lif_switch)



cbind(round(Li_switch1[index_t11, ] / rowSums(Li_switch1[index_t11, ]), 3), Z1[index_t11, 1], rowSums(Z1[index_t11, -1]) > 0)










