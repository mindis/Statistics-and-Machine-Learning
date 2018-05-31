#####Mixture Hidden Marcov Model#####
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
k1 <- 5   #������
k2 <- 7   #HMM�̏�Ԑ�
hh <- 3000   #���[�U�[��
item <- 1000   #���A�C�e����
pt <- rpois(hh, rgamma(hh, 13, 0.5))   #���[�U�[���Ƃ̍w���@�
pt[pt < 5] <- ceiling(runif(sum(pt < 5), 5, 10))
hhpt <- sum(pt)   #�����R�[�h��
w <- rpois(hhpt, rgamma(hhpt, 10, 0.75))   #���Ԃ��Ƃ̍w����
w[w < 5] <- ceiling(runif(sum(w < 5), 5, 10))
f <- sum(w)   #���w����

#ID�̐ݒ�
u_id <- rep(1:hh, pt)
t_id <- c()
for(i in 1:hh){t_id <- c(t_id, 1:pt[i])}

#�C���f�b�N�X���쐬
id_list <- list()
for(i in 1:hh){id_list[[i]] <- which(u_id==i)}


##�p�����[�^�̐ݒ�
#�f�B���N�����z�̃p�����[�^��ݒ�
alpha01 <- rep(10.0, k1)
alpha02 <- rep(1.5, k2)
alpha03 <- matrix(0.2, nrow=k2, ncol=k2)
diag(alpha03) <- 2.0
beta <- rep(0.1, item)

#�f�B���N�����z����p�����[�^�𐶐�
theta1 <- thetat1 <- as.numeric(extraDistr::rdirichlet(1, alpha01))
theta2 <- thetat2 <- extraDistr::rdirichlet(k1, alpha02)
theta3 <- thetat3 <- array(0, dim=c(k2, k2, k1))
phi <- phit <- array(0, dim=c(k2, item, k1))
for(i in 1:1000){
  for(j in 1:k1){
    theta3[, , j] <- thetat3[, , j] <- extraDistr::rdirichlet(k2, alpha03)
    phi[, , j] <- phit[, , j] <- extraDistr::rdirichlet(k2, beta)
  }
  if(sum(phi==0)==0) break
}

##HHMM���f���Ɋ�Â��f�[�^�𐶐�
Z1 <- matrix(0, nrow=hh, ncol=k1)
Z2_list <- list()
Data <- matrix(0, nrow=hhpt, ncol=item)

for(i in 1:hh){
  #���ݕϐ����T���v�����O
  Z1[i, ] <- rmnom(1, 1, theta1)
  z1_vec <- which.max(Z1[i, ])
  
  #�������z���}���R�t��Ԑ��ڂ𐶐�
  z2 <- matrix(0, nrow=pt[i], ncol=k2)
  z2_vec <- rep(0, pt[i])
  data <- matrix(0, nrow=pt[i], ncol=item)
  
  for(j in 1:pt[i]){
    if(j==1){
   
      #�����̏�Ԃ𐶐�
      z2[j, ] <- rmnom(1, 1, theta2[z1_vec, ])
      z2_vec[j] <- which.max(z2[j, ])
    
    } else {
      
      #��Ԑ��ڂ𐶐�
      z2[j, ] <- rmnom(1, 1, theta3[z2_vec[j-1], , z1_vec[i]])
      z2_vec[j] <- as.numeric(z2 %*% 1:k2)
    }
    #�}���R�t��Ԑ��ڂɂ��ƂÂ��A�C�e���w���𐶐�
    data[j, ] <- as.numeric(rmnom(1, w[u_id==i & t_id==j], phi[z2_vec[j], , z1_vec]))
  }
  #���������f�[�^���i�[
  Z2_list[[i]] <- z2
  Data[id_list[[i]], ] <- data
}

#���X�g��ϊ�
Z2 <- do.call(rbind, Z2_list)
storage.mode(Data) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")
sparse_data_T <- t(sparse_data)


####�}���R�t�A�������e�J�����@��Mixture HMM���f���𐄒�####
##�A���S���Y���̐ݒ�
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##�p�����[�^�̐^�l
beta1 <- betat1
theta1 <- thetat1
theta2 <- thetat2
theta3 <- thetat3
phi <- phit


##HHMM���f���̏����l��ݒ�
#�������z�̖��x�֐��̑ΐ��ޓx�̒萔
const <- lfactorial(w) - rowSums(lfactorial(Data))   

#�p�����[�^�̏����l��ݒ�
beta1 <- rep(0.3, hh)
theta1 <- extraDistr::rdirichlet(1, rep(2.0, k1))
alpha <- matrix(0.5, nrow=k2, ncol=k2)
diag(alpha) <- 1.5
theta2 <- as.numeric(extraDistr::rdirichlet(1, rep(10.0, k1)))
theta3 <- array(0, dim=c(k2, k2, k1))
phi <- array(0, dim=c(k2, item, k1))
for(j in 1:k1){
  theta3[, , j] <- thetat3[, , j] <- extraDistr::rdirichlet(k2, alpha)
  phi0 <- extraDistr::rdirichlet(k2, colSums(Data)/sum(Data) * item) + 0.0001
  phi[, , j] <- phi0 / rowSums(phi0)
}

##���O���z�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z
alpha01 <- 1.0
alpha02 <- 0.1
beta01 <- 0.1

##�T���v�����O���ʂ̕ۑ��p�z��
THETA1 <- matrix(0, nrow=R/keep, ncol=k1)
THETA2 <- matrix(0, nrow=R/keep, ncol=k2)
THETA3 <- array(0, dim=c(k2, k2, k1, R/keep))
PHI <- array(0, dim=c(k2, item, k1, R/keep))
SEG1 <- rep(0, hhpt)
SEG2 <- matrix(0, nrow=hhpt, ncol=k1)
SEG3 <- matrix
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"
storage.mode(SEG3) <- "integer"


##MCMC�p�C���f�b�N�X���쐬
max_time <- max(t_id)
index_t11 <- which(t_id==1)
index_t21 <- list()
index_t22 <- list()
for(j in 2:max_time){
  index_t21[[j]] <- which(t_id==j)-1
  index_t22[[j]] <- which(t_id==j)
}
index_k2 <- matrix(1:(k2*k1), nrow=k1, ncol=k2, byrow=T)


#�ΐ��ޓx�̊�l
LLst <- sum(dmnom(Data, w, colSums(Data)/sum(Data), log=TRUE))


####MCMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�p�^�[�����Ƃɑΐ��ޓx���v�Z
  #�z��̐ݒ�
  rf12 <- matrix(0, nrow=k1, ncol=k1)
  rf22 <- array(0, dim=c(k2, k2, k1))
  Zi1 <- rep(0, hhpt)
  Zi2 <- matrix(0, nrow=hhpt, ncol=k1)
  z2_vec <- rep(0, hhpt)
  LLi1 <- matrix(0, nrow=hhpt, ncol=k1*k2)
  LLm <- matrix(0, nrow=hhpt, ncol=k1*k2)
  
  for(j in 1:k1){
    log_theta3 <- matrix(log(theta3[j, ]), nrow=hhpt, ncol=k2, byrow=T)
    LLm[, index_k2[j, ]] <- as.matrix(const + sparse_data %*% t(log(phi[, , j])))
    LLi1[, index_k2[j, ]] <- log_theta3 + LLm[, index_k2[j, ]]
  }
  
  for(pd in 1:max_time){
    if(pd==1){
      ##1���ڂ̍w�����Ԃ̐��ݏ�Ԃ𐶐�
      ##��ʊK�w�̐��ݏ�Ԃ𐶐�
      #�ΐ��ޓx����ʊK�w�̊��Җޓx�ɕϊ�
      LLt1 <- LLi1[index_t11, ]
      par0 <- exp(LLt1 - rowMaxs(LLt1))
      par1 <- matrix(0, nrow=hh, ncol=k1)
      for(j in 1:k1){
        par1[, j] <- rowSums(par0[, index_k2[j, ]])
      }
      
      #�������z�����ʊK�w�̏�Ԃ𐶐�
      r <- matrix(theta1, nrow=hh, ncol=k1, byrow=T)
      par_rate1 <- r * par1 / rowSums(r * par1)   #���ݕϐ��̊����m��
      Zi2[index_t11, ] <- rmnom(hh, 1, par_rate1)   #��Ԃ𐶐�
      z2_vec[index_t11] <- as.numeric(Zi2[index_t11, ] %*% 1:k1)
      
      
      ##���ʊK�w�̐��ݏ�Ԃ𐶐�
      #��ʊK�w�̏�Ԃɉ����Ėޓx���v�Z
      LLt2 <- matrix(0, nrow=hh, ncol=k2)
      zi2 <- Zi2[index_t11, ]
      for(j in 1:k1){
        LLt2 <- LLt2 + LLt1[, index_k2[j, ]] * Zi2[index_t11, j] * zi2[, j]
      }
      
      #�������z���牺�ʊK�w�̏�Ԃ𐶐�
      Zi3 <- matrix(0, nrow=hhpt, ncol=k2)
      z3_vec <- rep(0, hhpt)
      par2 <- exp(LLt2 - rowMaxs(LLt2))
      par_rate2 <- par2 / rowSums(par2)   #���ݕϐ��̊����m��
      Zi3[index_t11, ] <- rmnom(hh, 1, par_rate2)   #��Ԃ𐶐�
      z3_vec[index_t11] <- as.numeric(Zi3[index_t11, ] %*% 1:k2)
      
    } else {
      
      ##2���ڈȍ~�̃}���R�t���ڂɊ�Â����ݏ�Ԃ𐶐�
      ##���ݏ�Ԃ��؂�ւ�邩�ǂ����𐶐�
      #�f�[�^�̐ݒ�
      index1 <- index_t21[[pd]]
      index2 <- index_t22[[pd]]
      n <- length(index2)
      zi2_j <- Zi2[index1, , drop=FALSE]
      zi3_j <- Zi3[index1, , drop=FALSE]
      z3_j <- z3_vec[index1]
      
      #���Җޓx���v�Z
      LLz1 <- matrix(0, nrow=n, ncol=k1*k2)
      LLi2 <- matrix(0, nrow=n, ncol=k1*k2) 
      LLt2 <- matrix(0, nrow=n, ncol=k2)
      theta4_par <- matrix(0, nrow=n, ncol=k2) 
      
      for(l in 1:k1){
        theta4_par <- theta4_par + theta4[z3_j, , l] * zi2_j[, l]
        LLz1[, index_k2[l, ]] <- LLi1[index2, index_k2[l, ]] * (1-zi2_j[, l])
        LLo <- (log(theta4[z3_vec[index1], , l]) + LLm[index2, index_k2[l, ]]) * zi2_j[, l] 
        LLt2 <- LLt2 + LLo
        LLi2[, index_k2[l, ]] <- LLz1[, index_k2[l, ]] + LLo
      }
      
      #��ʊK�w�̐��ݏ�Ԑ؊����ϐ��𐶐�
      max_z <- rowMaxs(LLi2)
      LLz0 <- exp(LLz1 - max_z)
      LLz <- matrix(0, nrow=n, ncol=k1+1)
      for(l in 1:k1){
        LLz[, l] <- rowSums(LLz0[, index_k2[l, ], drop=FALSE])
      }
      
      LLz[, ncol(LLz)] <- rowSums(exp(LLt2 - max_z))
      LLz[, -ncol(LLz)] <- LLz[, -ncol(LLz)] * (1-zi2_j)
      r <- beta1[u_id[index2]]   #������
      beta_par <- r * rowSums(LLz[, 1:k1, drop=FALSE])
      beta_rate <- beta_par / (beta_par + (1-r)*LLz[, ncol(LLz), drop=FALSE])   #���ݕϐ��̊����m��
      Zi1[index2] <- rbinom(n, 1, beta_rate)   #�񍀕��z����؊����ϐ��𐶐�
      index_z1 <- which(Zi1[index2]==1)
      
      
      ##�؊����ϐ��������Ƃ��ď�ʊK�w�𐶐�
      if(length(index_z1)==0){
        Zi2[index2, ] <- Zi2[index1, ]
      } else {
        r <- theta2[(Zi2 %*% 1:k1)[index1], ]
        Zi2[index2, ] <- Zi2[index1, ]   #1���O�̐��ݏ�Ԃ��J��z��
        par_rate1 <- r * LLz[, -ncol(LLz), drop=FALSE] / rowSums(r * LLz[, -ncol(LLz), drop=FALSE])   #���ݏ�Ԃ̊����m��
        
        if(nrow(par_rate1)==1){
          Zi2[index2, ] <- rmnom(length(index_z1), 1, par_rate1[index_z1, ])   #���ݏ�Ԃ𐶐�
          z2_vec[index2] <- as.numeric(Zi2[index2, ] %*% 1:k1)
        } else {
          Zi2[index2, ][index_z1, ] <- rmnom(length(index_z1), 1, par_rate1[index_z1, ])   #���ݏ�Ԃ𐶐�
          z2_vec[index2][index_z1] <- as.numeric(Zi2[index2, ][index_z1, ] %*% 1:k1)
        }
        #��ʊK�w�̃}���R�t���ڍs����X�V
        rf12 <- rf12 + t(Zi2[index1, , drop=FALSE]) %*% (Zi2[index2, , drop=FALSE] * Zi1[index2])
      }
      
      
      ##���ʊK�w�̐��ݏ�Ԃ𐶐�
      #��ʊK�w�̏�Ԃɉ����Ėޓx���v�Z
      LLo <- matrix(0, nrow=n, ncol=k2)
      for(j in 1:k1){
        LLo <- LLo + (LLi1[index2, , drop=FALSE] * Zi1[index2])[, index_k2[j, ], drop=FALSE] * Zi2[index2, j]
      }
      LLt3 <- LLo + LLt2*(1-Zi1[index2])
      
      #�������z���牺�ʊK�w�̏�Ԃ𐶐�
      par2 <- exp(LLt3 - rowMaxs(LLt3))
      par_rate2 <- par2 / rowSums(par2)   #���ݕϐ��̊����m��
      Zi3[index2, ] <- rmnom(n, 1, par_rate2)   #��Ԃ𐶐�
      z3_vec[index2] <- as.numeric(Zi3[index2, ] %*% 1:k2)
      
      #�}���R�t���ڍs����X�V
      for(j in 1:k1){
        rf22[, , j] <- rf22[, , j] + t(Zi3[index1, , drop=FALSE]) %*% 
          (Zi3[index2, , drop=FALSE] * Zi2[index2, j] * (1-Zi1[index2]))
      }
    }
  }
  
  ##�p�����[�^�̍X�V
  #�x�[�^���z���獬�������T���v�����O
  par1 <- as.numeric(tapply(Zi1[-index_t11], u_id[-index_t11], sum))
  par2 <- pt-1-par1
  beta1 <- rbeta(hh, par1+alpha01, par2+beta01)
  
  #�f�B�N�������z�����ʊK�w�̍��������T���v�����O
  wf1 <- colSums(Zi2[index_t11, ]) + beta02
  wf2 <- rf12 + beta02
  theta1 <- extraDistr::rdirichlet(1, wf1)
  theta2 <- extraDistr::rdirichlet(k1, wf2)
  diag(theta2) <- 0
  theta2 <- theta2 / rowSums(theta2)
  
  #�f�B�N�������z���牺�ʊK�w�̍��������T���v�����O
  for(j in 1:k1){
    wf3 <- colSums(Zi3[index_t11, ] * Zi2[index_t11, j]) + colSums(Zi3 * Zi1 * Zi2[, j]) + beta02
    wf4 <- rf22[, , j] + beta02
    theta3[j, ] <- extraDistr::rdirichlet(1, wf3)
    theta4[, , j] <- extraDistr::rdirichlet(k2, wf4)
  }
  
  #�f�B�N�������z����A�C�e���w���m���̃p�����[�^���T���v�����O
  for(j in 1:k1){
    vf <- t(sparse_data_T %*% (Zi3 * Zi2[, j])) + beta03
    phi[, , j] <- extraDistr::rdirichlet(k2, vf)
  }
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA1[mkeep, ] <- theta1
    THETA2[, , mkeep] <- theta2
    THETA3[, , mkeep] <- theta3
    THETA4[, , , mkeep] <- theta4
    PHI[, , , mkeep] <- phi
    
    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(mkeep >= burnin & rp%%keep==0){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
      SEG3 <- SEG3 + Zi3
    }
    
    #�T���v�����O���ʂ��m�F
    if(rp%%disp==0){
      print(rp)
      #�ΐ��ޓx���v�Z
      LLi <- matrix(0, nrow=hhpt, ncol=k1)
      for(j in 1:k1){
        LLi[, j] <- rowSums(sparse_data %*% log(t(phi[, , j])) * Zi3 * Zi2[, j])
      }
      LL <- sum(const + rowSums(LLi))
      
      #�p�����[�^��\��
      print(c(LL, LLst))
      print(round(rbind(theta1, thetat1), 3))
      print(round(cbind(theta2, thetat2), 3))
      print(round(rbind(theta3, thetat3), 3))
    }
  }
}


