#####Collaborative Topic Model#####
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
library(ggplot2)

#set.seed(93441)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
k <- 10
hh <- 1500   #���r���[�l��
item <- 500   #�A�C�e����
s <- rpois(hh, rgamma(hh, 7.5, 0.6))   #1�l������̃��r���[��
s[s==0] <- ceiling(runif(length(s[s==0]), 1, 5))
d <- sum(s)   #��������
v <- 1000   #��b��
w <- rpois(d, rgamma(d, 40, 0.8))   #����������̒P�ꐔ
f <- sum(w)   #���P�ꐔ

##�A�C�e���w����ݒ�
u_list <- list()
par <- as.numeric(extraDistr::rdirichlet(1, rep(3.0, item)))   #�A�C�e���w���m��
for(rp in 1:1000){
  for(i in 1:hh){
    for(j in 1:1000){
      pi <- rmnom(s[i], 1, par)
      if(max(colSums(pi))==1){
        break
      }
    }
    u_list[[i]] <- pi
  }
  if(min(colSums(do.call(rbind, u_list)))==0){
    break
  }
}
U <- do.call(rbind, u_list)
u_vec <- as.numeric(U %*% 1:item)
colSums(U)

##ID��ݒ�
u_id <- rep(1:hh, s)
w_id <- rep(u_vec, w)


##�p�����[�^�̐ݒ�
#�g�s�b�N���f���̃p�����[�^
alpha01 <- rep(0.2, k)
alpha02 <- rep(0.15, v)
theta1 <- thetat1 <- extraDistr::rdirichlet(item, alpha01)
phi <- phit <- extraDistr::rdirichlet(k, alpha02)

#���݈��q���f���̃p�����[�^
tau1 <- 0.5
tau2 <- 0.025
sigma <- 0.2
theta2 <- thetat2 <- mvrnorm(hh, rep(0, k), diag(tau1, k))
psi <- psit <- t(theta1 + mvrnorm(item, rep(0, k), diag(tau2, k)))


##Collavorative topic model�̃f�[�^�𐶐�
Z_list <- list()
WX <- matrix(0, nrow=d, ncol=v)
word_list <- list()
score <- rep(0, d)

for(i in 1:d){
  #�A�C�e���g�s�b�N�ƒP��𐶐�
  z <- rmnom(w[i], 1, theta1[u_vec[i], ])   #�g�s�b�N�𐶐�
  z_vec <- as.numeric(z %*% 1:k)
  words <- rmnom(w[i], 1, phi[z_vec, ])   #�P��𐶐�
  words_vec <- colSums(words)
  
  #�X�R�A�𐶐�
  r <- as.numeric(theta2[u_id[i], ] %*% psi[, u_vec[i]]) + rnorm(1, 0, sigma)
  
  #�p�����[�^���i�[
  Z_list[[i]] <- z
  WX[i, ] <- words_vec
  word_list[[i]] <- as.numeric(words %*% 1:v)
  score[i] <- r
}

#���X�g��ϊ�
Z <- do.call(rbind, Z_list)
words_vec <- unlist(word_list)
storage.mode(WX) <- "integer"


####EM�A���S���Y����Collaborative Topic Model�𐄒�####
##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, w, k){
  #���S�W�����v�Z
  Bur <- theta1[w_id, ] * t(phi)[words_vec, ]   #�ޓx
  Br <- Bur / rowSums(Bur)   #���S��
  r <- colSums(Br) / sum(Br)   #������
  bval <- list(Br=Br, Bur=Bur, r=r)
  return(bval)
}

##�C���f�b�N�X�̍쐬
#�g�s�b�N���f���̃C���f�b�N�X���쐬
doc_list <- list()
word_list <- list()
doc_ones <- list()
word_ones <- list()
doc_n <- rep(0, item)
word_n <- rep(0, v)

for(i in 1:item){
  doc_list[[i]] <- which(w_id==i)
  doc_n[i] <- length(doc_list[[i]])
  doc_ones[[i]] <- rep(1, doc_n[i])
}
for(i in 1:v){
  word_list[[i]] <- which(words_vec==i)
  word_n[i] <- length(word_list[[i]])
  word_ones[[i]] <- rep(1, word_n[i])
}

#���݈��q���f���̃C���f�b�N�X���쐬
user_list <- list()
user_n <- rep(0, hh)
item_list <- list()
item_n <- rep(0, item)
for(i in 1:hh){
  user_list[[i]] <- which(u_id==i)
  user_n[i] <- length(user_list[[i]])
}
for(i in 1:item){
  item_list[[i]] <- which(u_vec==i)
  item_n[i] <- length(item_list[[i]])
}


#�����l�̐ݒ�
theta1 <- extraDistr::rdirichlet(item, rep(0.5, k))
phi <- extraDistr::rdirichlet(k, rep(1.0, v))
theta2 <- extraDistr::rdirichlet(hh, rep(0.5, k))
sigma <- 0.2   #�덷�W���΍��̏����l
tau1 <- 0.025   #�A�C�e�����q�̃n�C�p�[�p�����[�^
tau2 <- 0.25   #���[�U�[���q�̃n�C�p�[�p�����[�^


#EM�A���S���Y���̍X�V�X�e�[�^�X
iter <- 1
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 0.1
LLo <- -1000000000   #�ΐ��ޓx�̏����l
LLw <- c()

##EM�A���S���Y���Ńp�����[�^���X�V
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  
  ##�g�s�b�N���f���ŒP��̃g�s�b�N�𐄒�
  ##�ޓx�Ɛ��ݕϐ�z�𐄒�
  par <- burden_fr(theta1, phi, words_vec, w_id, k)
  z <- par$Br   #���ݕϐ�z���o��
  r <- par$r   #������r���o��
  
  ##�g�s�b�N���f���̃p�����[�^���X�V
  #�g�s�b�N���ztheta1���X�V
  wsum <- matrix(0, nrow=item, ncol=k)
  for(i in 1:item){
    wsum[i, ] <- t(z[doc_list[[i]], ]) %*% doc_ones[[i]]
  }
  theta1 <- wsum / doc_n
  theta1[is.nan(theta1)] <- rep(1/k, k)
  
  #�P�ꕪ�zphi���X�V
  vsum <- matrix(0, nrow=k, ncol=v)
  for(j in 1:v){
    vsum[, j] <- t(z[word_list[[j]], ]) %*% word_ones[[j]]
  }
  phi <- vsum / matrix(rowSums(vsum), nrow=k, ncol=v)
  
  
  ##���݈��q���f���ŃX�R�A�v���𐄒�
  ##�A�C�e�����qpsi�𐄒�
  #�f�[�^�̐ݒ�
  x0 <- thetat2[u_id, ]
  lambda <- diag(tau1, k)
  
  #�A�C�e�����ƂɈ��q�s����œK��
  for(j in 1:item){
    
    #�A�C�e�����Ȃ��ꍇ�͎��̃A�C�e����
    if(item_n[j]==0){
      psi[, j] <- theta1[j, ]
      next
    }
    #�A�C�e�����ƂɃf�[�^��ݒ�
    index <- item_list[[j]]
    x <- x0[index, ]
    y <- score[index]
    omega <- diag(sigma^2, item_n[j])
    
    #psi�𐳑����ŏ����@�ōœK��
    if(item_n[j] > 1){
      psi[, j] <- as.numeric(solve(t(x) %*% omega %*% x + lambda) %*% (t(x) %*% omega %*% y + lambda %*% theta1[j, ]))
    } else {
      psi[, j] <- as.numeric(solve(x %*% omega %*% x + lambda) %*% (x %*% omega %*% y + lambda %*% theta1[j, ]))
    }
  }
  psi[is.nan(psi)] <- 0
  
  ##���[�U�[���qtheta2�𐄒�
  #�f�[�^�̐ݒ�
  x0 <- t(psi)[u_vec, ]
  lambda <- diag(tau2, k)
  
  #���[�U�[���ƂɈ��q�s����œK��
  for(i in 1:hh){
    #���[�U�[���ƂɃf�[�^��ݒ�
    index <- user_list[[i]]
    x <- x0[index, ]
    y <- score[index]
    omega <- diag(sigma, user_n[i])

    #theta2�𐳑����ŏ����@�ōœK��
    if(user_n[i] > 1){
      theta2[i, ] <- as.numeric(solve(t(x) %*% omega %*% x + lambda) %*% t(x) %*% omega %*% y)
    } else {
      theta2[i, ] <- as.numeric(solve(x %*% omega %*% x + lambda) %*% x %*% omega %*% y)
    }
  }
  
  ##�덷sigma���X�V
  mu <- rowSums(theta2[u_id, ] * t(psi)[u_vec, ])
  er <- score - mu
  sigma <- sum(er^2) / (d-1)

  
  ##EM�A���S���Y���̍X�V
  #�ΐ��ޓx�̌v�Z
  LLS1 <- sum(log(rowSums(par$Bur)))   #�g�s�b�N���f���̑ΐ��ޓx
  LLS2 <- -sum(dnorm(score, mu, sqrt(sigma), log=TRUE))   #���݈��q���f���̑ΐ��ޓx
  LLS <- LLS1 + LLS2
  
  #��������
  iter <- iter+1
  dl <- LLS-LLo
  LLo <- LLS
  LLw <- c(LLw, LLo)
  print(LLo)
}

#####���茋�ʂ̗v��#####
##�p�����[�^����l
round(cbind(par$Br, Z %*% 1:k), 3)   #�����g�s�b�N
round(cbind(theta1, thetat1), 2)   #�g�s�b�N���z
round(cbind(t(phi), t(phit)), 3)   #�P�ꕪ�z
round(data.frame(n=user_n, x=theta2, y=thetat2), 2)   #���[�U�[���q
round(data.frame(n=item_n, x=t(psi), y=t(psit)), 2)   #�A�C�e�����q

##�X�R�A�̗\���l
pred_mu0 <- theta2 %*% psi
pred_mu1 <- thetat2 %*% psit
round(cbind(pred_mu0[, 1:10], pred_mu1[, 1:10]), 2)


