#####Multinomial Topic Regression Model####
options(warn=2)
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

#set.seed(5723)

####�f�[�^�̔���####
dataset <- 2   #�f�[�^�Z�b�g��
k <- 15   #�g�s�b�N��
hh <- 2000   #���[�U�[��
page <- 400   #�y�[�W��
item <- 1000   #�A�C�e����
w1 <- rpois(hh, rgamma(hh, 35, 0.3))   #1�l�������̃y�[�W�{����
w2 <- rpois(hh, rgamma(hh, 32, 0.3))   #1�l������̃A�C�e���w����
f1 <- sum(w1)   #���y�[�W�{����
f2 <- sum(w2)   #���A�C�e���w����


#ID�̐ݒ�
u_id1 <- rep(1:hh, w1)
u_id2 <- rep(1:hh, w2)
t_id1 <- c()
t_id2 <- c()
for(i in 1:hh){
  t_id1 <- c(t_id1, 1:w1[i])
  t_id2 <- c(t_id2, 1:w2[i])
}

##�f�U�C���s��̐ݒ�
Data <- cbind(rep(1, hh*dataset), rep(0:(dataset-1), rep(hh, dataset)))
index_logit <- rep(1:hh, dataset)
index_data <- matrix(1:(hh*dataset), nrow=hh, ncol=dataset)

##�p�����[�^�̐ݒ�
for(rp in 1:100){
  print(rp)
  #�g�s�b�N���z�̃p�����[�^
  alpha01 <- alphat01 <- cbind(mvrnorm(hh, rep(0, k-1), diag(4.0, k-1)), 0)
  alpha02 <- alphat02 <- cbind(mvrnorm(hh, rep(0, k-1), diag(1.0, k-1)), 0)
  alpha11 <- rep(0.15, page)
  alpha12 <- rep(0.1, item)
  
  #�p�����[�^�𐶐�
  logit <- Data[, 1]*alpha01[index_logit, ] + Data[, 2]*alpha02[index_logit, ]   #���W�b�g�̌v�Z
  theta1 <- thetat1 <- exp(logit[index_data[, 1], ]) / rowSums(exp(logit[index_data[, 1], ]))
  theta2 <- thetat2 <- exp(logit[index_data[, 2], ]) / rowSums(exp(logit[index_data[, 2], ]))
  phi0 <- t(extraDistr::rdirichlet(page, rep(0.01, k))) * 
    (matrix(extraDistr::rdirichlet(1, rep(2.0, page)), nrow=k, ncol=page, byrow=T))
  phi <- phit <- phi0 / rowSums(phi0)
  gamma0 <- t(extraDistr::rdirichlet(item, rep(0.01, k))) * 
    (matrix(extraDistr::rdirichlet(1, rep(2.0, item)), nrow=k, ncol=item, byrow=T))
  gamma <- gammat <- gamma0 / rowSums(gamma0)
  
  ##���f���ɂ��ƂÂ��P��𐶐�����
  WX1 <- matrix(0, nrow=hh, ncol=page)
  WX2 <- matrix(0, nrow=hh, ncol=item)
  Z1_list <- list()
  Z2_list <- list()
  wd1_list <- list()
  wd2_list <- list()
  
  for(i in 1:hh){
    #�������z����g�s�b�N�𐶐�
    z1 <- rmnom(w1[i], 1, theta1[i, ])
    z1_vec <- as.numeric(z1 %*% 1:k)
    z2 <- rmnom(w2[i], 1, theta2[i, ])
    z2_vec <- as.numeric(z2 %*% 1:k)
    
    #�g�s�b�N���z���痚���𐶐�
    data1 <- rmnom(w1[i], 1, phi[z1_vec, ])
    data2 <- rmnom(w2[i], 1, gamma[z2_vec, ])
    
    #�f�[�^���i�[
    WX1[i, ] <- colSums(data1)
    WX2[i, ] <- colSums(data2)
    Z1_list[[i]] <- z1
    Z2_list[[i]] <- z2
    wd1_list[[i]] <- as.numeric(data1 %*% 1:page)
    wd2_list[[i]] <- as.numeric(data2 %*% 1:item)
  }
  if(min(colSums(WX1)) > 0 & min(colSums(WX2)) > 0) break
}

#���X�g��ϊ�
wd1 <- unlist(wd1_list)
wd2 <- unlist(wd2_list)
z1 <- unlist(Z1_list)
z2 <- unlist(Z2_list)
sparse_data1 <- as(WX1, "CsparseMatrix")
sparse_data2 <- as(WX2, "CsparseMatrix")


##�C���f�b�N�X���쐬
user_list1 <- list()
user_list2 <- list()
page_list <- list()
item_list <- list()
for(i in 1:hh){
  user_list1[[i]] <- which(u_id1==i)
  user_list2[[i]] <- which(u_id2==i)
}
for(i in 1:page){page_list[[i]] <- which(wd1==i)}
for(i in 1:item){item_list[[i]] <- which(wd2==i)}


####�}���R�t�A�������e�J�����@��Multinomial Topic Regression Model�𐄒�####
##�P�ꂲ�Ƃɖޓx�ƕ��S�����v�Z����֐�
burden_fr <- function(theta, phi, wd, w, k){
  Bur <-  matrix(0, nrow=length(wd), ncol=k)   #���S�W���̊i�[�p
  for(j in 1:k){
    #���S�W�����v�Z
    Bi <- rep(theta[, j], w) * phi[j, wd]   #�ޓx
    Bur[, j] <- Bi   
  }
  Br <- Bur / rowSums(Bur)   #���S���̌v�Z
  bval <- list(Br=Br, Bur=Bur)
  return(bval)
}

##���W�b�g���f���̑ΐ��ޓx�֐�
loglike <- function(Z, Data, beta1, beta2, index){
  
  #���W�b�g�Ɗm���̌v�Z
  logit <- Data[, 1]*beta1[index, ] + Data[, 2]*beta2[index, ]
  Pr <- exp(logit) / rowSums(exp(logit))
  
  LLi <- rowSums(Z * log(Pr))
  LL <- sum(LLi)
  val <- list(LLi=LLi, LL=LL)
  return(val)
}


##�A���S���Y���̐ݒ�
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##�p�����[�^�̐^�l
alpha1 <- alphat01
alpha2 <- alphat02
theta1 <- thetat1
theta2 <- thetat2
oldcov <- diag(dataset*(k-1))
diag(oldcov) <- c(rep(4.0, k-1), rep(1.0, k-1)) 
inv_cov <- solve(oldcov)
phi <- phit
gamma <- gammat


##�����l��ݒ�
alpha1 <- cbind(mvrnorm(hh, rep(0, k-1), diag(3.0, k-1)), 0)
alpha2 <- cbind(mvrnorm(hh, rep(0, k-1), diag(0.75, k-1)), 0)
logit <- Data[, 1]*alpha1[index_logit, ] + Data[, 2]*alpha2[index_logit, ]   #���W�b�g�̌v�Z
theta1 <- exp(logit[index_data[, 1], ]) / rowSums(exp(logit[index_data[, 1], ]))
theta2 <- exp(logit[index_data[, 2], ]) / rowSums(exp(logit[index_data[, 2], ]))
phi0 <- extraDistr::rdirichlet(k, colSums(WX1)/sum(WX1)*10) + 0.0001
phi <- phi0 / rowSums(phi0)
gamma0 <- extraDistr::rdirichlet(k, colSums(WX2)/sum(WX2)*10) + 0.0001
gamma <- gamma0 / rowSums(gamma0)

##���O���z�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z
alpha01 <- 0.5
alpha02 <- 0.5

#�ϗʌ��ʂ̎��O���z
nu <- k-1   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(rep(1, (dataset*k)-dataset))


##�p�����[�^�̊i�[�p�z��
ALPHA1 <- array(0, dim=c(hh, k, R/keep))
ALPHA2 <- array(0, dim=c(hh, k, R/keep))
SIGMA <- matrix(0, nrow=R/keep, ncol=(dataset*k)-dataset) 
PHI <- array(0, dim=c(k, page, R/keep))
GAMMA <- array(0, dim=c(k, item, R/keep))
SEG1 <- matrix(0, nrow=f1, ncol=k)
SEG2 <- matrix(0, nrow=f2, ncol=k)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"

#�ΐ��ޓx�̊�l
LLst <- sum(WX1 %*% log(colSums(WX1)/f1) + WX2 %*% log(colSums(WX2)/f2))




####�M�u�X�T���v�����O��HTM���f���̃p�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�f�[�^�Z�b�g���ƂɃg�s�b�N�𐶐�
  #�g�s�b�N�ޓx���v�Z
  page_par <- burden_fr(theta1, phi, wd1, w1, k)
  item_par <- burden_fr(theta2, gamma, wd2, w2, k)
  page_rate <- page_par$Br
  item_rate <- item_par$Br
  
  #�g�s�b�N�̊����m������g�s�b�N�𐶐�
  Zi1 <- rmnom(f1, 1, page_rate)
  Zi2 <- rmnom(f2, 1, item_rate)
  z1_vec <- as.numeric(Zi1 %*% 1:k)
  z2_vec <- as.numeric(Zi2 %*% 1:k)
  
  
  ##�o�����̕��z��phi�����gamma���X�V
  #�f�B�N�������z�̃p�����[�^���v�Z
  vf0 <- matrix(0, nrow=k, ncol=page)
  gf0 <- matrix(0, nrow=k, ncol=item)
  for(j in 1:page){
    vf0[, j] <- colSums(Zi1[page_list[[j]], , drop=FALSE])
  }
  for(j in 1:item){
    gf0[, j] <- colSums(Zi2[item_list[[j]], , drop=FALSE])
  }
  vf <- vf0 + alpha01
  gf <- gf0 + alpha02
  
  #�f�B�N�������z����p�����[�^���T���v�����O
  phi <- extraDistr::rdirichlet(k, vf)
  gamma <- extraDistr::rdirichlet(k, gf)
  
  
  ##�K�w�x�C�Y�������W�b�g���f���Ńg�s�b�N���z�̃p�����[�^���T���v�����O
  #�g�s�b�N�������烍�W�b�g���f���̉����ϐ��𐶐�
  item_sum <- page_sum <- matrix(0, nrow=hh, ncol=k)
  for(i in 1:hh){
    page_sum[i, ] <- rep(1, w1[i]) %*% page_rate[user_list1[[i]], ]
    item_sum[i, ] <- rep(1, w2[i]) %*% item_rate[user_list2[[i]], ]
  }
  Z <- rbind(page_sum, item_sum)
  
  ##MH�@�ŉ�A�p�����[�^���T���v�����O
  #�V�����p�����[�^���T���v�����O
  alphad1 <- alpha1; alphad2 <- alpha2
  alphan1 <- alphad1 + cbind(mvrnorm(hh, rep(0, k-1), diag(0.03, k-1)), 0)
  alphan2 <- alphad2 + cbind(mvrnorm(hh, rep(0, k-1), diag(0.01, k-1)), 0)
  er_new <- cbind(alphan1[, -k], alphan2[, -k]) - 0
  er_old <- cbind(alphad1[, -k], alphad2[, -k]) - 0
  
  #�ΐ��ޓx�Ƒΐ����O���z�̌v�Z
  lognew0 <- loglike(Z, Data, alphan1, alphan2, index_logit)$LLi
  logold0 <- loglike(Z, Data, alphad1, alphad2, index_logit)$LLi
  logpnew <- -0.5 * rowSums(er_new %*% inv_cov * er_new)
  logpold <- -0.5 * rowSums(er_old %*% inv_cov * er_old)
  
  #ID�ʂɑΐ��ޓx�̘a�����
  lognew <- lognew0[index_data[, 1]] + lognew0[index_data[, 2]]
  logold <- logold0[index_data[, 1]] + logold0[index_data[, 2]]
  
  #MH�T���v�����O
  rand <- runif(hh)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew + logpnew - logold - logpold)   #�̑𗦂��v�Z
  tau <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #tau�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((tau >= rand)*1 + (tau < rand)*0), nrow=hh, ncol=k)
  alpha1 <- flag*alphan1 + (1-flag)*alphad1   #alpha��rand�������Ă�����̑�
  alpha2 <- flag*alphan2 + (1-flag)*alphad2
  alpha <- cbind(alpha1[, -k], alpha2[, -k])
  
  #�������W�b�g���f���Ŋm���ɕϊ�
  theta1 <- exp(alpha1) / rowSums(exp(alpha1))
  theta2 <- exp(alpha1+alpha2) / rowSums(exp(alpha1+alpha2))
  
  ##�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^
  
  V_par <- V + t(alpha) %*% alpha
  Sn <- nu + hh
  
  #�t�E�B�V���[�g���z���番�U�����U�s��𔭐�
  oldcov <- rwishart(Sn, solve(V_par))$IW
  inv_cov <- solve(oldcov)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    ALPHA1[, , mkeep] <- alpha1
    ALPHA2[, , mkeep] <- alpha2
    SIGMA[mkeep, ] <- diag(oldcov)
    PHI[, , mkeep] <- phi
    GAMMA[, , mkeep] <- gamma
    
    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(mkeep >= burnin & rp%%keep==0){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
    }
    
    #�T���v�����O���ʂ��m�F
    if(rp%%disp==0){
      print(rp)
      print(mean(tau))
      print(round(Z[c(1, 2001), ], 3))
      print(c(sum(log(rowSums(page_par$Bur)))+sum(log(rowSums(item_par$Bur))), LLst))
      print(round(rbind(diag(oldcov)[5:24], c(rep(5, 10), rep(1.0, 10))), 3))
      print(round(cbind(phi[, 1:10], phit[, 1:10]), 3))
    }
  }
}
var(alpha)

####�T���v�����O���ʂ̉����Ɨv��####
burnin <- 2000/keep   #�o�[���C������
RS <- R/keep

##�T���v�����O���ʂ̉���
#HMM�̏������z�̃T���v�����O����
matplot(THETA1, type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

matplot(SIGMA, type="l")

#HMM�̃p�����[�^�̃T���v�����O����
matplot(t(ALPHA1[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA2[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA2[15, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

#�����̃g�s�b�N���z�̃T���v�����O����
matplot(t(THETA3[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA3[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA3[15, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

#�P��̏o���m���̃T���v�����O����
matplot(t(PHI[, 1, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N1�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[, 100, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N2�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[, 200, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N2�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[, 300, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N2�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[, 400, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N3�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[, 500, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N4�̒P��̏o�����̃T���v�����O����")


##�T���v�����O���ʂ̗v�񐄒��
#�g�s�b�N���z�̎��㐄���
topic_mu <- apply(THETA[, , burnin:(R/keep)], c(1, 2), mean)   #�g�s�b�N���z�̎��㕽��
round(cbind(topic_mu, thetat), 3)
round(topic_sd <- apply(THETA[, , burnin:(R/keep)], c(1, 2), sd), 3)   #�g�s�b�N���z�̎���W���΍�

#�P��o���m���̎��㐄���
word_mu <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #�P��̏o�����̎��㕽��
word <- round(t(rbind(word_mu, phit)), 3)
colnames(word) <- 1:ncol(word)
word

##�g�s�b�N�̎��㕪�z�̗v��
round(cbind(z1, seg1_mu <- SEG1 / length(burnin:RS)), 3)
round(cbind(z2, seg2_mu <- SEG2 / rowSums(SEG2)), 3)






