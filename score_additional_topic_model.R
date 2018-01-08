#####���t����g�s�b�N���f��####
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


####�f�[�^�̔���####
#set.seed(423943)
#�f�[�^�̐ݒ�
s <- 5   #���t��
k1 <- 3   #���t���Ƃ̃g�s�b�N��
k2 <- 10   #�����W���S�̂̃g�s�b�N��
d <- 2500   #������
v1 <- 300   #���t�Ɋ֌W�̂����b��
v2 <- 100   #���t�Ɋ֌W�̂Ȃ���b��
v <- v1 + v2   #��b��
w <- rpois(d, rgamma(d, 55, 0.50))   #1����������̒P�ꐔ
f <- sum(w)

#�C���f�b�N�X���쐬
id_vec <- rep(1:d, w)
index_hh <- list()
for(i in 1:d){
  index_hh[[i]] <- which(id_vec==i)
}

#���t�f�[�^�𐶐�
pr <- runif(s, 1, 5)
Y <- rmnom(d, 1, pr)
y <- as.numeric(Y %*% 1:s)
y_vec <- as.numeric(y)[id_vec] 
index_y <- list()
for(j in 1:s) {index_y[[j]] <- which(y==j)}
y_freq <- as.numeric(table(y))


##�p�����[�^�̐ݒ�
#���t���Ƃ̃g�s�b�N�𐶐�
alpha0 <- rep(0.3, k1)   #���t�f�[�^�̕����̃f�B���N�����O���z�̃p�����[�^
alpha1 <- list()
for(j in 1:s){
  alpha1[[j]] <- c(rep(0.4, v1), rep(0.005, v2))   #���t�Ɋ֌W�̂���P��̃f�B���N�����O���z�̃p�����[�^
}
alpha2 <- c(rep(0.1, v1), rep(5, v2))   #���t�Ɋ֌W�̂Ȃ��P��̃f�B���N�����O���z�̃p�����[�^

#�f�B���N�������̔���
thetat <- theta <- rdirichlet(d, alpha0)   #�����̃g�s�b�N���z���f�B���N���������甭��
phit <- phi <- list()
for(j in 1:s) {phit[[j]] <- phi[[j]] <- rdirichlet(k1, alpha1[[j]])}   #�]�_�Ɋ֌W�̂���P�ꕪ�z���f�B���N���������甭��
gammat <- gamma <- rdirichlet(1, alpha2)   #�]�_�Ɋ֌W�̂Ȃ��P�ꕪ�z���f�B���N���������甭��
betat <- beta <- rbeta(sum(f), 15, 15)   #�P�ꂪ���t�Ɗ֘A���邩�ǂ����̃p�����[�^


##�������z�̗�������f�[�^�𔭐�
WX <- matrix(0, nrow=d, ncol=v)
x_list <- list()
x <- rep(0, f)
Z <- list()
index_v1 <- 1:v1
index_v2 <- (v1+1):v

for(i in 1:d){
  
  #�����̃g�s�b�N�𐶐�
  z <- rmnom(w[i], 1, theta[i, ])   #�����̃g�s�b�N���z�𔭐�
  z_vec <- z %*% c(1:k1)   #�g�s�b�N�������x�N�g����
  
  #��ʌꂩ�ǂ����𐶐�
  x_list[[i]] <- rbinom(w[i], 1, beta[index_hh[[i]]])
  
  phi[[y[i]]][z_vec[x_list[[i]]==1], ]
  phi
  #���������g�s�b�N����P��𐶐�
  wn <- rmnom(sum(x_list[[i]]), 1, phi[[y[i]]][z_vec[x_list[[i]]==1], ])   #�����̃g�s�b�N����P��𐶐�
  an <- rmnom(sum(1-x_list[[i]]), 1, gammat)
  wdn <- colSums(wn) + colSums(an)   #�P�ꂲ�Ƃɍ��v����1�s�ɂ܂Ƃ߂�
  WX[i, ] <- wdn
  Z[[i]] <- z
  x[index_hh[[i]]] <- x_list[[i]]
  print(i)
}

####�g�s�b�N���f������̂��߂̃f�[�^�Ɗ֐��̏���####
##���ꂼ��̕������̒P��̏o������ѕ⏕���̏o�����x�N�g���ɕ��ׂ�
##�f�[�^����pID���쐬
ID_list <- list()
wd_list <- list()

#�������Ƃɕ���ID����ђP��ID���쐬
for(i in 1:nrow(WX)){
  print(i)
  
  #�P���ID�x�N�g�����쐬
  ID_list[[i]] <- rep(i, w[i])
  num1 <- (WX[i, ] > 0) * (1:v)
  num2 <- which(num1 > 0)
  W1 <- WX[i, (WX[i, ] > 0)]
  number <- rep(num2, W1)
  wd_list[[i]] <- number
}

#���X�g���x�N�g���ɕϊ�
ID_d <- unlist(ID_list)
wd <- unlist(wd_list)

##�C���f�b�N�X���쐬
doc_list <- list()
word_list <- list()
y_list <- list()
index_s <- list()
words_list <- list()
docs_list <- list()
docs_vec <- list()
wd0 <- list()
w0 <- list()

for(i in 1:length(unique(ID_d))) {doc_list[[i]] <- which(ID_d==i)}
for(i in 1:length(unique(wd))) {word_list[[i]] <- which(wd==i)}
for(j in 1:s) {y_list[[j]] <- which(y==j)}
for(j in 1:s) {index_s[[j]] <- which(y_vec==j)}

for(i in 1:s){
  words_list[[i]] <- list()
  wds <- wd[index_s[[i]]]
  
  for(j in 1:length(unique(wd))){
    words_list[[i]][[j]] <- which(wds==j)
  }
}

for(j in 1:s){
  wd0[[j]] <- wd[index_s[[j]]]
  w0[[j]] <- w[index_y[[j]]]
  docs_list[[j]] <- list()
  dcs <- ID_d[index_s[[j]]]
  vec <- unique(dcs)
  for(i in 1:length(vec)){
    docs_list[[j]][[i]] <- which(dcs==vec[i])
    docs_vec[[j]] <- vec
  }
}


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
burnin <- 1000/keep

##���O���z�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z
alpha01 <- 1.0
alpha02 <- 0.5


##�p�����[�^�̏����l
#�g�s�b�N���z�̏����l
theta <- extraDistr::rdirichlet(d, rep(1, k1))   #�����g�s�b�N�̃p�����[�^�̏����l


#��ʌ�̒P�ꕪ�z�̏����l
inv_idf <- colSums(WX > 0)/d
gamma <- inv_idf / sum(inv_idf)

#���t�֘A�P�ꕪ�z�̏����l
phi <- list()
for(j in 1:s) {
  M <- colSums(WX[y==j, ])
  phi[[j]] <- extraDistr::rdirichlet(k1, (M+1)*log(1/inv_idf))
}


#�������̏����l
rd <- matrix(0, nrow=f, ncol=2)
rd[, 1] <- 0.5
rd[, 2] <- 0.5


##�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(d, k1, R/keep))
PHI1 <- PHI2 <- PHI3 <- PHI4 <- PHI5 <- array(0, dim=c(k1, v, R/keep))
GAMMA <- matrix(0, nrow=R/keep, ncol=v)
SEG1 <- matrix(0, nrow=f, ncol=k1)
SEG2 <- rep(0, nrow=f)
storage.mode(SEG1) <- "integer"
gc(); gc()

##MCMC����p�z��
wsum0 <- matrix(0, nrow=d, ncol=k1)
vf0 <- matrix(0, nrow=k1, ncol=v)
zf0 <- rep(0, v)
switching_rate <- rep(0, f)
Zi1_full <- matrix(0, nrow=f, ncol=k1)
Zi2_zeros <- rep(0, f)
vec1 <- 1/1:k1


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  #�p�����[�^�̊i�[�p�z����X�V
  Brs <- list()
  Zi1 <- list()
  Zi2 <- list()
  
  for(j in 1:s){
    ##���t�P�ʂł̃g�s�b�N���z���T���v�����O
    #�P�ꂲ�ƂɃg�s�b�N���z�̃p�����[�^�𐄒�
    theta0 <- theta[index_y[[j]], ]
    Brs[[j]] <- burden_fr(theta0, phi[[j]], wd0[[j]], w0[[j]], k1)
    word_rate <- Brs[[j]]$Br
    
    #�������z���g�s�b�N���T���v�����O
    n <- nrow(word_rate)
    rand1 <- matrix(runif(n), nrow=n, ncol=k1)
    user_cumsums <- rowCumsums(word_rate)
    Zi1[[j]] <- ((k1+1) - (user_cumsums > rand1) %*% rep(1, k1)) %*% vec1   #�g�s�b�N���T���v�����O
    Zi1[[j]][Zi1[[j]]!=1] <- 0
    Zi1_full[index_s[[j]], ] <- Zi1[[j]]
    
    
    ##�P�ꂲ�Ƃɋ��t�g�s�b�N�Ɗ֌W�����邩�ǂ������T���v�����O
    #�񍀕��z�̃p�����[�^�𐄒�
    LLz1 <- rd[index_s[[j]], 1] * rowSums(Brs[[j]]$Bur)
    LLz2 <- rd[index_s[[j]], 2] * gamma[wd0[[j]]]
    switching_rate[index_s[[j]]] <- LLz1 / (LLz1+LLz2)
    
    #�񍀕��z����֌W�L���𐶐�
    Zi2[[j]] <- rbinom(length(index_s[[j]]), 1, switching_rate[index_s[[j]]])
  }
  
  ##���t���ƂɃp�����[�^���T���v�����O
  for(j in 1:s){
    #�����̃g�s�b�N�̃p�����[�^�𐄒�
    zi <- Zi1[[j]] * matrix(Zi2[[j]], nrow=length(index_s[[j]]), ncol=k1)
    Zi2_zeros[index_s[[j]]] <- 1-Zi2[[j]]
    
    for(i in 1:y_freq[j]){
      wsum0[index_y[[j]][i], ] <- colSums(zi[docs_list[[j]][[i]], ])
    }
    
    #�P�ꕪ�z�̃p�����[�^�𐄒�
    for(l in 1:v){
      vf0[, l] <- colSums(zi[words_list[[j]][[l]], ,drop=FALSE])
    }
    #�f�B�N�������z����P�ꕪ�z���T���v�����O
    vf <- vf0 + alpha02
    phi[[j]] <- extraDistr::rdirichlet(k1, vf)
  }
  
  #�f�B�N�������z����g�s�b�N���z���T���v�����O
  wsum <- wsum0 + alpha01
  theta <- extraDistr::rdirichlet(d, wsum)
  
  
  ##��ʌ�̃p�����[�^���T���v�����O
  #��ʌ�p�����[�^gamma���T���v�����O
  for(j in 1:v){
    zf0[j] <- sum(Zi2_zeros[word_list[[j]]])
  }
  zf <- zf0 + alpha02
  gamma <- extraDistr::rdirichlet(1, zf)
  
  #���������T���v�����O
  Zi2_ones <- 1-Zi2_zeros
  for(i in 1:d){
    rd[i, 1] <- mean(Zi2_ones[doc_list[[i]]])
  }
  rd[, 2] <- 1-rd[, 1]
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    PHI1[, , mkeep] <- phi[[1]]
    PHI1[, , mkeep] <- phi[[2]]
    PHI1[, , mkeep] <- phi[[3]]
    PHI1[, , mkeep] <- phi[[4]]
    PHI1[, , mkeep] <- phi[[5]]
    GAMMA[mkeep, ] <- gamma
    
    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(rp%%keep==0 & rp >= burnin){
      SEG1 <- SEG1 + Zi1_full      
      SEG2 <- SEG2 + Zi2_ones
    }
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(round(cbind(theta[y==1, ][1:7, ], thetat[y==1, ][1:7, ], theta[y==2, ][1:7, ], thetat[y==2, ][1:7, ],
                      theta[y==3, ][1:7, ], thetat[y==3, ][1:7, ]), 3))
    print(round(cbind(phi[[1]][, 296:305], phit[[1]][, 296:305]), 3))
    print(round(cbind(phi[[2]][, 296:305], phit[[2]][, 296:305]), 3))
    print(round(cbind(phi[[3]][, 296:305], phit[[3]][, 296:305]), 3))
    print(round(rbind(gamma[290:309], gammat[290:309]), 3))
  }
}



####�T���v�����O���ʂ̉����Ɨv��####
burnin <- 1000/keep   #�o�[���C������
RS <- R/keep

##�T���v�����O���ʂ̉���
#�����̃g�s�b�N���z�̃T���v�����O����
matplot(t(THETA[1, , ]), type="l", ylab="�p�����[�^", main="����1�̃g�s�b�N���z�̃T���v�����O����")
matplot(t(THETA[2, , ]), type="l", ylab="�p�����[�^", main="����2�̃g�s�b�N���z�̃T���v�����O����")
matplot(t(THETA[3, , ]), type="l", ylab="�p�����[�^", main="����3�̃g�s�b�N���z�̃T���v�����O����")
matplot(t(THETA[4, , ]), type="l", ylab="�p�����[�^", main="����4�̃g�s�b�N���z�̃T���v�����O����")

#�P��̏o���m���̃T���v�����O����
matplot(t(PHI[1, 1:10, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N1�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[2, 11:20, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N2�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[3, 21:30, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N3�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[4, 31:40, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N4�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[5, 41:50, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N5�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[6, 51:60, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N6�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[7, 61:70, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N7�̒P��̏o�����̃T���v�����O����")
matplot(t(PHI[8, 71:80, ]), type="l", ylab="�p�����[�^", main="�g�s�b�N8�̒P��̏o�����̃T���v�����O����")

##�T���v�����O���ʂ̗v�񐄒��
#�g�s�b�N���z�̎��㐄���
topic_mu <- apply(THETA[, , burnin:(R/keep)], c(1, 2), mean)   #�g�s�b�N���z�̎��㕽��
round(cbind(topic_mu, thetat), 3)
round(topic_sd <- apply(THETA[, , burnin:(R/keep)], c(1, 2), sd), 3)   #�g�s�b�N���z�̎���W���΍�

#�P��o���m���̎��㐄���
word_mu <- apply(PHI[, , burnin:(R/keep)], c(1, 2), mean)   #�P��̏o�����̎��㕽��
round(rbind(word_mu, phit)[, 1:50], 3)