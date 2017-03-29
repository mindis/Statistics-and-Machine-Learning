#####���q����#####
library(polycor)
library(psych)
library(MASS)
library(reshape2)
library(plyr)
####�������������f�[�^�𔭐�####
#set.seed(489)
n <- 1000   #�T���v����
k <- 7   #�ϐ���
ff <- 3   #���q��
mu <- rnorm(7, 5, 5)   #���σx�N�g��
A <- matrix(runif(k*ff, 0, 0.6), k, ff)   #���q���׍s�� 
F <- matrix(rnorm(n*ff, 0, 1), n, ff)   #���ʈ��q�s��

##�����������ݒ肵�āA�ϑ��ϐ����쐬
MU <- matrix(mu, n, k, byrow=T)
X <- F %*% t(A) + matrix(rnorm(n*k, 0, 0.4), n, k, byrow=T)   #�ϑ��ϐ�
A

#�����������f�[�^�̗v��
summary(X)
colMeans(X)
cor(X)
plot(as.data.frame(X), col=4)


####���ϗʐ��K���z����̃f�[�^����####
val <- 7   #�����ϐ��̐�
n <- 1000   #�T���v����

##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

#�Q���Ƃ̑��֍s����쐬(�Q�ł��ׂē���)
corM <- corrM(col=val, lower=-0.3, upper=0.8)
eigen(corM)

##���֍s�񂩂番�U�����U�s����쐬����֐����`
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
    diag(c) <- m   #�Ίp�s������̕��U�ɖ߂�
  }
  diag(c) <- m
  cc <- c * corM
  #�ŗL�l�����ŋ����I�ɐ���l�s��ɏC������
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

#���U�����U�s����쐬(�Q�ł��ׂē���)
Sigma <- covmatrix(col=val, corM=corM, lower=1, upper=3)

#�Q���Ƃ̕ϐ��̕��ς��쐬
mu <- c(rnorm(val, 3, 2))

##���ϗʐ��K���z����̗����𔭐�������
val; n
X <- mvrnorm(n=n, mu, Sigma$covariance)

##�f�[�^��v��
round(colMeans(X), 2)   #�ϐ����Ƃ̕���
round(var(X), 2)   #���U�����U�s��
round(cor(X), 2)   #���֍s��
summary(X)   #�v��

##�U�z�}�s��̍쐬
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#�U�z�}�s��
pairs(as.data.frame(X[, 1:7]), panel=panel.smooth, bg="lightblue", diag.panel=panel.hist,
      upper.panel=panel.cor)



##�Ŗސ���ɂ����q���׍s��̐���
##EM�A���S���Y��
#�����l�̐ݒ�
R <- var(X)
Aold <- matrix(rep(0.5, k*ff), k, ff)
Aold[1, 2:3] <- 0
Aold[2, 3] <- 0
D2old <- diag(R - Aold %*% t(Aold))
D2old <- diag(D2old)
I <- diag(1, ff)

#EM�A���S���Y���̐ݒ�
max.iter <- 1000   #�ő�J��Ԃ���
iter <- 1   #�J��Ԃ��J�E���^�[
tol <- 10^(-1)   #����l�̕ω��̋��e�x
S.zz <- R   #S.zz�̏����t�����Ғl�̓f�[�^�̕��U�����U�s��

##EM�A���S���Y���ɂ�鐄��
while(iter < max.iter){
  #E�X�e�b�v
  Sigma <- Aold %*% t(Aold) + D2old
  delta <- t(Aold) %*% solve(Sigma)
  S.zf <- S.zz %*% t(delta)
  S.ff <- delta %*% S.zz %*% t(delta) + (I - delta %*% Aold)
  #M�X�e�b�v
  Anew <- S.zf %*% solve(S.ff)
  D2new <- diag(S.zz - S.zf %*% solve(S.ff) %*% t(S.zf))
  #����l�̕ω�
  diff <- max(abs(Anew - Aold), abs(D2new - diag(D2old)))
  if(diff < tol) break;
  Aold <- Anew
  D2old <- diag(D2new)
  iter <- iter+1
  print(diff)
}

#���ʂ��m�F
delta <- round(1-D2new, 2)
Anew
promax(Anew)

##�֐���p���Đ���
fa(X, nfactors=3, rotate="promax", fm="ml", cor="cov")
promax(Anew)
help(fa)