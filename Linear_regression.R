#####���`��A���f��#####
library(MASS)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

####�f�[�^�̔���####
k <- 4   #���ϗʉ�A���f���̐�
col <- 15   #�ϐ���
cont <- 7   #�A���ϐ���
bin <- 3   #��l�ϐ���
multi <- 5   #���l�ϐ��̃J�e�S����
n <- 10000   #�T���v����

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

#���ϗʉ�A���f���̑��֍s����쐬
corM <- corrM(col=k, lower=-0.5, upper=0.6)

##���֍s�񂩂番�U�����U�s����쐬����֐����`
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
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

#���U�����U�s����쐬
Sigma <- covmatrix(col=k, corM=corM, lower=15, upper=23)
Cov <- Sigma$covariance

##�����ϐ��̔���
##7�̘A���ϐ��A3�̓�l�ϐ��A1�̑��l�ϐ�(�J�e�S����5)
#�A���ϐ�
X_cont <- matrix(runif(n*cont, 0, 1), nrow=n, ncol=cont, byrow=T)

#��l�ϐ�
X_bin <- matrix(0, n, bin)
pr_b <- runif(3, 0.3, 0.7)
for(i in 1:bin){
  bin <- rbinom(n, 1, pr_b[i])
  X_bin[, i] <- bin
}

#���l�ϐ�
pr_mm <- runif(multi, 0.1, 1.0)
pr_m <- pr_mm/sum(pr_mm)
X_multim <- t(rmultinom(n, 1, pr_m))
X_multi <- X_multim[, -which.min(colSums(X_multim))]   #�璷�ȕϐ����폜
colSums(X_multim); colSums(X_multi) 

##�f�[�^�̌���
X <- data.frame(cont=round(X_cont, 3), bin=X_bin, multi=X_multi)
XM <- as.matrix(X)
summary(X)

##��A�W���̐ݒ�
beta1 <- runif(col-1, -3, 8)
beta01 <- runif(1, 12, 16)
beta2 <- runif(col-1, -3.5, 9)
beta02 <- runif(1, 13, 18)
beta3 <- runif(col-1, -2.8, 10.5)
beta03 <- runif(1, 10, 20)
beta4 <- runif(col-1, -4, 11)
beta04 <- runif(1, 12, 22)

#��A���f���̕��ύ\����ݒ�
z1 <- XM %*% beta1 + beta01
z2 <- XM %*% beta2 + beta02
z3 <- XM %*% beta3 + beta03
z4 <- XM %*% beta4 + beta04
summary(z1); summary(z2); summary(z3); summary(z4)

##���ϗʐ��K�������牞���ϐ��𔭐�������
Z <- cbind(z1, z2, z3, z4)
Y <- trunc(Z + mvrnorm(n=n, rep(0, 4), Cov))   #���ϗʐ��K�����ŉ����ϐ��𔭐�
Y <- ifelse(Y <= 0, 1, Y)   #Y��0�ȉ��Ȃ�1�ɒu��������

#�����ϐ��̗v��
summary(Y)
cov(Y)
cor(Y)

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

#pairs(as.data.frame(Y), panel=panel.smooth, bg="lightblue", diag.panel=panel.hist,
#      upper.panel=panel.cor)


####���`��A���f���Ő���####
##�e��A���f�����Ɨ��Ƃ��čŏ����@�Ő���
XM1 <- as.matrix(data.frame(inter=1, XM))   #�ؕЂ�����f�U�C���s��

betan <- matrix(0, k, col)
Sigman <- numeric()
for(i in 1:k){
   coef <- solve(t(XM1) %*% XM1) %*% t(XM1) %*% Y[, i]
   sig <- sum((Y[, i] - XM1 %*% coef)^2) / (n-1)
   betan[i, ] <- coef
   Sigman <- c(Sigman, sig)
}

##���茋�ʂƗv��
round(betan, 3)   #��A�W���̐��茋��
round(beta_true <- rbind(c(beta01, beta1), c(beta02, beta2), c(beta03, beta3), c(beta04, beta4)), 3)   #�^�̉�A�W��
round(Sigman, 3)   #���U�̐��茋��
round(diag(Cov), 3)   #�^�̕��U

#��A���f�����Ƃ̌덷
p1 <- trunc(XM1 %*% betan[1, ])
p2 <- trunc(XM1 %*% betan[2, ])
p3 <- trunc(XM1 %*% betan[3, ])
p4 <- trunc(XM1 %*% betan[4, ])

data.frame(Y=Y, p1, p2, p3, p4, e1=Y[, 1]-p1, e2=Y[, 2]-p2, e3=Y[, 3]-p3, e4=Y[, 4]-p4)

##���v�ʂ̌v�Z
#R2���v��
ytotal1 <- Y[, 1] - mean(Y[, 1])
ypred1 <- XM1 %*% betan[1, ] - mean(Y[, 1])
round(R2_1 <- sum(ypred1^2) / sum(ytotal1^2), 3) 

ytotal2 <- Y[, 2] - mean(Y[, 2])
ypred2 <- XM1 %*% betan[2, ] - mean(Y[, 2])
round(R2_2 <- sum(ypred2^2) / sum(ytotal2^2), 3)

ytotal3 <- Y[, 3] - mean(Y[, 3])
ypred3 <- XM1 %*% betan[3, ] - mean(Y[, 3])
round(R2_3 <- sum(ypred3^2) / sum(ytotal3^2), 3)

ytotal4 <- Y[, 4] - mean(Y[, 4])
ypred4 <- XM1 %*% betan[4, ] - mean(Y[, 4])
round(R2_4 <- sum(ypred4^2) / sum(ytotal4^2), 3)

#AIC
aic <- c()
for(i in 1:k){
  a <- n*(log(2*pi)+1) + n*log(Sigman[i]) + 2*(col+2)
  aic <- c(aic, a)
}
aic

##�֐����g���ĉ�A���f���𐄒�
res1 <- list()
for(i in 1:k){
  rm(YX)
  YX <- data.frame(Y=Y[, i], XM) 
  l <- lm(Y ~ ., data=YX)
  res1[[i]] <- l
}
summary(res1[[1]])
summary(res1[[2]])
summary(res1[[3]])
summary(res1[[4]])


####���ϗʉ�A���f���ɂ�鐄��####
##�ŏ����@�ő��ϗʉ�A���f���𐄒�
round(BETA <- solve(t(XM1) %*% XM1) %*% t(XM1) %*% Y, 3)   #���ϗʉ�A���f���̉�A�W��
round(SIGMA <- (t(Y - XM1 %*% BETA) %*% (Y - XM1 %*% BETA)) / (n-1), 3)   #���U�����U�s��̐���l

##���ʂ̊m�F�Ɠ��v��
#�^�̌W���Ƃ̔�r
#��A�W��
round(t(BETA), 3)   
round(beta_true, 3)

#���U�����U�s��
round(SIGMA, 3)
round(Cov, 3)

#���֍s��
round(cov2cor(SIGMA), 3)
round(cov2cor(Cov), 3)
round(cor(Y), 3)   #�ϑ��f�[�^�̑��֍s��

##AIC���v�Z����
#���ϗʐ��K���z�̖ޓx�֐����v�Z
dmv <- function(y, x, beta, s){
  LLo  <-  1/(sqrt((2 * pi)^nrow(s) * det(s))) * 
    exp(-t(as.matrix(y) - t(beta) %*% as.matrix(x)) %*%
          solve(as.matrix(s)) %*%
          (as.matrix(y) - t(beta) %*% as.matrix(x)) / 2)
  return(LLo)
}
Li <- apply(cbind(Y, XM1), 1, function(x) dmv(y=x[1:k], x=x[(k+1):(col+k)], beta=BETA, s=SIGMA))

LLs <- sum(log(Li))   #�ΐ��ޓx�֐����v�Z 
(AIC <- -2*LLs +2*(col*k+sum(1:4)+1))   #AIC���v�Z

##�Ɨ���A���f���Ƒ��ϗʉ�A���f����AIC�̔�r
aic; AIC
sum(aic); AIC

#�\���l�Ǝ����l�̔�r
fit <- trunc(XM1 %*% BETA)   #�\���l
data.frame(Y=Y, p=fit, e=Y-fit)   #�����ϐ��Ɨ\���l�Ƃ̔�r

#�Ɨ���A���f���Ƒ��ϗʉ�A���f���̓��\���덷�̔�r
#�Ɨ����f���̓��\���덷���v�Z
sq_error <- c()
for(i in 1:k){
 error <- sum((Y[, i] - XM1 %*% betan[i, ])^2)
 sq_error <-  c(sq_error, error)
}

#���ϗʃ��f���̓��\���덷���v�Z
SQ_error <- sum((Y - XM1 %*% BETA)^2)

#�덷���r
sq_error; SQ_error
sum(sq_error); SQ_error

data.frame(Y=Y, P=fit, E=Y-fit)   #���ϗʉ�A���f���̊ϑ��f�[�^���Ƃ̌덷
data.frame(y=Y, p1, p2, p3, p4, e1=Y[, 1]-p1, e2=Y[, 2]-p2, e3=Y[, 3]-p3, e4=Y[, 4]-p4)   #�Ɨ����f���̌덷
data.frame(E=Y-fit, e1=Y[, 1]-p1, e2=Y[, 2]-p2, e3=Y[, 3]-p3, e4=Y[, 4]-p4)   #�덷�̔�r

#�^�̉�A�W���Ƃ̌덷
round(t(BETA), 3)
round(betan, 3)
round(beta_true, 3)
round(abs(beta_true - t(BETA)), 3)   #���ϗʃ��f���̉�A�W���̌덷
round(abs(beta_true - betan), 3)   #�Ɨ����f���̉�A�W���̌덷
round(sum(abs(beta_true - t(BETA))), 3)   #���ϗʃ��f���̌덷�̍��v
round(sum(abs(beta_true - betan)), 3)   #�Ɨ����f���̌덷�̍��v

##�֐����g���đ��ϗʉ�A���f���𐄒�
res2 <- lm(Y ~ XM)
summary(res2)


