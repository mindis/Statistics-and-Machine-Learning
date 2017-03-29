#####k-means�@######
library(MASS)
library(dplyr)
library(reshape2)
####�f�[�^�̔���#####
#set.seed(493)
k <- 4   #�N���X�^�[��
val <- 6   #�����ϐ��̐�
n <- 500   #�Z�O�����g���Ƃ̃f�[�^��

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
corM <- corrM(col=val, lower=-0.2, upper=0.2)
eigen(corM)

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
Sigma1 <- covmatrix(col=val, corM=corM, lower=25, upper=35)
Sigma2 <- covmatrix(col=val, corM=corM, lower=30, upper=40)
Sigma3 <- covmatrix(col=val, corM=corM, lower=20, upper=30)
Sigma4 <- covmatrix(col=val, corM=corM, lower=27, upper=38)
Sigma <- list(Sigma1$covariance, Sigma2$covariance, Sigma3$covariance, Sigma3$covariance)

#�Q���Ƃ̕ϐ��̕��ς��쐬
mu1 <- c(rnorm(val, 14, 10))
mu2 <- c(rnorm(val, 19, 10))
mu3 <- c(rnorm(val, 9, 8))
mu4 <- c(rnorm(val, 23, 13))
mu <- list(mu1, mu2, mu3, mu4)

##���ϗʐ��K���z����̗����𔭐�������
k; n
X <- matrix(0, 0, val)
for(kk in 1:k){
  xx <- mvrnorm(n=n, mu[[kk]], Sigma[[kk]])
  X <- rbind(X, xx)
}

#�N���X�^�[�̃x�N�g�����쐬
y <- rep(1:4, rep(n, 4))
table(y)

#�f�[�^���������ăf�[�^��v��
YX <- as.data.frame(cbind(y, X))
by(YX, YX[, 1] , function(x) summary(x))   #�Z�O�����g���Ƃ̗v��֐�
by(YX[, 2:7], YX[, 1], colMeans)   #�Z�O�����g���Ƃ̕���
by(YX[, 2:7], YX[, 1] , var)   #�Z�O�����g���Ƃ̕��U�����U�s��

####k-means�@�ŃZ�O�����g�ɕ�����####
##����������ݒ�
#�ϐ���W�������čs�̍��v�l�������ɕ��ׂď��������Ƃ���
Xscale <- scale(YX[, -1])
Xsums <- cbind(YX, apply(Xscale, 1, sum))
names(Xsums)[8] <- "Xscale"
Xseg <- Xsums[order(Xsums[, 8], decreasing = TRUE), c(-1, -8)]
trueSeg <- Xsums[order(Xsums[, 8], decreasing = TRUE), 1]   #�^�̃Z�O�����g�����Ԃɕ��ׂĂ���
ns <- dim(Xseg)[1]/k   #1�Z�O�����g������̕�����

#�Z�O�����g���Ƃɕ������ĕ��σx�N�g�������߂�
X1mean <- colMeans(Xseg[1:ns, ])
X2mean <- colMeans(Xseg[(ns+1):(ns*2), ])
X3mean <- colMeans(Xseg[(ns*2+1):(ns*3), ])
X4mean <- colMeans(Xseg[(ns+1):(ns*4), ]) 
Xmean <- list(X1mean, X2mean, X3mean, X4mean)

##���ꂼ��̃T���v���ɑ΂��ăZ�O�����g���Ƃ̃��[�N���b�h���������߂āA
##���̋������ŏ��ƂȂ�Z�O�����g�ɃT���v�������������������������
segsign <- rep(1:4, rep(500, 4))
Sw <- numeric()
for(m in 1:k){
  ss <- apply(Xseg[segsign==m, ], 1, function(x) sum((x-Xmean[[m]])^2))
  (ssm <- sum(ss)/table(segsign)[[m]])
  Sw <- c(Sw, ssm)  
}
(Sold <- sum(Sw))   #�N���X�^�[�����a�̍��v���v�Z
diff <- 100   #�����a�̍��̏����l
tol <- 1   #��~�


###k-means�@�̃A���S���Y��
while(abs(diff) >= tol){
  ##�����Z�O�����g�����肷��A���S���Y��
  Xdis <- matrix(0, nrow=nrow(Xseg), ncol=k)
  for(kk in 1:k){
    xxd <- apply(Xseg, 1, function(x) sum((x-Xmean[[kk]])^2))
    Xdis[, kk] <- xxd
  }
  #��������Z�O�����g�̐ݒ�
  seg <- apply(Xdis, 1, which.min)   #�����a���ŏ��̃Z�O�����g�ɏ���������
  segN <- table(seg)
  YXseg <- cbind(seg, Xseg)
  
  ##�Z�O�����g�����̕]�������߂�A���S���Y��
  #�Z�O�����g���Ƃɕ��σx�N�g�������߂�
  Xmean <- list()
  for(m in 1:k){
    xm <- apply(Xseg[seg==m, ], 2, mean)
    Xmean[[m]] <- xm
  }
  
  #�N���X�^�[�������a���v�Z
  Sw <- numeric()
  for(m in 1:k){
    ss <- apply(Xseg[seg==m, ], 1, function(x) sum((x-Xmean[[m]])^2))
    (ssm <- sum(ss)/segN[m])
    Sw <- c(Sw, ssm)  
  }
  S <- sum(Sw)   #�N���X�^�[�����a�̍��v���v�Z
  diff <- abs(Sold - S)   #�ȑO�̃N���X�^�[�����a�Ƃ̍��̐�Βl���v�Z
  Sold <- S
  print(S)
}

##���ʂ̗v��
resultX <- cbind(trueSeg, seg, Xseg)   #�^�̃Z�O�����g�A���肳�ꂽ�Z�O�����g�A�f�[�^������
sortlist <- order(resultX[, 2])
resultX <- resultX[sortlist, ]
rownames(resultX) <- c(1:nrow(resultX))   #�s�ԍ���U�蒼�� 
table(resultX[, 2], resultX[, 1])   #�^�̃Z�O�����g�Ɛ��肳�ꂽ�Z�O�����g�̌딻�ʕ\

#���肳�ꂽ�Z�O�����g���Ƃ̗v��
by(resultX[, 3:8], resultX[, 2] , function(x) summary(x))   #���肳�ꂽ�Z�O�����g���Ƃ̗v��֐�
by(resultX[, 3:8], resultX[, 2], colMeans)   #�Z�O�����g���Ƃ̕���
by(resultX[, 3:8], resultX[, 2] , var)   #�Z�O�����g���Ƃ̕��U�����U�s��
plot(resultX[, 3:8], col=resultX[, 2])

#�^�̃Z�O�����g���Ƃ̗v��
by(resultX[, 3:8], resultX[, 1] , function(x) summary(x))   #���肳�ꂽ�Z�O�����g���Ƃ̗v��֐�
by(resultX[, 3:8], resultX[, 1], colMeans)   #�Z�O�����g���Ƃ̕���
by(resultX[, 3:8], resultX[, 1] , var)   #�Z�O�����g���Ƃ̕��U�����U�s��
plot(resultX[, 3:8], col=resultX[, 1])

####�֐����g��####
res <- kmeans(x=resultX[, 3:8], centers=k)
resX <- cbind(res$cluster, resultX)
names(resX)[1] <- "cluster_f"
table(resX$seg, resX$cluster_f)   #�����kmeans�@��R��kmeans�@���r

