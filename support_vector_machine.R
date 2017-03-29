#####�T�|�[�g�x�N�^�[�}�V��#####
library(MASS)
library(kernlab)
library(mlbench)
library(dplyr)
library(reshape2)
library(quadprog)
####�f�[�^�̔���#####
k <- 2   #�Q�̐�
val <-12   #�����ϐ��̐�
n <- 500   #�Q���Ƃ̊w�K�Ɏg�����߂̃f�[�^��
nt <- 500   #�Q���Ƃ̃e�X�g�Ɏg�����߂̃f�[�^��

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

#���U�����U�s����쐬(�Q�ł��ׂē���)
Sigma1 <- covmatrix(col=val, corM=corM, lower=25, upper=35)
Sigma2 <- covmatrix(col=val, corM=corM, lower=30, upper=40)
Sigma <- list(Sigma1$covariance, Sigma2$covariance)
Sigma
#�Q���Ƃ̕ϐ��̕��ς��쐬
mu1 <- c(rnorm(val, 10, 10))
mu2 <- c(rnorm(val, 14, 10))
mu <- list(mu1, mu2)

##���ϗʐ��K���z����̗����𔭐�������
k; n; nt
X <- matrix(0, 0, val)
for(kk in 1:k){
  xx <- mvrnorm(n=n+nt, mu[[kk]], Sigma[[kk]])
  X <- rbind(X, xx)
}

#���t�f�[�^�̃x�N�g�����쐬
y <- rep(c(1, 0), rep(n+nt, 2))
table(y)

#�ꕔ�̕ϐ��𗣎U�ϐ��ɕϊ�
mu_rbind <- colMeans(rbind(mu[[1]], mu[[2]]))   #2�Q�̓�������
x2 <- t(apply(X[, 8:val], 1, function(x) ifelse(x > mu_rbind[8:val], 1, 0)))   #�������ς�荂���v�f��1�Ƃ���
mu_rbind[8:val]
head(cbind(X[, 8:val], x2), 20)
head(YX)

#�A���ϐ���W����
Xscale <- scale(X[, 1:7])

#�f�[�^������
YX <- as.data.frame(cbind(y, Xscale, x2))
by(YX, YX[, 1] , function(x) summary(x))   #�Q���Ƃ̗v��֐�
by(YX[, 2:8], YX[, 1] , function(x) cor(x))   #�Q���Ƃ̘A���ϐ��̑���
by(YX[, 9:val+1], YX[, 1], colMeans)   #�Q���Ƃ̗��U�ϐ��̏o����

#�A���ϐ��̐����ϐ����v���b�g
yplotcolor <- YX[, 1]
yplotcolor[yplotcolor==0] <- 2
plot(YX[, 2:8], col=yplotcolor)   #�A���ϐ��̎U�z�}
YX[YX[, 1]==0, 1] <- -1   #0�̌Q��-1�ɂ���

#�w�K�f�[�^�ƃe�X�g�f�[�^�ɕ�����
YXl <- YX[c(1:500, 1001:1500), ]
YXt <- YX[c(501:1000, 1501:2000), ]

####�T�|�[�g�x�N�^�[�}�V���̎���####
##���`�T�|�[�g�x�N�^�[�}�V��
#�ő剻����֐�
dim(YXl)
yy <- matrix(YXl[, 1], nrow=nrow(YXl), ncol=ncol(YXl[, -1])) 
xx <- yy * YXl[, -1]
Q <- as.matrix(xx[, -1]) %*% t(xx[, -1])
c <- matrix(rep(-1, nrow(YXl)))

#�������
sm <- 0.0001   #�\�t�g�}�[�W��
A <- t(YXl[, 1])
b <- 0
l <- matrix(rep(0, nrow(YXl)))
u <- matrix(rep(sm), nrow(YXl))
r <- 0

#�ʓ񎟌v���������
sv <- ipop(c, Q, A, b, l, u, r)
sv
dual(sv)

#beta�����߂�
a <- primal(sv)
aa <- a>0   #a>0�̂ݎ��o��
alpha <- matrix(a[aa], nrow=nrow(xx[aa, ]), ncol=ncol(xx))
(beta <- colSums(alpha * xx[aa, ]))   #�T�|�[�g�x�N�^�[��p����beta���Z�o

#beta����b�����߂�
b <- mean(YXl[aa, 1] - as.matrix(YXl[aa, -1]) %*% as.matrix(beta))

##���ތ���
#�w�K�f�[�^�̏ꍇ
resl <- as.matrix(YXl[, -1]) %*% as.matrix(beta) + b    #���苫�E
resll <- sign(resl)
sum(as.numeric(resll[1:500]==1))/500
sum(as.numeric(resll[501:1000]==-1))/500

#�e�X�g�f�[�^�̏ꍇ
rest <- as.matrix(YXt[, -1]) %*% as.matrix(beta) + b   #���苫�E
restt <- sign(rest)
sum(as.numeric(restt[1:500]==1))/500
sum(as.numeric(restt[501:1000]==-1))/500

#���ʌ��ʂ̃N���X�W�v
testerror1 <- cbind(YXt[, 1], restt)
table(testerror1[, 1], testerror1[, 2])
 
##�֐����g���ꍇ
model <- ksvm(factor(YXl[, 1]) ~ ., data=YXl, kernel="vanilladot", prob.model=FALSE)
model
print(model@alpha)   #alpha�̌W��
print(model@b)   #beta�̌W��
print(model@SVindex)   #�T�|�[�g�x�N�^-

#���ʌ���
pre <- predict(model, as.data.frame(YXt), type="response")   #�\������
pre <- as.numeric(pre)
pre[pre==1] <- -1 
pre[pre==2] <- 1
testerror2 <- cbind(YXt[, 1], pre)
table(testerror2[, 1], testerror2[, 2])

####�J�[�l���@��p�����T�|�[�g�x�N�^�[�}�V��####
##�J�[�l���֐��̒�`
##�O�����s����쐬����
#�������J�[�l��
gram1 <- (3 + as.matrix(Xnew) %*% t(as.matrix(Xnew))) + (3 + as.matrix(Xnew) %*% t(as.matrix(Xnew)))^2
round(gram1[1:15, 1:15], 3)
round(gram1[985:1000, 985:1000], 3)
round(eigen(gram1)$value, 3)   #������l���ǂ����m�F

#�K�E�X�J�[�l��
sigma <- 1/2
kf_gauss <- function(x1, sigma){
  x1 <- as.matrix(x1)
  g1 <- matrix(t(x1), nrow(x1)^2, ncol(x1), byrow=T)
  g2 <- matrix(rep(x1, c(rep(nrow(x1), nrow(x1)*ncol(x1)))), nrow(x1)^2, ncol(x1))
  g3 <- (g2 - g1)^2
  gram <- exp(-sigma*sqrt(matrix(apply(g3, 1, sum), nrow(x1), nrow(x1), byrow=T)))
  return(gram)
}
gram2 <- kf_gauss(x1=Xnew, sigma=sigma)
round(gram2[1:15, 1:15], 3)
round(gram2[985:1000, 985:1000], 3)
round(eigen(gram2)$value, 3)   #������l���ǂ����m�F

#�֐����g��
L <- kernelMatrix(rbfdot(sigma=1/2), Xnew)   #�K�E�X�J�[�l���ŕϊ�

#forloop�ŃO�����s����쐬����ƌv�Z���Ԃ�����
#for(i in 1:n){
#  for(j in 1:n){
#    r <- sum(Xnew[i, ] - Xnew[j, ])
#    gm[i, j] <- r
#  }
#  print(i)
#}
