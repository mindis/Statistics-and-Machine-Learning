#####�m���p�����g���b�N��A#####
library(MASS)
library(fda)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(324)
##�f�[�^�̐ݒ�
x <- c(1:(52*3))   #�����P�ʂ�3�N

##�o�̓f�[�^�𔭐�
tb <- 2.0   #�����l
trend <- numeric()
s <- seq(0.85, 0.2, length=length(x))
for(i in 1:length(x)){
  r <- rnorm(5, tb, 0.02)
  sort <- sort(r)
  bi <- rbinom(1, 1, s[i])
  bb <- ifelse(bi == 1, sort[4], sort[2])
  tb <- bb
  trend <- c(trend, bb)
}
plot(trend, lwd=1, xlab="half_month", ylab="p")

####3���X�v���C�����####
##�ߓ_�����߂�
t <- seq(x[1], x[length(x)], length=16)   #�ߓ_

##�ߓ_����Ǐ��I�ȑ����������
b <- create.bspline.basis(rangeval=c(t[1], t[length(t)]), breaks=t)
plot(b)

#�����������s
sm <- smooth.basis(x, trend, b); sm   #�����������s

##���ʂ��v���b�g
plot(trend, lwd=1, xlab="half_month", ylab="p")
lines(sm, lwd=2, col="red")

