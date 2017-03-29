#####�o�X���f��#####
library(reshape2)
library(plyr)
####�f�[�^�̔���####
#set.seed(2134)
(m <- runif(1, 60, 75))   #���y���̏���l
(p <- runif(1, 0.0005, 0.005))   #�v�V�ҌW��
(q <- runif(1, 0.001, 0.025))   #�͕�ҌW��
t <- 1:365   #����

F <- m*(1-exp(-(p+q)*t)) / (1+q/p*exp(-(p+q)*t))   #�^�̊֐�
Fe <- m*(1-exp(-(p+q)*t)) / (1+q/p*exp(-(p+q)*t)) + rnorm(365, 0, 4)   #�ϑ��l

#���ʂ��v���b�g
plot(t, F, type="l", lwd=2, ylim=c(0, 85), xlab="�o�ߓ���", ylab="���y��", main="�o�X���f���ɂ�镁�y����")
lines(t, Fe, lty=2, lwd=1)

####�o�X���f���̃p�����[�^����####
##����`�ŏ����@�Ő���
M2 <- nls(Fe ~ m*(1-exp(-(p+q)*t)) / (1+q/p*exp(-(p+q)*t)), 
          start=list(p=0.001, q=0.01, m=max(Fe)), m=max(Fe), trace=TRUE)

para <- M2$m$getPars()

#���肳�ꂽ�p�����[�^�ɂ��\��
(p.hat <- para[1])
(q.hat <- para[2])
(m.hat <- para[3])

Mhat <- m.hat*(1-exp(-(p.hat+q.hat)*t)) / (1+q.hat/p.hat*exp(-(p.hat+q.hat)*t))   #�\���l

#���ʂ��v���b�g
plot(t, Mhat, type="l", lwd=2, ylim=c(0, 85), xlab="�o�ߓ���", ylab="���y��", main="�o�X���f���ɂ�镁�y���\��")
lines(t, Fe, lty=2, lwd=1)   #�����l
legend("topleft", legend=c("Predicted", "observed"), lty=1:2)

summary(M2)
