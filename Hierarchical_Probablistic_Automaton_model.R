#####Hierarchical Probablistic Automaton model#####
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
k1 <- 6   #��ʊK�w
k2 <- 6   #���ԊK�w
k3 <- 10   #���ʊK�w
d <- 3000   #������
v1 <- 700   #��ʌ�
v2 <- 200   #�Ɨ���
v3 <- 200   #�@�\��
v4 <- 25   #��؂�L��
v <- v1 + v2 + v3 + v4   #����b��
w <- rpois(d, rgamma(d, 65, 0.4))   #��b��
f <- sum(w)   #����b��

#ID�̐ݒ�
d_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}

##�p�����[�^�̎��O���z�̐ݒ�
#�������z�̎��O���z
pi01 <- rep(5.0, k1-1)

#�}���R�t���ڍs��̎��O���z
alpha01 <- c(rep(50.0, k1-1), 30)   #�ŏ�ʊK�w�̎��O���z















