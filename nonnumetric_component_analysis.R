#####��v�ʎ听������#####
library(homals)
library(MASS)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)


####�f�[�^�̓ǂݍ���####
data(galo)
Data <- galo

####�f�[�^�N�����W���O####
##�J�e�S���J���ϐ����_�~�[�ϐ��ɕϊ�
#advice���_�~�[�s��ɕϊ�
advice_table <- table(1:nrow(Data), Data$advice)
advice_name <- colnames(advice_table)
advice <- matrix(as.numeric(advice_table), nrow=nrow(Data), ncol=ncol(advice_table))
colnames(advice) <- advice_name

#SES���_�~�[�s��ɕϊ�
SES_table <- table(1:nrow(Data), Data$SES)
SES_name <- colnames(SES_table)
SES <- matrix(as.numeric(SES_table), nrow=nrow(Data), ncol=ncol(SES_table))
colnames(SES) <- SES_name

#���ʂ�ϊ�
sex <- ifelse(Data$gender=="M", 1, 0)   

##�f�[�^�̌���
X = cbind(sex, IQ=Data$IQ, advice, SES, school=Data$School)


####���ݍŏ����@�Ŕ�v�ʎ听�����͂𐄒�####
##�p�����[�^�̏����l��ݒ�
