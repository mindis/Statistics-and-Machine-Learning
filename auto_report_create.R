#' # �������|�[�g�쐬
# /*
library(knitr)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)
# */
#'
#' ## �f�[�^�̔���

# /* 
#set.seed(4543) 
# */

#'
#' ### �����ϐ��ƃp�����[�^�̐ݒ�
#' #### �f�[�^�ƃp�����[�^�̐ݒ�
#+ data
col <- 10   #�p�����[�^��
x <- rnorm(5000, 0, 3)
alpha <- 0.5
beta1 <- rnorm(10, 0, 0.7)
Xm <- matrix(x, nrow=length(x)/col, ncol=col, byrow=T)  
X <- as.data.frame(Xm)
X1 <- cbind(1, Xm)
n <- nrow(X)   #�T���v����

#' ### �I���m���̌v�Z�Ɖ����ϐ��̔���
#' #### �I���m���̌v�Z
p <- exp(X1 %*% as.vector(beta)) / (1 + exp(X1 %*% as.vector(beta)))

# /* 
#p <- plogis(alpha + beta1 * x)   #����ł�OK
# */

#' - �I���m���̕���
round(mean(p), 3)
#' - �I���m���̕��ʓ_
round(quantile(p), 3)

#' - �����������m���m���̕��z
hist(p, breaks=25 , col="#0000ff40", border = "#0000ff", main="�����������m���̕��z", xlab="�m��")


#' ### �I�����ʂ𔭐�������
choice <- rbinom(n, 1, p)
#' - �I�����ʂ̒P���W�v
mean(choice)
table(choice)

# /*
#�f�[�^������
Xy <- cbind(choice, X) 
Xy <- as.data.frame(Xy)
round(head(Xy, 10), 3)
# */

#' ## ���W�X�e�B�b�N��A���f���𐄒�
#' ### �ΐ��ޓx�̐ݒ�
#+ analysis
fr <- function(b, x, y){
  #' - �p�����[�^�̐ݒ�
  alpha <- b[1]
  beta <- b[2:11]
  
  #' - �ޓx���`���č��v����
  Xb <- alpha + as.matrix(X) %*% as.vector(beta) 
  p <- exp(Xb) / (1 + exp(Xb))
  LLS <- y*log(p) + (1-y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

#' ### �ΐ��ޓx���ő剻����
#' #### ���j���[�g���@�Ő���
b0 <- c(rep(0, ncol(X)+1))   #�����p�����[�^�̐ݒ�
res <- optim(b0, fr, gr=NULL, x=X, y=choice, method="BFGS", hessian=TRUE, control=list(fnscale=-1))

#' #### ���茋�ʂƓK���x
#' - ���肳�ꂽ�p�����[�^
round(b <- res$par, 3)   
#' - �^�̃p�����[�^
round(beta, 3)   
#' - t�l
round(tval <- b/sqrt(-diag(solve(res$hessian))), 3)   
#' - AIC
round(AIC <- -2*res$value + 2*length(res$par), 3)   
#' -BIC
round(BIC <- -2*res$value + log(n)*length(b), 3)   #BIC

#' #### �K���x
prob <- round(exp(as.matrix(X1) %*% as.vector(b)) / (1 + exp(as.matrix(X1) %*% as.vector(b))), 3)   #���肳�ꂽ�m��
#' - �^�̊m���Ƃ̌덷
cbind(round(p, 3), prob , round(p-prob, 3))[1:10, ]   

# /*
rbind(beta, b, beta-b)   #�p�����[�^�̌덷
# */

#' ## �֐����g���Ȃ�
#+ func
res2 <- glm(choice ~ V1 + V2 + V3 + V4 +V5 + V6 + V7 + V8 + V9 + V10, data = Xy, family=binomial(link=logit))
#' #### �֐�glm�̐��茋��
summary(res2)

# /*
#��������o�͂���܂���
glmb <- coef(res2)
error <- rbind(b, glmb, beta, b-beta, glmb-beta)   #�p�����[�^�̌덷
rownames(error) <- c("beta", "glmb", "btrue", "beta-btrue", "glmb-btrue")
#�����܂ŏo�͂���܂���
# */

#' ## ���茋�ʂ̕\
#+ results="asis", echo=FALSE
kable(round(error, 3))

# /*
spin("D:/Statistics/statistics_master/Simulation/auto_report_create.R")
# */

