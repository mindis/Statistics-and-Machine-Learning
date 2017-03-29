#####���W�X�e�B�b�N��A���f��#####
####�f�[�^�̔���####
#set.seed(4543)
#�����ϐ��ƃp�����[�^�̐ݒ�
col <- 10   #�p�����[�^��
x <- rnorm(5000, 0, 3)
alpha <- 0.5
beta1 <- rnorm(10, 0, 0.7)
(beta <- c(alpha, beta1))
Xm <- matrix(x, nrow=length(x)/col, ncol=col, byrow=T)  
X <- as.data.frame(Xm)
X1 <- cbind(1, Xm)
n <- nrow(X)   #�T���v����

#�I���m���̌v�Z
p <- exp(X1 %*% as.vector(beta)) / (1 + exp(X1 %*% as.vector(beta)))
#p <- plogis(alpha + beta1 * x)   #����ł�OK
mean(p)

#�I�����ʂ𔭐�������
choice <- rbinom(n, 1, p)
mean(choice)

#�f�[�^������
Xy <- cbind(choice, X) 
Xy <- as.data.frame(Xy)
round(head(Xy, 20), 4)

####���W�X�e�B�b�N��A���f���𐄒�####
##�ΐ��ޓx�̐ݒ�
fr <- function(b, x, y){
  #�p�����[�^�̐ݒ�
  alpha <- b[1]
  beta <- b[2:11]
  
  #�ޓx���`���č��v����
  Xb <- alpha + as.matrix(X) %*% as.vector(beta) 
  p <- exp(Xb) / (1 + exp(Xb))
  LLS <- y*log(p) + (1-y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

##�ΐ��ޓx���ő剻����
b0 <- c(rep(0, 11))   #�����p�����[�^�̐ݒ�
res <- optim(b0, fr, gr=NULL, x=X, y=choice, method="BFGS", hessian=TRUE, control=list(fnscale=-1))

#����
(b <- res$par)   #���肳�ꂽ�p�����[�^
beta   #�^�̃p�����[�^
(tval <- b/sqrt(-diag(solve(res$hessian))))   #t�l
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(n)*length(b))   #BIC

#�K���x
prob <- round(exp(as.matrix(X1) %*% as.vector(b)) / (1 + exp(as.matrix(X1) %*% as.vector(b))), 3)   #���肳�ꂽ�m��
cbind(round(p, 3), prob , round(p-prob, 3))   #�^�̊m���Ƃ̌덷
rbind(beta, b, beta-b)   #�p�����[�^�̌덷

#�֐����g���Ȃ�
res2 <- glm(choice ~ V1 + V2 + V3 + V4 +V5 + V6 + V7 + V8 + V9 + V10, data = Xy, family=binomial(link=logit))
summary(res2)
glmb <- coef(res2)
round(rbind(beta, b, glmb, beta-b, glmb-b), 3)   #�p�����[�^�̌덷