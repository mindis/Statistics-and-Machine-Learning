####���W�X�e�B�b�N��A���f���ɂ��y�[�W�����N���####
library(MASS)
library(mlogit)
library(nnet)
library(reshape2)
library(dplyr)
library(caret)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
h <- 50   #�T���v����
n <- 300   #�������[�h��

#ID�̐ݒ�
id <- rep(1:h, rep(n, h))
w <- rep(1:n, h)
ID <- data.frame(no=1:length(w), id, w)

##�����N�f�[�^�̔���
#�����Ɉ��������邩�ǂ����̃f�[�^
#�p�����[�^�̐ݒ�
#�o���L���̃f�[�^�̃p�����[�^
theta0 <- rnorm(h, 0, 1.0)   #�T���v���̃p�����[�^
theta0[h] <- 0   #��ϐ���0�ɒu������
theta0.b <- -1.0   #�ؕ�
rho <- 0.4   #���O�T���ϐ��̃p�����[�^

#�����N�f�[�^�̃p�����[�^
theta1 <- theta0 + runif(h, -1.25, 1.25)   #�����N�f�[�^�̃p�����[�^
theta1[h] <- 0

#�o���L���̃��W�b�g�Ɗm�����v�Z
logit0 <- theta0.b + theta0 + rho*log(exp(theta1)) 
p0 <- exp(logit0)/(1+exp(logit0))
Pr0 <- rep(p0, rep(n, h))

#�񍀕��z�ɂ�艞���ϐ��𔭐�
y0 <- c()
for(i in 1:h){
  y0 <- c(y0, rbinom(n, 1, p0[i]))
}

##�������ʂ𔭐�������
y1 <- rep(0, nrow(ID)) 

#�I���W���̎擾
for(i in 1:n){
  print(i)
  in_w <- subset(ID$id, ID$w==i & y0==1)   #�I���W�����擾
  c_cnt <- length(in_w)
  
  #20�Ԗڂ̏��ʂ܂ŋL�^����
  if(c_cnt < 20){
    c_cnt <- c_cnt 
    } else {
      c_cnt <- 20
    }
  
  if(c_cnt < 2) next   #�I�𐔂�2�ȉ��Ȃ玟��

#�����N�f�[�^�̔���
  for(j in 1:(c_cnt-1)){
    #�m���̌v�Z�Ɖ����ϐ��̔���
    Pr <- exp(theta1[in_w])/sum(exp(theta1[in_w]))
    choise <- t(rmultinom(1, 1, Pr))
    
    #�����ϐ��������đI���W�����X�V
    index.z <- subset(1:length(choise), choise==1)
    y1[ID$w==i & ID$id==in_w[index.z]] <- j
    
    in_w <- in_w[-index.z]   #�I���W�����X�V
  }
}
YX <- data.frame(ID, y0, y1)   #�f�[�^�̌���
 
#id���Ƃɏ��ʂ̏W�v
agg <- YX %>%  
          dplyr::filter(y0==1, y1!=0) %>% 
          dplyr::group_by(id) %>%
          dplyr::summarize(mean=mean(y1), max=max(y1), min=min(y1))

round(data.frame(agg), 1)
xtabs(~ ID$id[y0==1])
xtabs(~ ID$id[y1>0])


####�y�[�W�����N���W�b�g���f�����Ŗސ���####
##�����L���O�f�[�^�̐ݒ�
X.rank <- matrix(0, nrow=sum(y1 > 0), ncol=h)
Y.rank <- matrix(0, nrow=sum(y1 > 0), ncol=h)
len <- c()
rank_sum <- 0
rank <- c()
w.id <- c()

for(i in 1:n){
  print(i)
  rank_len <- length(YX[ID$w==i & y1!=0, "y1"])
  rank_sum <- rank_sum + rank_len
  len <- c(len, rank_len)
  r1 <- (rank_sum - rank_len) + 1
  
  #�����N�W���̐ݒ�
  x.id <- YX[ID$w==i & y0>0, "id"]
  x.rank <- rep(0, h)
  x.rank[x.id] <- 1
  X.rank[r1, ] <- x.rank
  
  #�I�������N�̐ݒ�
  y.id <- YX[ID$w==i & y1==1, "id"]
  y.rank <- rep(0, h)
  y.rank[y.id] <- 1
  Y.rank[r1, ] <- y.rank
  
  #id�̐ݒ�
  w.id <- c(w.id, rep(i, rank_len))   #word_id�̍X�V
  rank <- c(rank, 1)   #�����N�̍X�V
  
  for(j in 2:rank_len){
    r2 <- (rank_sum - rank_len) + j
    rank <- c(rank, j)   #�����N�̍X�V
   
    #�I�����ꂽ�f�[�^�͑I���W�������菜��
    #�����N�W���̍X�V
    x.rank[subset(1:length(y.rank), y.rank==1)] <- 0
    X.rank[r2, ] <- x.rank 
    
    #�I�������N�̍X�V
    y.id <- YX[ID$w==i & y1==j, "id"]
    y.rank <- rep(0, h)
    y.rank[y.id] <- 1
    Y.rank[r2, ] <- y.rank
  }
}


#�v�񓝌v��
len   #�����N�f�[�^��
colSums(Y.rank)   #�����N�I���̏o����
colSums(X.rank)   #�����N�W���̏o����


##�o���L���̃f�U�C���s��̐ݒ�
X.page <- kronecker(diag(h), rep(1, n))
X.page <- X.page[, -h] 


##�y�[�W�����N���W�X�e�B�b�N��A���f���̑ΐ��ޓx
loglike <- function(x, Y.rank, Y.page, X.rank, X.page, r, h){
  #�p�����[�^�̐ݒ�
  theta.rank <- x[r[1]:r[2]]
  b.page0 <- x[r[3]]
  theta.page <- x[r[4]:r[5]]
  rho <- x[r[6]]
  
  ##�ΐ��ޓx�̐ݒ�
  #�����N���W�b�g���f���̑ΐ��ޓx
  #�����N�f�[�^�̌��p�֐�
  logit.rank <- X.rank * cbind(matrix(theta.rank, nrow=nrow(X.rank), ncol=h-1, byrow=T), 0)
  U.rank <- exp(logit.rank) + (X.rank-1)�@
  
  #�ΐ��ޓx�̘a���v�Z
  LLs.rank<- rowSums(Y.rank * logit.rank) - log(rowSums(U.rank))
  LL.rank <- sum(LLs.rank)
  
  #�o���L���̓�l���W�X�e�B�b�N��A���f��
  #�o���L���̌��p�֐�
  logsum <- log(exp(theta.rank))   #���O�T���ϐ��̐ݒ�
  logsum.v <- X.page %*% logsum
  
  logit.page <- b.page0 + X.page %*% theta.page + rho * logsum.v   #���p�֐�
  Pr.page <- exp(logit.page) / (1 + exp(logit.page))   #�m���̌v�Z
  
  #�ΐ��ޓx�̘a���v�Z
  LLs.page <- Y.page*log(Pr.page) + (1-Y.page)*log(1-Pr.page)
  LL.page <- sum(LLs.page)
  
  #�ΐ��ޓx�����v
  LL <- sum(LL.rank + LL.page)
  return(LL)
}

##���j���[�g���@�Ńy�[�W�����N���W�b�g���f���𐄒�
#�p�����[�^�̐ݒ�
r.cum <- c(1, ncol(X.rank)-2, 1, 1, ncol(X.page)-1, 1)
r.par <- cumsum(r.cum)

#���j���[�g���@�őΐ��ޓx���ő剻
#�����l�ݒ�ŃG���[���o���ꍇ�͏����l��ݒ肵�Ȃ����悤�ݒ�
for(i in 1:1000){
  #�����p�����[�^�̐ݒ�
  x <- c(rnorm(h-1, 0, 1), runif(1, -2.0, 0), rnorm(h-1, 0, 1), runif(1, 0.2, 0.6))

  #���j���[�g���@�ōő剻
  res <- try(optim(x, loglike, gr=NULL, Y.rank=Y.rank, Y.page=y0, X.rank=X.rank, X.page=X.page, r=r.par, h=h,
                   method="BFGS", hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #�G���[����
}

####���茋�ʂƗv��####
##���肳�ꂽ�p�����[�^�Ɛ^�̃p�����[�^�̔�r
#���肳�ꂽ�p�����[�^
beta <- res$par

#�y�[�W�����N�̉�A�W��
round(beta1 <- beta[r.par[1]:r.par[2]], 3)
round(theta1, 3)

#�o���L���̉�A�W��
round(beta0 <- beta[r.par[3]:length(beta)], 3)
round(c(theta0.b, theta0[-h], rho), 3)

##�v�񓝌v��
res$value   #�ő剻���ꂽ�ΐ��ޓx
round(tval <- beta/sqrt(-diag(ginv(res$hessian))), 3)   #t�l
round(AIC <- -2*res$value + 2*length(beta), 3)   #AIC


##�K���x
#���肳�ꂽ�o���L���̊m��
logsum <- log(exp(beta1))   #���O�T���ϐ��̐ݒ�
logsum.v <- X.page %*% logsum
logit.page <- cbind(1, X.page, logsum.v) %*% beta0   #���p�֐�
Pr.page <- exp(logit.page) / (1 + exp(logit.page))   #�m���̌v�Z
(Pr.page0 <- round(data.frame(id=1:h, Prt=unique(Pr0), Pr0=unique(Pr.page)), 3))   #�^�̊m���Ƃ̔�r

#���肳�ꂽ�y�[�W�����N�̊m��
logit.rank <- X.rank * cbind(matrix(beta1, nrow=nrow(X.rank), ncol=h-1, byrow=T), 0)
U.rank <- exp(logit.rank) + (X.rank-1)�@
round(Pr.rank <- U.rank/rowSums(U.rank), 2)
(Pr.rank1 <- round(data.frame(w.id, rank, Y=Y.rank %*% 1:h, Pr=Pr.rank), 2))

