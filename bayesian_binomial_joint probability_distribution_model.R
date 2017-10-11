####�񍀕��z���f���̃p�����[�^�̓������胂�f��####

####�f�[�^�̔���####
n <- 800
p <- 0.6
x <- rbinom(1, n, p)

####�}���R�t�A�������e�J�����@�̐ݒ�####
##�A���S���Y���̐ݒ�
R <- 10000
keep <- 2

##���O���z�̐ݒ�
a1 <- 1
b1 <- 1
num <- (x+100):1000
p1 <- rep(1/length((x+100):1000), length((x+100):1000))

##�����l�̐ݒ�
N <- x * 2
Pr <- x / N

##�T���v�����O���ʂ̕ۑ��p�z��
Sample <- rep(0, R/keep)
Prob <- rep(0, R/keep)

####MCMC�Ńp�����[�^���T���v�����O#####
for(rp in 1:R){

  ##�T���v����N���T���v�����O
  old_N <- N
  new_N <- sum((t(rmultinom(1, 1, p1)) * num))
  
  lognew <- dbinom(x, new_N, Pr, log=T)
  logold <- dbinom(x, old_N, Pr, log=T)
  
  #MH�T���v�����O
  alpha <- min(1, exp(lognew - logold))
  if(alpha == "NAN") alpha <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V�����Œ����beta���̑�
  if(u < alpha){
    N <- new_N
    
    #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
  } else {
    N <- old_N
  }
  
  ##�x�[�^���z�����������T���v�����O
  Pr <- rbeta(1, a1+x, b1+N-x)

  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep   
    Sample[mkeep] <- N
    Prob[mkeep] <- Pr
    print(rp)
  }
}


####�T���v�����O���ʂ̊m�F�Ɖ���
burnin <- 1000

#MCMC�T���v�����O���ʂ̉���
matplot(1:(R/keep), Sample, type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(1:(R/keep), Prob, type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
plot(Sample[burnin:(R/keep)], Prob[burnin:(R/keep)], xlab="�T���v����", ylab="�����")

#�p�����[�^�̎��㕽��
mean(Prob[burnin:(R/keep)])
mean(Sample[burnin:(R/keep)])

