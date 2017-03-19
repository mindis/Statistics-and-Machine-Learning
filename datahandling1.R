####データハンドリング練習用####
options(digits=10)
(cum <- cumsum(runif(100, 0, 0.05)))
(cumprob <- subset(cum, cum < 1))
p <- c(cumprob[1])
for(i in 2:length(cumprob)){
  pp <- cumprob[i] - cumprob[i-1]
  p <- c(p, pp)
}

p_last <- 1 - cumprob[length(cumprob)]
(prob <- c(p, p_last))

choice <- t(rmultinom(10000000, 1, prob)) %*% c(1:length(prob))
dummy <- matrix(0, nrow(choice), length(table(choice)))
for(i in 1:nrow(choice)){
  iid <- choice[i]
  dummy[i, choice[i]] <- 1
}
colMeans(dummy)
colSums(dummy)
vector()

vec <- vector()
dosu <- rmultinom(1, n[1], p)
sa <- max(dosu)/sum(dosu) - min(dosu)/sum(dosu) 
vec <- c(vec, sa) 
vec

a1 <- runif(10)
p <- rep(0.1, 10)
n <- c(10, 50, 100, 500, 1000, 5000, 10000, 1000000, 10000000, 100000000)
diff1 <- matrix(0, 10000, length(n))
for(i in 1:nrow(diff1)){
  vec <- vector()
  for(j in 1:length(n)){
    dosu <- rmultinom(1, n[j], p)
    sa <- max(dosu/sum(dosu)) - min(dosu/sum(dosu)) 
    vec <- c(vec, sa) 
  }
  diff1[i, ] <- vec
}
round(diff1, 5)
diff1

a <- round(apply(diff1, 2, max), 6)
b <- round(apply(diff1, 2, min), 6)
c <- a - b
rbind(a, b, c)

diff1
system.time(pp <- rmultinom(1000, 100, p))
t(pp)

subset(1:length(table(choice)), choice[1] == 1)


load("cheese.rda")
Q <- log(cheese$VOLUME)
P <- log(cheese$PRICE)
Prom <- cheese$DISP
Stores <- unique(cheese[, 1])
cheese
(freq <- data.frame(table(cheese[, 1])))
(flag <- as.vector(unique(cheese[, 1])[freq[, 2] < 65]))

data <- data.frame()
for(i in 1:length(flag)){
  d <- subset(cheese, (cheese[, 1] == flag[i]) == 1)
  data <- rbind(data, d)
}
(d1 <- as.matrix(data))
(d1 <- as.data.frame(d1)) 
unique(d1[, 1])
unique(as.vector(d1[, 1]))
head(d1, 100)
d1[2, ]

L <- matrix(0, nrow=nrow(cheese), ncol=88)
for(i in 1:nrow(cheese)){
  ind <- subset(1:88, (cheese[i, 1]==Stores)==1)
  L[i, ind] <- 1
}