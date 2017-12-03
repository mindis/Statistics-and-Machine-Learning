####ログサム変数実験####
library(matrixStats)
library(gtools)

seg <- 5
par <- 10000
p <- rdirichlet(seg, rep(0.75, par))
x <- as.numeric(t(rmultinom(1, 200, p[1, ])))

#通常の尤度の積
pi01 <- dmultinom(x, 200, p[1, ])
pi02 <- dmultinom(x, 200, p[2, ])
pi03 <- dmultinom(x, 200, p[3, ])
pi04 <- dmultinom(x, 200, p[4, ])
pi05 <- dmultinom(x, 200, p[5, ])
pi00 <- c(pi01, pi02, pi03, pi04, pi05)

#潜在変数zの割当確率
1/seg*pi00 / sum(1/seg*pi00)


#logsumexpの尤度の積
pi11 <- dmultinom(x, 200, p[1, ], log=T)
pi12 <- dmultinom(x, 200, p[2, ], log=T)
pi13 <- dmultinom(x, 200, p[3, ], log=T)
pi14 <- dmultinom(x, 200, p[4, ], log=T)
pi15 <- dmultinom(x, 200, p[5, ], log=T)
pi11 <- c(pi11, pi12, pi13, pi14, pi15)

expl <- 1/seg * exp(pi11 - max(pi11))
exp(log(expl) - (log(sum(exp(log(expl) - max(log(expl))))) + max(log(expl))))


(pi0) - (log(sum(exp(pi0 - max(pi0)))) + max(pi0))

    