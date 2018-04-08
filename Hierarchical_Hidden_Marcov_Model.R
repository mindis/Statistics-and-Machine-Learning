#####Hierarchical Hidden Marcov Model#####
library(MASS)
library(bayesm)
library(HMM)
library(matrixStats)
library(Matrix)
library(extraDistr)
library(reshape2)
library(qrmtools)
library(slfm)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)


hh <- 2000
hist(rpois(hh, rgamma(hh, 25, 1.5)))
