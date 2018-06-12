#####Modeling the Stock and Threshold Effect of Personal Selling#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(actuar)
library(extraDistr)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

####ÉfÅ[É^ÇÃî≠ê∂####
hh <- 2000  
pt <- rpois(hh, rgamma(hh, 20, 0.5))

      