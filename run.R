rm(list = ls())
# library(odeintr)
library(rstudioapi)
library(Rcpp)
library(deSolve)
# library(sundialr)

setwd(dirname(getActiveDocumentContext()$path))
# set.seed(15)

source('FW_functions.R')
source("parameters.R")
source("create_model.R")
sourceCpp("define_parameter_class.cpp") # O2 optimisation by default (to check)

# deSolve ask the ODE functions to return a list
lsoda.wrapper = function(t, y, parms){
  return(list(m$ODE(y, t)))
}

# define a community body mass structure
# repeat the creation of the L matrix while I have isolated basal species
# I thus need to redraw community body masses
isolated = TRUE
comp = -1
while (isolated){
  # not explicit in the paper, but I think BM are drawn with a uniform dist
  # and then raised to 10 power of 
  BM = c(sort(runif(nb_b, 1, 6)), sort(runif(nb_a, 2, 9)))
  BM = 10^BM
  
  L = create.L2(BM, Ropt, gamma)
  # Lbis = create.L(BM, Ropt, gamma)
  # herbivores are ... herbivores
  L[,1:nb_b] = 0.0
  # Lbis[,1:nb_b] = 0.0
  # check the presence of isolated species or consumers without prey:
  isolated = any(colSums(L) + rowSums(L) == 0) | any(colSums(L[,(nb_b+1):nb_s])==0)
  comp = comp + 1
  # colSums(L[,(nb_b+1):nb_s])
}
cat('number of rejected webs = ', comp, '\n')


# create the binary FW structure
FW = L
FW[L>0] = 1
# create the model object
m = create_model(nb_s, nb_b, nb_n, BM, FW, L)
time_vec = seq(0,1000000,10000)
parms = c(0,1) # useless, but mandatory?
  # run the model
df =  tryCatch(lsoda(m$bioms, time_vec, lsoda.wrapper, parms), 
               warning= function(cond){return(rbind(rep(NA, nb_s + nb_n), rep(NA, nb_s + nb_n)))}
)

length(m$dB) == length(m$bioms)
col = c(rep('black', 2), rep('green', nb_b), rep('red', nb_a))
ymax = max(df[,2:ncol(df)])

# df= df[1:20, ]
plot(df[,2]~ df[,1], type = 'n', ylim = c(0,ymax))
for (i in 2:ncol(df)){
  lines(df[,i]~ df[,1], col = col[i-1])
}
# selecting the last row (final biomasses)
nrow(df)
df[nrow(df), ]


sum(df < 0.0000001) >1
# 
tail(df)
# library(rbenchmark)
# benchmark(odeintr::integrate_sys(sys = m$ODE, init = m$bioms, duration = 1500), replications = 5)

