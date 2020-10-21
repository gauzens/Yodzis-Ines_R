# look to the effect of the number of herbivore species on food web persistence
rm(list = ls())
library(odeintr)
library(rstudioapi)
library(Rcpp)
library(deSolve)
setwd(dirname(getActiveDocumentContext()$path))
set.seed(123)

source('FW_functions.R')
source("parameters.R")
sourceCpp("define_parameter_class.cpp")
source('create_model_with_temp.R')

lsoda.wrapper = function(t, y, parms){
  return(list(m$ODE(y, t)))
}

# redifine the number of species
nb_s = 50
nb_b = 15
nb_a = nb_s - nb_b
# create an object to save values
saves = list()
time_vec = seq(0,1000000,10000)
parms = c(0,1) # useless, but mandatory?

persistence = c()
sigmas = seq(1:20) # temperature
temp = 293.15
for (sigma in sigmas){
  cat('starting calculations at ', temp, ' degrees\n')

  # first generate a food web 
  isolated = TRUE
  comp = -1
  while (isolated){
    BM = c(sort(runif(nb_b, 1, 6)), sort(runif(nb_a, 2, 9)))
    BM = 10^BM
    
    L = create.L2(BM, Ropt, gamma)
    
    # herbivores are ... herbivores
    L[,1:nb_b] = 0.0
    
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
  m = create_model_temp(nb_s, nb_b, nb_n, BM, FW, L, 293.15, sigma)
  # m$b = 
  # run the model
  df =  tryCatch(lsoda(m$bioms, time_vec, lsoda.wrapper, parms), 
                 warning= function(cond){return(rbind(rep(NA, nb_s + nb_n), rep(NA, nb_s + nb_n)))}
  )
  # save the simulartion results:
  saves[[nb_b]] = df
  # save(saves, file = "one.file.Rdata")
  # calculate food web persistence:
  persistence = c(persistence, sum(df[nrow(df), ] > 0.0000001) / nb_s)
  
  # destroy the model object (don't know if this is needed, but at least safer)
  rm(m)
}


plot(persistence ~ sigmas)






