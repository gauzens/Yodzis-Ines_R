# look to the effect of the number of herbivore species on food web persistence
rm(list = ls())
library(odeintr)
library(rstudioapi)
library(Rcpp)
setwd(dirname(getActiveDocumentContext()$path))

set.seed(2)

source('FW_functions.R')
source("parameters.R")
source("create_model.R")
sourceCpp("define_parameter_class.cpp")

lsoda.wrapper = function(t, y, parms){
  return(list(m$ODE(y, t)))
}


# redifine the number of species
nb_s = 15

# create an object to save values
saves = list()

persistence = c()
nb.basal.sp = seq(1:10)
time_vec = seq(0,1000000,10000)
parms = c(0,1) # useless, but mandatory?

for (nb_b in nb.basal.sp){
  cat('starting calculations for ', nb_b, ' basal species\n')
  # adapt the number of ANIMAL species
  nb_a = nb_s - nb_b
  # first generate a food web with nb basal species
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
  m = create_model(nb_s, nb_b, nb_n, BM, FW, L)
  # run the model
  df =  tryCatch(lsoda(m$bioms, time_vec, lsoda.wrapper, parms), 
                 warning= function(cond){return(rbind(rep(NA, nb_s + nb_n), rep(NA, nb_s + nb_n)))}
  )
  # save the simulartion results:
  saves[[nb_b]] = df
  # save(saves, file = "one.file.Rdata")
  # calculate food web persistence:
  # remove first column from df because it contains time
  # and also the 2nd and 3rd because they contains nutrient biomass 
  to.remove = c(1,2,3)
  not.extinct = df[nrow(df), -to.remove] > 0.0000001
  web.persistence = sum(not.extinct) / nb_s
  persistence = c(persistence, web.persistence)
  
  # destroy the model object (don't know if this is needed, but at least safer)
  rm(m)
}


plot(persistence ~ nb.basal.sp)






