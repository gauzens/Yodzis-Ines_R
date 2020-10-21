# look to the effect of the number of herbivore species on food web persistence
rm(list = ls())
library(odeintr)
library(rstudioapi)
library(Rcpp)
setwd(dirname(getActiveDocumentContext()$path))
set.seed(123)

source('FW_functions.R')
source("parameters.R")
sourceCpp("define_parameter_class.cpp")
source("create_model.R")

# redefine the number of species
nb_s = 15

# number of replicates:
nreps = 4
# create an object to save values
saves = list()

return_persistence = function(nb_b){
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
  df = odeintr::integrate_sys(sys = m$ODE, init = m$bioms, duration = 150000)
  # save the simulartion results:
  saves[[nb_b]] = df
  # save(saves, file = "one.file.Rdata")
  # calculate food web persistence:
  # remove first column from df because it contains time
  # and also the 2nd and 3rd because they contains nutrient biomass 
  to.remove = c(1,2,3)
  persistence = sum(df[nrow(df), -to.remove] > 0.0000001) / nb_s
  
  # destroy the model object (don't know if this is needed, but at least safer)
  rm(m)
  results = c(nb_b, persistence)
  names(results) = c("nb_b", "persistence")
  return(results)
}

nreps = 4

n_basals = rep(seq(1:8), nreps)

results = sapply(n_basals, return_persistence)
# how long you need to run 1 replicate on your computer? 
nreps = 1
n_basals = rep(seq(1:8), nreps)
library(rbenchmark)
benchmark(sapply(n_basals, return_persistence), replications = 1)


# multiprocessing seems to not be fine with rcpp classes
# library("future")
# library(future.apply)
# plan(multiprocess)
# future_sapply(n_basals, return_persistence, parameters)

results = t(results)
plot(results[,2] ~ results[,1])
boxplot(results[,2] ~ as.factor(results[,1]))

library(ggplot2)
results = as.data.frame(results)
ggplot(results, aes(x = nb_b, y = persistence))+
  geom_point()+
  geom_smooth()
  # geom_boxplot(aes(group=nb_b))



