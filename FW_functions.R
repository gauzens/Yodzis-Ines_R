


# create the L.matrix
create.L = function(BM, Ropt, gamma){
  fill = function(x){
    # return(   (  (BM/(BM[prey]*Ropt)) * exp(1-BM/(BM[prey]*Ropt))  )^gamma   )
    # cat(gamma)
    return(   (  (x/(BM*Ropt)) * exp(1-x/(BM*Ropt))  )^gamma   )
  }
  nb_s = length(BM)
  # matrix.L = sapply(seq_len(nb_s), fill)
  matrix.L = t(sapply(BM, fill))
  matrix.L[matrix.L < 0.01] = 0

  return(matrix.L)
}

create.L2 = function(BM, Ropt, gamma){
  s = length(BM)
  L = matrix(NA, nrow = s, ncol = s)
  for (prey in 1:s){
    for (pred in 1:s){
      L[prey, pred] = (  (BM[pred]/(BM[prey]*Ropt)) * exp(1-BM[pred]/(BM[prey]*Ropt))  )^gamma
    }
  }
  L[L < 0.01] = 0
  return(L)
}


# create the matrix of b values
create.b = function(BM, b0, bprey, bpred, L){
  nb_s = length(BM)
  M = matrix(1, nrow = nb_s, ncol = nb_s)
  prey = (M*BM)^bprey
  prey[1:nb_b] = b0*20
  return(b0 * prey * t(M*BM)^bpred * L)
  # alternatively: 
  # fill = function(prey){
  #   return(b0*(BM[prey]^bprey)*(BM^bpred)*L[prey,])
  # }
  # matrix.L = sapply(seq_len(nb_s), fill)
}

# create_model = function(nb_s, nb_b, nb_n, BM, FW, L){
#   create.matrix.param = function(BM, b0, bprey, bpred){
#     nb_s = length(BM)
#     M = matrix(1, nrow = nb_s, ncol = nb_s)
#     return(b0 * (M*BM)^bprey * t(M*BM)^bpred)
#   }
#   model = new(parameters,nb_s,nb_b,nb_n)
#   model$BM = BM
#   # model$FW = FW
#   model$K1 = runif(nb_b, nut_up_min, nut_up_max)
#   model$K2 = runif(nb_b, nut_up_min, nut_up_max)
#   model$D = D
#   model$S = rnorm(nb_n, mu_nut, sd_nut)
#   model$r = BM[1:nb_b]^-0.25
#   model$X = c(rep(x_P, nb_b), rep(x_A, nb_s - nb_b)) * BM^-0.25
#   model$e = c(rep(e_P, nb_b), rep(e_A, nb_s - nb_b))
#   model$w = colSums(FW)
#   model$b = create.b(BM,b0, bprey, bpred, L)
#   model$c =  rnorm(1, mu_c, sd_c)
#   model$h = create.matrix.param(BM,h0, hprey, hpred)
#   model$q = rnorm(1, 0.5, 0.2)
#   model$v1 = v1
#   model$v2 = v2
#   model$r = BM[1:nb_b]^-0.25
#   model$F = matrix(0.0, nrow = nb_s, ncol = nb_s - nb_b)
#     # initial biomasses: 
#   bioms = c(0.8*model$S, runif(nb_s, 30,40))
#   model$bioms = bioms
#   return(model)
# }



# calculate species' trophic levels
TL <- function(FW){
  FW <- t(FW)
  nn <- rowSums(FW); nn[nn==0] <- 1
  ww <- diag(1/nn)
  L1 <- ww %*% FW
  L2 <- L1 - diag(rep(1,length(nn)))
  b <- -1*rep(1,length(nn))
  Tro.lev <- solve(L2) %*% b
  
  return(Tro.lev)
}


Omnivory = function(TL, FW){
  nb_s = length(TL)
  sapply(seq_len(nb_s), function(x) sd(TL[FW[,x]!=0]))
}

F.rate = function(bioms, BM, b, c, q, h, w){
  fill = function(pred){
    res = (w[pred]*b[,pred]*bioms^(1+q)) / (1 + c*bioms[pred] + w[pred]*h[pred]*sum(FW[,pred]*b[,pred]*bioms^(1+q)))
  }
}

  