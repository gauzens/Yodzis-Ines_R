create_model_temp = function(nb_s, nb_b, nb_n, BM, FW, L, temp, sigma){
  create_matrix_parameter = function(BM, b0, bprey, bpred){
    nb_s = length(BM)
    M = matrix(1, nrow = nb_s, ncol = nb_s)
    return(b0 * (M*BM)^bprey * t(M*BM)^bpred)
  }
  #initial biomasses: 
  
  model = new(parameters,nb_s,nb_b,nb_n)
  model$BM = BM
  model$K1 = runif(nb_b, nut_up_min, nut_up_max)
  model$K2 = runif(nb_b, nut_up_min, nut_up_max)
  model$D = D
  model$S = rnorm(nb_n, mu_nut, sd_nut)
  model$r = BM[1:nb_b]^-0.25
  model$X = c(rep(x_P, nb_b), rep(x_A, nb_s - nb_b)) * BM^-0.25
  model$e = c(rep(e_P, nb_b), rep(e_A, nb_s - nb_b))
  model$w = colSums(FW)
  model$b = create_matrix_parameter(BM,b0, bprey, bpred) * L * exp((temp - 293.45)^2/(2*sigma))
  # model$FW = FW
  model$c =  rnorm(1, mu_c, sd_c)
  model$h = create_matrix_parameter(BM,h0, hprey, hpred)
  model$q = rnorm(1, 0.5, 0.2)
  model$v1 = v1
  model$v2 = v2
  model$r = BM[1:nb_b]^-0.25
  model$F = matrix(0.0, nrow = nb_s, ncol = nb_s - nb_b)
  # initial biomasses:
  bioms = c(0.8*model$S, runif(nb_s, 30,40))
  model$bioms = bioms
  return(model)
  
  
}
