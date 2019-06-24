
#'Function to compute A
#'
#' This function make a Dirichlet distribution
#'
#'@param n The number of patients to simulate.
#'@param prop The different proportion of cell types in a vector of size k (number of cell types).
#'@param alpha0 The variance between the patients, 1 is an high variance and 10,000 is a very low variance.
#'
#'@return This function return a matrix of size k*n with the cell type proportions.
#'
#'@export
compute_A = function(n = 100, prop = c(0.3, 0.6, 0.1), alpha0 = 1){
  t(gtools::rdirichlet(n = n, alpha = prop * alpha0))
}

#'Function to create the experimental group
#'
#' This function make a random experimental group.
#'
#'@param n The number of patients
#'@param plates NA by default, if you want a plate effect put a vector of plates names.
#'@param sex False by default, True if you want a sex parameter.
#'@param age False by defaul, True if you want an age parameter.
#'
#'@return This function return a dataframe with different random covariables.
#'
#'@export
create_exp_grp = function(n, plates = NA, sex = F, age = F){
  expgrp = data.frame(organism = rep("Homo sapien", times = n))
  if(sum(!is.na(plates))){
    plates = sample(plates, n, replace = TRUE)
    expgrp$plates = plates
  }
  if(sex == T){
    id_femme = sample(1:n, n/2)
    sexe = c(1:n)  
    sexe[id_femme] = "F"
    sexe[-id_femme] = "M"
    expgrp$sex = sexe
  }
  if (age == T){
    expgrp$age = sample(45:87, n, replace = TRUE)
  }
  return(expgrp)
}


#'Function to compute the female version of the T matrix
#'
#' This function makes a T matrix with probes modified by sex
#'
#'@param sites The names of probes affected by the sex effect.
#'@param coeff The coefficient of variation for probes affected by the sex.
#'@param T_brut The original T matrix.
#'
#'@return This function return a matrix modified for the sex.
#'
compute_TF = function(sites, coeff, T_brut){
  T_F = T_brut
  T_F[sites,] = (T_F[sites,] - coeff[sites])
  T_F[which(T_F < 0)] = 0
  T_F[which(T_F > 1)] = 1
  return(T_F)
}


#' Computation of the D matrix.
#' 
#' This function compute the D matrix with an effect on age and / or plates
#'
#'@param A_mat The A matrix.
#'@param T_mat The T matrix.
#'@param exp_grp The matrix of experimental data for all the patients.
#'@param sites_sex NA by default, set sites affected by sex if you want a sex effect
#'@param coeff_sex NA by default, set coeff of sex if you want a sex effect.
#'@param plate_effect, Na by default, put plate_effect if you want one
#'
#'@return This function return a matrix of size probes*n with possibly an effect of confounding factors
#'
#'@export
compute_D = function(A_mat, T_mat, exp_grp, sites_sex = NA, coeff_sex = NA,
                           plate_effect = NA){
  #Effect of sex
  if (!sum(is.na(sites_sex)) & sum(!is.na(coeff_sex))){
    T_F = compute_TF(sites_sex, coeff_sex, T_mat)
  }
  D_brut = c()
  for (i in 1:ncol(A_mat)){
    if (exp_grp[i, "sex"] == "F"){
      T_cas = T_F
    } else {
      T_cas = T_mat
    }
    D_cas = T_cas %*% A_mat[, i]
    D_brut = cbind(D_brut, D_cas)
  }
  
  if(sum(!is.na(plate_effect))){
    for(p in 1:ncol(A_mat)){
      plaque = exp_grp$plates[p]
      D_brut[p,] = D_brut[p,] * plate_effect[plaque, 2]
    }
    D_brut[which(D_brut > 1)] = 1
  }
  return(D_brut)
}

#' Add noise to one matrix.
#'
#'@param data The matrix to disturb.
#'@param mean The mean of the normal distribution of noises.
#'@param sd The standard deviation of the normal distribution of noises.
#'@param val_min, the minimum value of the matrix
#'@param val_max, the maximal value of the matrix
#'
#'@return This function return the matrix with noise.
#'
#'@export
add_noise = function(data, mean = 0, sd = 0.2, val_min = 0, val_max = 1){
  noise = matrix(stats::rnorm(prod(dim(data)), mean = mean, sd = sd), 
                 nrow = nrow(data))
  datam = data + noise
  datam[datam < val_min] = data[datam < val_min]
  datam[datam > val_max] = data[datam > val_max]
  return(datam)
}

#' Computing the Mean Absolute Error
#'
#' This function compute the mean absolute error between 2 matrices
#'
#' @param M1 First matrix.
#' @param M2 Second matrix
#' 
#' @return This function return the mean absolute error between the 2 matrices.
#' 
MAE <- function(M1, M2) {
  mean(abs(M1 - M2))
}

#' Compare the A matrices
#'
#' This function compute errors between the two A matrices,
#'   after guessing the permutation from the minimum MAE.
#'
#' @param A_r The real A matrix used for the simulation..
#' @param A_est The calculated A matrix.
#'   Permutation is guessed from the minimum MAE.
#'
#'@return This function return the matrix with noise.
#' @export
compare_A <- function(A_r, A_est) {
  
  N <- ncol(A_r)
  K <- nrow(A_r)
  stopifnot(K > 1)
  
  stopifnot(ncol(A_est) == N)
  stopifnot(nrow(A_est) < ncol(A_est))
  stopifnot(nrow(A_est) <= 10)
  stopifnot(!anyNA(A_est))
  
  # if not supplying enough types (make sure that {nrow(A_est) >= K})
  if (nrow(A_est) < K) A_est <- rbind(A_est, matrix(0, K - nrow(A_est), N))
  
  comb <- unlist(
    combinat::combn(nrow(A_est), K, fun = combinat::permn, simplify = FALSE),
    recursive = FALSE
  )
  comb_MAE <- sapply(comb, function(perm) {
    MAE(A_r, A_est[perm, , drop = FALSE])
  })
  
  perm <- comb[[which.min(comb_MAE)]]
  A_estim_perm <- A_est[perm, , drop = FALSE]
  
  return(MAE  = MAE(A_r, A_estim_perm))
}

