#'Function to run RefFreeEwas with euclidean dist.
#'
#' This function run RefFreeCellMix and select A and T.
#'
#'@param D The matrix patients*probes.
#'@param nbcell The number of cell types
#'
#'@return This function return a list with the A and the T matrices computed.
#'
#'@export
RFE <- function(D, nbcell = 5){
  results = RefFreeEWAS::RefFreeCellMix(D, K = nbcell, iters = 9)
  #Extraction des matrices T et A
  T_est = results$Mu
  A_est = t(results$Omega)
  return(list(A = A_est, T = T_est))
}


#'Function to run RefFreeEwas "V2" with an svd initialization.
#'
#' This function run RefFreeCellMix and select A and T.
#'
#'@param D The matrix patients*probes.
#'@param nbcell The number of cell types
#'
#'@return This function return a list with the A and the T matrices computed.
#'
#'@export
RFE_with_svd <- function(D, nbcell = 5){
  init = RefFreeEWAS::RefFreeCellMixInitializeBySVD(D)
  mu = init[,1:nbcell]

  results = RefFreeEWAS::RefFreeCellMix(D, mu0 = mu,  K = nbcell, iters=9)
  #Extraction des matrices T et A
  T_est = results$Mu
  A_est = t(results$Omega)
  return(list(A = A_est, T = T_est))
}
