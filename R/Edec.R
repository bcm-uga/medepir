#'Function to run Edec.
#'
#' This function run Edec. and select A and T.
#'
#'@param D The matrix patients*probes.
#'@param nbcell The number of cell types
#'@param infloci The list of informative probes.
#'
#'@return This function return a list with the A and the T matrices computed.
#'
#'@example 
#'
#'@export

Edec = function(D, nbcell = 5, infloci){
  
  if (!requireNamespace("Edec", quietly = TRUE))
    stop("Please install package {Edec} (https://github.com/BRL-BCM/EDec).")
  
  #Application de la méthode pour cette liste
  deconv = EDec::run_edec_stage_1(meth_bulk_samples = D, 
                            informative_loci = infloci, num_cell_types = nbcell)
  #Extraction des matrices T et A
  T_est = as.matrix(deconv$methylation)
  A_est = as.matrix(deconv$proportions)
  A_est = t(A_est)
  return(list(A = A_est, T = T_est))
}
