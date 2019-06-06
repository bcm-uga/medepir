#'Function to run MeDeCom.
#'
#' This function run MeDeCom and select A and T.
#'
#'@param D The matrix patients*probes.
#'@param nbcell The number of cell types.
#'@param lambdas Values of parameter lambda to be tested.
#'
#'@return This function return a list with the A and the T matrices computed.
#'
#'@example 
#'
#'@export

MDC = function(D, nbcell = 5, lambdas = c(0,10^(-5:-1))){
  #Run MeDeCom for all the lambda values
  result_medecom = MeDeCom::runMeDeCom(D, nbcell, lambdas, NINIT = 10, NFOLDS = 10, ITERMAX = 300, NCORES = 5)
  #Extraction of each results
  res_lambdas = round(as.numeric(MeDeCom::getStatistics(result_medecom, nbcell, lambdas)))
  #Choice of the best lambda
  lambda_max = which.min(res_lambdas)
  best_lambda = lambdas[lambda_max]
  #Extraction of T and A matrices
  A_est = MeDeCom::getProportions(result_medecom, K = nbcell, best_lambda)
  T_est = MeDeCom::getLMCs(result_medecom, K = nbcell, best_lambda)  
  return(list(A = A_est, T = T_est))
}
