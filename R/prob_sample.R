#' Sample randomization
#'
#' Probabilistic randomization
#' @return Random matrix
#' @export
psample <- function(pmat, n, random, links){
  random[links[sample(1:nrow(links), n, replace=FALSE, prob=pmat[links]),]] <- 1
  return(random)
}
