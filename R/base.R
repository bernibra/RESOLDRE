#' Main program
#'
#' Informed randomizations main script
#' @param m A matrix to be randomized.
#' @param cormed A correlation-informed matrix used to bias the randomization.
#' @return The result of the randomization in some way...
#' @export
resoldre <- function(mat = NULL, pmat = NULL, cormed = NULL, randomizations = 1, verbose = FALSE){
  if(is.null(mat)){stop("You must specify a target matrix mat to be randomized")}
  if(is.matrix(mat)){
    if(any(mat!=1 & mat!=0)){
      stop("The target matrix mat must be an interaction matrix")
    }
  }else{
      stop("mat must be a matrix")
  }
  if(is.null(cormed) & is.null(pmat)){
    pmat <- matrix(1, nrow(mat), ncol(mat))
  }else{
    if(!is.null(cormed) & !is.null(pmat)){
      if(is.matrix(pmat)){
        warning("cormed will be ignored because you have already provided pmat")
      }else{
        stop("pmat must be a matrix")
      }
    }
    if(is.null(cormed)){
      if(!is.matrix(pmat)){stop("pmat must be a matrix")}
    }else{
      if(is.matrix(cormed)){
#         pmat <- pglmm(mat, cormed)
        pmat <- matrix(1, nrow(mat), ncol(mat))
      }else{
        stop("cormed must be a matrix")
      }
    }
  }
  return(list(mat=mat, pmat=pmat))
}
