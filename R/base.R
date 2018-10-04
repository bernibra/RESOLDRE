#' Informed Randomizations
#'
#' This functions performs different type of randomization of an adjacency matrix given either a probability or a correlation matrix.
#' @usage resoldre(mat = NULL, pmat = NULL, cormed = NULL, randomizations = 1, bipartite = FALSE, seed=NULL, nperm=1, perspective="rows", degree=NULL, maxit=40, ss=0.1, tolpql=10^-6, maxitpql=200, f=NULL, ...)
#' @param mat
#' an adjacency matrix to be randomized.
#' @param pmat
#' an optional probability matrix defining the probability of encountering any of the possible interactions in the adjacency matrix. This should be unset if a correlation matrix is specified. If 'pmat' and 'cormed' are unset, this probability matrix is considered to be an all-ones matrix.
#' @param cormed
#' a correlation matrix used to bias the randomization. The dimensions of the matrix need to agree with the 'perspective' set for the randomization as well as the adjacency matrix 'mat'. Perfectly correlated pairs should contain the lowest elements of the matrix. If 'pmat' is set, 'cormed' will be ignored.
#' @param randomizations
#' an integer value specifying total number of randomizations.
#' @param bipartite
#' a logical value specifying whether or not 'mat' describes a unipartite or bipartite graph. It is set as FALSE by default.
#' @param seed
#' a single value setting a seed for the whole randomization process.
#' @param nperm
#' an integer value specifying the number of link swaps made by each randomization.
#' @param perspective
#' this parameter can be set to either "rows" or "columns" . It specifies the perspective taken for the estimation of the probability matrix if 'pmat' is NULL and 'cormed' is set (see details).
#' @param degree
#' this parameter determines if the randomizations preserve the degree of "columns", "rows", "both" or "sample".
#' @param maxit
#' a control parameter dictating the maximum number of iterations in the optimization.
#' @param ss
#' initial estimates of the parameter for the random effect that scales the variance (it's kind of unnecessary).
#' @param tolpql
#' a control parameter dictating the tolerance for convergence in the PQL estimates of the mean components of the binomial GLMM (it's kind of unnecessary).
#' @param maxitpql
#' a control parameter dictating the maximum number of iterations in the PQL estimates of the mean components of the binomial GLMM (it's kind of unnecessary).
#' @param f
#' a function to be applied to each randomized matrix.
#' @param ...
#' optional arguments to f.
#' @details
#' Given an n x m adjacency matrix and an m x m correlation matrix relating its columns, resoldre will estimate the probability of encountering any of possible interactions based on the columns correlation matrix if 'perspective' is set as "columns". Similarly, if 'cormed' is otherwise an n x n correlation matrix relating the columns of the adjacency matrix, 'perspective' needs to be set as "rows".
#' @return The result of the randomization in some way...
#' @export
resoldre <- function(mat = NULL, pmat = NULL, cormed = NULL, randomizations = 1,
                     bipartite = FALSE, seed=NULL, nperm=1, perspective="rows", degree="sample",
                     maxit=40, ss=0.1, tolpql=10^-6, maxitpql=200, f=NULL, ...){


  # This was originally an option. It basically chooses between different ways of computing the rewiring proability
  fprob <- 0
  ##

  #Input checks
  if(is.null(mat)){stop("You must specify a target matrix 'mat' to be randomized.")}
  if(perspective!="rows" & perspective!="columns"){
    stop(gettextf("'perspective' should be one of %s", paste(dQuote(c("rows","columns")), collapse = ", ")), domain = NA)
  }

  if (!is.null(degree)) {
    if(degree!="rows" & degree!="columns" & degree!="both" & degree!="sample"){
      stop(gettextf("'degree' should be one of %s", paste(dQuote(c("rows","columns", "both", "sample")), collapse = ", ")), domain = NA)
    }
  }

  if(is.matrix(mat)){
    if(any(mat!=1 & mat!=0)){stop("The target matrix 'mat' must be an interaction matrix.")}
    if(!any(mat>0)){stop("'mat' can't be a null matrix.")}
    if(nrow(mat)!=ncol(mat) & !bipartite){stop("'mat' must be a square matrix unless you specify 'bipartite=TRUE'.")}
  }else{
      stop("'mat' must be a matrix")
  }
  if(is.null(cormed) & is.null(pmat)){
    pmat <- matrix(1, nrow(mat), ncol(mat))
  }else{
    if(!is.null(cormed) & !is.null(pmat)){
      if(is.matrix(pmat) & nrow(pmat)==nrow(mat) & ncol(pmat)==ncol(mat)){
        if (any(pmat>0)) {
          warning("'cormed' will be ignored because you have already provided 'pmat'.")
        } else{
          warning("'pmat' can't be a null matrix.")
        }
      }else{
        stop("'pmat' must be a matrix with the same dimensions as 'mat'.")
      }
    }
    if(is.null(cormed)){
      if(!is.matrix(pmat) | nrow(pmat)!=nrow(mat) | ncol(pmat)!=ncol(mat)){
        stop("'pmat' must be a matrix with the same dimensions as 'mat'.")
      } else{
        if (!any(pmat>0)) {warning("'pmat' can't be a null matrix.")}
      }
    }else{
      if(is.matrix(cormed) & nrow(cormed)==ncol(cormed)){
        if ((perspective=="rows" & nrow(cormed)!=nrow(mat)) | (perspective=="columns" & nrow(cormed)!=ncol(mat))){
          stop("The dimensions of 'cormed' must agree with 'perspective'.")
        }else{
          pmat <- probability_estimation(mat=mat, vcv=cormed, degree=degree, perspective=perspective, maxit=maxit, ss = ss, tolpql = tolpql, maxitpql = maxitpql)
        }
      }else{
        stop("'cormed' must be a square matrix.")
      }
    }
  }
  if(degree=="sample"){
    stop("'pmat' must be a matrix with the same dimensions as 'mat'.")
  }

  if(!is.null(seed)){
    set.seed(seed)
  }


  #######################################################################################################
  # Randomizations
  #######################################################################################################

  #I should do this for all the probability functions
  pmat <- pmat*1./max(pmat)

  if(degree=="sample"){

    random <- matrix(0,nrow(mat), ncol(mat))
    links <- which(random==0, arr.ind = TRUE)
    pmat <- (rowSums(mat) %*% t(colSums(mat))) * pmat

    if(!is.null(f)){
      results <- lapply(1:randomizations, function(x) f(psample(pmat, sum(mat), random, links),...))
    }else{
      results <- lapply(1:randomizations, function(x) psample(pmat, sum(mat), random, links))
    }

  }else{

    links <- which(mat>0, arr.ind = TRUE)

    if(bipartite){

      if(degree=="both"){
        type <- 1
        degree <- matrix(0,0,0)
      }else if(degree=="rows"){
        type <- 2
        degree <- colSums(mat)
      }else{
        type <- 3
        degree <- rowSums(mat)
      }

      links <- links-1
      if(!is.null(f)){
        results <- lapply(1:randomizations, function(x) f(randomize(mat, pmat, links, NULL, nperm, type, fprob, degree = degree),...))
      }else{
        results <- lapply(1:randomizations, function(x) randomize(mat, pmat, links, NULL, nperm, type, fprob, degree = degree))
      }

    }else{

      type <- 0
      if(!is.null(degree)){
        warning("'degree' will be ignored because 'bipartite=FALSE'")
      }

      unilinks <- c()
      bilinks <- c()

      for(i in 1:nrow(links)){
        if(links[i,1]==links[i,2]){
          next
        }else if(mat[links[i,1],links[i,2]]==mat[links[i,2],links[i,1]] & links[i,1]<links[i,2]){
          bilinks <- rbind(bilinks,c(links[i,1],links[i,2]))
        }else if(mat[links[i,1],links[i,2]]!=mat[links[i,2],links[i,1]]){
          unilinks <- rbind(unilinks,c(links[i,1],links[i,2]))
        }
      }

      if(length(unilinks)==0){unilinks <- NULL}else{unilinks <- unilinks-1}
      if(length(bilinks)==0){bilinks <- NULL}else{bilinks <- bilinks-1}

      if(!is.null(f)){
        results <- lapply(1:randomizations, function(x) f(randomize(mat, pmat, unilinks, bilinks, nperm, type, fprob, degree = matrix(0,0,0)),...))
      }else{
        results <- lapply(1:randomizations, function(x) randomize(mat, pmat, unilinks, bilinks, nperm, type, fprob, degree = matrix(0,0,0)))
      }
    }
  }
  ######################################################################################################

  return(results)
}


