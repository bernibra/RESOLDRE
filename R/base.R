#' Main program
#'
#' Informed randomizations main script
#' @param m A matrix to be randomized.
#' @param cormed A correlation-informed matrix used to bias the randomization.
#' @return The result of the randomization in some way...
#' @export
resoldre <- function(mat = NULL, pmat = NULL, cormed = NULL, randomizations = 1,
                     bipartite = FALSE, seed=NULL, nperm=100, perspective="rows", fprob=0,
                     maxit=40, ss=0.1, tolpql=10^-6, maxitpql=200, f=NULL, ...){

  #Input checks
  if(is.null(mat)){stop("You must specify a target matrix 'mat' to be randomized.")}
  if(perspective!="rows" & perspective!="columns"){
    stop(gettextf("'perspective' should be one of %s", paste(dQuote(c("rows","columns")), collapse = ", ")), domain = NA)
  }
  if(is.matrix(mat)){
    if(any(mat!=1 & mat!=0)){stop("The target matrix 'mat' must be an interaction matrix.")}
    if(nrow(mat)!=ncol(mat) & !bipartite){stop("'mat' must be a square matrix unless you specify 'bipartite=TRUE'.")}
  }else{
      stop("'mat' must be a matrix")
  }
  if(is.null(cormed) & is.null(pmat)){
    pmat <- matrix(1, nrow(mat), ncol(mat))
  }else{
    if(!is.null(cormed) & !is.null(pmat)){
      if(is.matrix(pmat) & nrow(pmat)==nrow(mat) & ncol(pmat)==ncol(mat)){
        warning("'cormed' will be ignored because you have already provided 'pmat'.")
      }else{
        stop("'pmat' must be a matrix with the same dimensions as 'mat'.")
      }
    }
    if(is.null(cormed)){
      if(!is.matrix(pmat) | nrow(pmat)!=nrow(mat) | ncol(pmat)!=ncol(mat)){
        stop("'pmat' must be a matrix with the same dimensions as 'mat'.")
      }
    }else{
      if(is.matrix(cormed) & nrow(cormed)==ncol(cormed)){
        if ((perspective=="rows" & nrow(cormed)!=nrow(mat)) | (perspective=="columns" & nrow(cormed)!=ncol(mat))){
          stop("The dimensions of 'cormed' must agree with 'perspective'.")
        }else{
          pmat <- probability_estimation(mat=mat, vcv=cormed, perspective = perspective, maxit = maxit, ss = ss, tolpql = tolpql, maxitpql = maxitpql)
        }
      }else{
        stop("'cormed' must be a square matrix.")
      }
    }
  }

  #I should do this for all the probability functions
  pmat <- pmat*1./max(pmat)


  if(!is.null(seed)){
    set.seed(seed)
  }


  links <- which(mat>0, arr.ind = TRUE)

  if(bipartite){
    links <- links-1
    if(!is.null(f)){
      results <- lapply(1:randomizations, function(x) f(randomize_bipartite(mat, pmat, links, nperm, fprob),...))
    }else{
      results <- lapply(1:randomizations, function(x) randomize_bipartite(mat, pmat, links, nperm, fprob))
    }
  }else{

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
      results <- lapply(1:randomizations, function(x) f(randomize_unipartite(mat, pmat, unilinks, bilinks, nperm, fprob),...))
    }else{
      results <- lapply(1:randomizations, function(x) randomize_unipartite(mat, pmat, unilinks, bilinks, nperm, fprob))
    }
  }

  return(results)
}


