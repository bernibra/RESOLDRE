#' Main program
#'
#' Informed randomizations main script
#' @param m A matrix to be randomized.
#' @param cormed A correlation-informed matrix used to bias the randomization.
#' @return The result of the randomization in some way...
#' @export
resoldre <- function(mat = NULL, pmat = NULL, cormed = NULL, randomizations = 1, verbose = FALSE, seed=NULL, nperm=100, fprob=0){

  #Input checks
  bipartite <- FALSE
  if(is.null(mat)){stop("You must specify a target matrix mat to be randomized.")}
  if(is.matrix(mat)){
    if(any(mat!=1 & mat!=0)){stop("The target matrix mat must be an interaction matrix.")}
    if(nrow(mat)!=ncol(mat)){bipartite <- TRUE}
  }else{
      stop("mat must be a matrix")
  }
  if(is.null(cormed) & is.null(pmat)){
    pmat <- matrix(1, nrow(mat), ncol(mat))
  }else{
    if(!is.null(cormed) & !is.null(pmat)){
      if(is.matrix(pmat)){
        if(nrow(pmat)==nrow(mat) & ncol(pmat)==ncol(mat)){
          warning("cormed will be ignored because you have already provided pmat.")
        }else{
          stop("pmat must have the same dimensions as matrix.")
        }
      }else{
        stop("pmat must be a matrix.")
      }
    }
    if(is.null(cormed)){
      if(!is.matrix(pmat) | nrow(pmat)!=nrow(mat) | ncol(pmat)!=ncol(mat)){
        stop("pmat must be a square matrix.")
        }
    }else{
      if(is.matrix(cormed) & nrow(cormed)==ncol(cormed)){
#         pmat <- pglmm(mat, cormed)
        pmat <- matrix(1, nrow(mat), ncol(mat))
      }else{
        stop("cormed must be a square matrix.")
      }
    }
  }

  #I should do this for all the probability functions
  if (fprob==0){
    pmat <- pmat*1./max(pmat)
  }

  if(!is.null(seed)){
    set.seed(seed)
  }


  #Input setup
  links <- which(mat>0, arr.ind = TRUE)

  if(bipartite){

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

  }

  newmat <- randomize(mat, unilinks, bilinks, 1)


  return(list(mat=newmat, pmat, unilinks=univec, bilinks=bilinks, nperms, fprob))
}


