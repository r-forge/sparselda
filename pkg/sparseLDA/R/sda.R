
sda <- function (x, ...) UseMethod("sda")

## todo: split the stop arg into two different and have one supercede
## todo: added scaling within the function?

sda.default <- function(x, y, lambda=1e-6, stop, maxIte=100, trace=FALSE, tol=1e-6, ...){
  ##
  ## sda performs Sparse Linear Disciminant Analysis
  ## Solving: argmin{|(y*theta-x*b)|_2^2 + t*|beta|_1 + lambda*|beta|_2^2}
  ##
  ## InPUT:
  ## x      : matrix of n observations down the rows and p variable columns. The
  ##          columns are assumed normalized
  ## y      : matrix initializing the dummy variables representing the groups or a factor
  ## lambda : the weight on the L2-norm for elastic net regression. Default: 1e-6.
  ## stop   : If STOP is negative, its 
  ##          absolute value corresponds to the desired number of variables. If STOP
  ##          is positive, it corresponds to an upper bound on the L1-norm of the
  ##          b coefficients. There is a one to one correspondence between stop
  ##          and t.
  ## maxIte : Maximum number of iterations. Default: 100.
  ## trace  : trace = FALSE turns printing of RSS off and trace = TRUE turns it on.
  ## tol    : Tolerance for the stopping criterion (change in RSS). Default is 1e-6.
  ##
  ## OUTPUT:
  ## $beta    : The sparse loadings
  ## $theta   : Optimal scores
  ## $rss     : Residual Sum of Squares at each itearation
  ##
  ## Author: Line H. Clemmensen, IMM, DTU, lhc@imm.dtu.dk
  ## Based on the elastic net algorithm by Hui Zou and Trevor Hastie
  ##

  ## this is stright from nnet:::formula
  class.ind <- function(cl) {
        n <- length(cl)
        x <- matrix(0, n, length(levels(cl)))
        x[(1:n) + n * (as.vector(unclass(cl)) - 1)] <- 1
        dimnames(x) <- list(names(cl), levels(cl))
        x
    }

  if(is.factor(y))
    {
      classes <- levels(y)
      factorY <- y
      y <- class.ind(y)
    } else {
      classes <- colnames(y)
      factorY <- factor(colnames(y)[apply(y, 1, which.max)])
    }
  if(!is.matrix(x)) x <- as.matrix(x)
  predNames <- colnames(x)
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  K <- dim(y)[2]
  RSS <- 1e6
  RSSold <- Inf
  ite <- 0


  Dpi <- t(y)%*%y/n ## diagonal matrix of class priors
  Dpi_inv <- diag(1/sqrt(diag(Dpi)))
  theta <- 1/sum(diag(Dpi))*diag(rep(1,K))[,1:K]/K
  Ytheta <- y%*%theta
  if (length(stop)< (K)){
    stop <- rep(stop[1],1,K)
  }
  if (stop[1]<0) sparse <- "varnum" else sparse <- "penalty" 
  Yhat <- matrix(0,n,K)
  b <- matrix(0,p,K)
  rss <- rep(0,maxIte)

  while (abs(RSSold-RSS)/RSS > tol & ite < maxIte){ 
    RSSold <- RSS
    ite <- ite + 1
    ## 1. Estimate beta:    
    for (j in 1:K){
      Yc <- Ytheta[,j] 
      beta<- solvebeta(x, Yc, paras=c(lambda, abs(stop[j])),sparse=sparse)
      b[,j] <- t(beta)
      Yhat[,j] <- x%*%b[,j] 
    }    
    
    ## 2. Optimal scores: (balanced Procrustes problem)
    B <- t(y)%*%Yhat
    sb <- svd(B,nu=K,nv=K)
    theta.old <- theta
    theta <- Dpi_inv%*%sb$u%*%t(sb$v)
    Ytheta <- y%*%theta
    RSS <- sum((Ytheta-Yhat)*(Ytheta-Yhat)) + lambda*sum(b*b)
    rss[ite] <- RSS
    if (trace){ 
      cat('ite: ', ite, ' ridge cost: ', RSS, ' |b|_1: ', sum(abs(b)),'\n')
    }
  }
  rss <- rss[1:ite]
## remove trivial directions
  Ik <- sb$d > 1e-6
  M <- sum(Ik)
  theta <- theta[,1:M,drop = FALSE]
  Ytheta <- y%*%theta
  b <- b[,1:M,drop = FALSE]
  
  for (j in 1:M){
    Yc <- Ytheta[,j]
    beta<- solvebeta(x, Yc, paras=c(lambda, abs(stop[j])),sparse=sparse)
    b[,j] <- t(beta)
    Yhat[,j] <- x%*%b[,j]
  }    
  if (trace){
    RSS <- sum((Ytheta-Yhat)*(Ytheta-Yhat)) + lambda*sum(b*b)
    cat('final update, ridge cost: ', RSS, ' |b|_1: ', sum(abs(b)),'\n')
  }

  notZero <- apply(b, 1, function(x) any(x != 0))
  b <- b[notZero,]
  origP <- ncol(x)
  x <- x[, notZero, drop = FALSE]
  varNames <- colnames(x)
  
  sl <- x %*% b
  colnames(sl) <- paste("score", 1:ncol(sl), sep = "")
  lobj<-lda(sl,factorY, ...)
  
  ## todo save fitted and probs (and data?) for predict

  structure(
            list(call = match.call(),
                 beta = b,
                 theta = theta,
                 varNames = varNames,
                 varIndex = which(notZero),
                 origP = origP,
                 rss = rss[1:ite],
                 fit = lobj,
                 classes = classes,
                 lambda = lambda,
                 stop = stop),
            class = "sda")
}

print.sda <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")


    cat("lambda =", format(x$lambda, digits = digits),
        "\nstop =", format(x$stop, digits = digits),
        "\n\n")

    top <- if(!is.null(x$varNames)) x$varNames else paste("Predictor", x$varIndex, sep = "")
    varOrder <- if(is.matrix(x$beta)) order(apply(abs(x$beta), 1, sum)) else order(abs(x$beta))
    top <- top[varOrder]
    top <- top[1:min(5, length(top))]
    top <- paste(top, collapse = ", ")
    
    if(length(x$beta) > 5)
      {
        cat("Top 5 predictors (out of ",
            length(x$varIndex),
            "):\n\t",
            top,
            sep = "")
      } else {
        cat("Predictors:\n\t",
            top,
            "\n",
            sep = "")
      }
    invisible(x)
  }

predict.sda <- function(object, newdata = NULL, ...)
  {
    if(!is.matrix(newdata)) newdata <- as.matrix(newdata)
    if(!is.null(object$varNames))
      {
        newdata <- newdata[, object$varNames, drop = FALSE]
      } else {
        if(ncol(newdata) != object$origP) stop("dimensions of training and testing X different")
        newdata <- newdata[, object$varIndex, drop = FALSE]
      }
    x <- newdata %*% object$beta
    predict(object$fit, newdata = x, ...)
  }



## todo: plot and summary 



