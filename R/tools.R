##' KL distance function
##' D(spec1 || spec2) = sum(spec1 * log(spec1 / spec2))
##' @param spec1 spectral signal 1
##' @param spec2 spectral signal 2 at identical wavelengths
##' 
##' @export
KLd <- function(spec1,spec2){

    if(length(spec1)==length(spec2)){
        return(sum(spec1 * log(spec1 / spec2)))
    } else {

        stop("spectra not of equal length")
    }
}


#' Delchambre's Weighted PCA function
#'
#' @param X Dataset matrix with columns corresponding to variables and rows
#' to oberservations (as is standard in R, but in constrast to Delchambre's code).
#' The matrix will be scaled after the weighted mean is subtracted in the code
#' no need to do this prior as in the original Delchambre code.
#' @param W Weight matrix, similar to X, but containing data weights. 
#' @param ncomp number of principle components to retrieve
#' @param niter number of iterations, a larger number allows more
#' precise estimation of the princomp. 
#' @param nrefine parameters that allows refinment. Advice is to
#' not touch this parameter unless you know what you are doing. 
#' @param xi Regularization factor: controls the influence of
#' uncertain (rare) variables.
#' @return Principle components arranged in columns
#'  \itemize{
#' \item [1] = P: Principle components
#' \item [2] = C: Coefficients 
#' \item [3] = PC: Projected original data (PC = P%*%C)
#' \item [4] = variance explained
#' }

#' @author Ludovic Delchambre (original matlab), Marco D. Visser (R translation).
#' @examples {
#'  set.seed(2020)
#'  X<-matrix(rnorm(10*10,0,4),ncol=10)
#'  W<-matrix(runif(10*10,0.999,1),ncol=10) ## pretty certain data
#'
#' 
#' pca<-princomp(scale(X))
#' dpca<-dPCA(scale(X),W)[[2]]
#'
#' 
#' } 
#' @export
#'
dPCA <- function(D,W,xi=rep(1,ncol(W)),
                 eigendecom=TRUE){

    D <- t(as.matrix(D))
    W <- t(as.matrix(W))

    if(!is.numeric(D)|!is.numeric(W)) stop("Data can only contain numeric values!")
    if(any(apply(W,2,function(X) all(X<=0)))) stop("weights must be > zero for observations")

    nobs <- ncol(D)
    nvar <- nrow(D)
    
    ## create matrix X centered on weighted mean
    wMu <- rowMeans(W*D)/rowMeans(W)
    X <- D - wMu

    ## scale by weighted sd in each dimension
    sdW <- sqrt(apply(W^2 * X^2,1,sum)/rowSums(W^2))
    for(i in 1:nrow(X)) X[i,] <- X[i,]/sdW[i] # devides each row by sdW

    
    ## weighted covariance matrix
    if(any(xi!=1)){
        ws <- rowSums(W)
        covar <- ((ws%*%t(ws))^xi)*((X*W)%*%t(X*W))/(W%*%t(W))
    } else {
        covar <- (W*X)%*%t(W*X) / (W%*%t(W)) 
    }

    covar[is.nan(covar)] <- 0    

    
    ##  calculate eigenvalues and eigenvectors
    EV <- eigen(covar) # ordered from the dominant lambda and decending 

    ## Calculate explained variance
    chi2 <- numeric(length=nvar) # explained variance

    ncomp <- ncol(EV$vectors)
    
    for (p in 1:ncomp){  

        ## build  principal components matrix
        P <- EV$vectors[, 1:p]
        
        ## Solve for coefficients matrix
        C <- matrix(data=NA, nrow=p, ncol=nobs)
        
        for (i in 1:nobs){
            u <- diag(W[, i]^2) 
            C[, i] <- solve(t(P)%*%u%*%P, t(P)%*%u%*%X[, i])
        }

        ## Use coefficients such to projections original data ( X = P%*%C)
        PC <- P%*%C

        ## explained variance
        chi2[p] <- sum(sum((W*PC)^2)) / sum(sum((W*X)^2))

    }

    ## order based on X^2 statistic 
    P <- EV$vectors
    C <- matrix(data=NA, nrow=ncol(P), ncol=nobs)

    for (i in 1:nobs){
        w <- diag(W[,i]^2) 
        C[, i] <- solve(t(P)%*%w%*%P, t(P)%*%w%*%X[, i])
    }

    PC <- P%*%C
    
    return(list(princomp=P,coefficients=C,proj=PC,varexp=chi2))
}




