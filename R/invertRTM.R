#' invert a requested RTM (internal function)
#'
#' List of aliases: prospect5, prospectd
#' @param modReq  model request object built in bRTM
#' @param pars the required parameters (vector or list),
#' and newdata
#' @return prediction from the requested model
invertRTM <- function(pars){
    result <- UseMethod("irtm",pars)
    return(result)
}

## leaf models

irtm.prospect5 <- function(pars){
  predict.bayes(model5,nn5b,eigenR5,rescaler5,plsr5,newdata=pars)
}

irtm.prospectd <- function(pars){
  predict.bayes(modeld,nnd,eigenRd,rescalerD,plsrd,newdata=pars)
}


################################################################################
## start: Supporting Tools
################################################################################
## plsr prediction function
## @param object fitted pls object
## @param newdata a dataframe to predict from
## @param ncomp number of components to use in prediction
## @param comps number of components used in prediction
predict.pls <- function (object, newdata, ncomp = 1:object$ncomp,
                         comps=object$ncomp) {
    if (missing(newdata) || is.null(newdata)){
        stop("data not supplied")
    }  else if (is.matrix(newdata)) {
        if (ncol(newdata) != length(object$Xmeans)) 
            stop("'newdata' does not have the correct number of columns")
        newX <- newdata
    }

    nobs <- dim(newX)[1]

    ## prep coefficients
    beta <- object$coefficients[, , comps, drop = FALSE]
    g1 <- which(comps > 1)
    beta[, , g1] <- beta[, , g1, drop = FALSE] - object$coefficients[, 
            , comps[g1] - 1, drop = FALSE]
    B <- rowSums(beta, dims = 2)
    B0 <- object$Ymeans - object$Xmeans %*% B
    pred <- newX %*% B + rep(B0, each = nobs)
    return(pred)
}

## MANN prediction function
## @param object fitted multivariate neural net model object
## @param newdata a dataframe to predict from
## @param eigenR eigenvectors for data reduction
## @ncomp number of components object needs
predict.nn <- function(object, newdata, eigenR,rescaler, ncomp = 40){

  if(nrow(newdata)==1){

    scaledR <- (newdata-rescaler[["scaled:center"]])/
      rescaler[["scaled:scale"]]
    pcaR <- scaledR%*%(eigenR$vectors)
    pcaR <- pcaR[,1:ncomp]
    biases <- matrix(1,ncol=ncol(object$w2))
    pred <- feedforward(c(biases,pcaR),object$w1,object$w2)$output

  } else {

    scaledR <- t(apply(newdata,1,
                       function(X) (X-rescaler[["scaled:center"]])/
                                   rescaler[["scaled:scale"]]))
    pcaR <- scaledR%*%(eigenR$vectors)
    pcaR <- pcaR[,1:ncomp]
    biases <- sapply(1:ncol(object$w2),function(X) rep(1,nrow(pcaR)))
    pred <- feedforward(cbind(biases,pcaR),object$w1,object$w2)$output
  }

  return(pred)
}


## Bayes prediction function
## @param object fitted pls object
## @param nnfit neural net fit
## @param plsfit PLS fit
## @param newdata a dataframe (spectra) to predict from
predict.bayes <- function(object,nnfit,eigenR,rescaler,plsfit,newdata){

  ## PLS predictions
  predpls <- predict.pls(plsfit,newdata) ## pls checks newdata

  ## NN predictions
  predNN <- predict.nn(nnfit, newdata, eigenR, rescaler,ncomp = 40)
  colnames(predNN) <- colnames(predpls)

  ord <- colnames(object)
  bias <- matrix(rep(object["b0",],nrow(predpls)),ncol=ncol(predpls),byrow=TRUE)
  uncert <- matrix(rep(object["sigma",],nrow(predpls)),ncol=ncol(predpls),byrow=TRUE)

  pred <- bias  + predNN[,ord]%*%diag(object["bnn",]) + predpls[,ord]%*%diag(object["bpls",])
  colnames(pred) <- ord

  mu <- pred
  low.ci <- pred-2*uncert
  upp.ci <- pred+2*uncert
  low.ci[low.ci<0] <- 0
  upp.ci[upp.ci<0] <- 0
  mu[mu<0] <- 0
  list(mu=mu,lower.ci=low.ci,upper.ci=upp.ci)
}

################################################################################
## basic NN functions for a general MANN
################################################################################
## activation function S3 method
act <- function(object,...){

  UseMethod("act")

}

act.linear <- function(object) object
act.sigmoid <- function(object) 1 / (1 + exp(-object))
act.exponential <- function(object) exp(object)

feedforward <- function(X, w1, w2,a.hid="sigmoid",a.out="linear") {
  ## forward propagation
  ## X design matrix
  ## w1, w2 = weight matrices
  z1 <- X %*% w1
  class(z1) <- a.hid
  h <- act(z1)
  z2 <- cbind(1, h ) %*% w2 ## add the bias term again
  class(z2) <- a.out
  list(output = act(z2), h = h)
}


################################################################################
## end: Supporting Tools
################################################################################

################################################################################
## start: Inversion datasets


#' Bayesian fitted weight matrix for PROSPECT5
#'
#' Weight coefficients for neural network
#' and plsr predictions .
#'
#' @docType data
#' @keywords datasets
#' @name model5
## @usage data(model5) ## not public
NULL

#' fitted weight matrix for PROSPECT5
#'
#' Weight matrices for a fit neural network
#' on simulated data from PROSPECT5.
#'
#' @docType data
#' @keywords datasets
#' @name nn5b
## @usage data(nn5b) ## not public
NULL

#' fitted PLSR for PROSPECT5
#'
#' A partial least squares model fit
#' on simulated data from PROSPECT5.
#'
#' @docType data
#' @keywords datasets
#' @name plsr5
## @usage data(plsr5) ## not public
NULL

#' eigen decomposition for PROSPECT5
#'
#' data reduction used on simulated
#' data from PROSPECT5 (for NN and PLSR)
#'
#' @docType data
#' @keywords datasets
#' @name eigenRb
## @usage data(eigenRb) ## not public
NULL

#' Bayesian fitted weight matrix for PROSPECTD
#'
#' Weight coefficients for neural network
#' and plsr predictions .
#'
#' @docType data
#' @keywords datasets
#' @name modeld
## @usage data(modeld) ## not public
NULL

#' fitted weight matrix for PROSPECTD
#'
#' Weight matrices for a fit neural network
#' on simulated data from PROSPECTD.
#' 
#' @docType data
#' @keywords datasets
#' @name nnd
## @usage data(nnd) ## not public
NULL

#' fitted PLSR for PROSPECTD
#'
#' A partial least squares model fit
#' on simulated data from PROSPECTD.
#'
#' @docType data
#' @keywords datasets
#' @name plsrd
## @usage data(plsrd) ## not public
NULL

#' eigen decomposition for PROSPECTD
#'
#' data reduction used on simulated
#' data from PROSPECTD (for NN and PLSR)
#'
#' @docType data
#' @keywords datasets
#' @name eigenRd
## @usage data(eigenRd) ## not public
NULL


