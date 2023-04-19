## All generic S3 methods

#' Plot RTM return spectra vs. wavelength
#'
#' @param x predictions from an RTM
#' @param \dots additional plot arguments
#' @return plots to the device a ccrtm standard spectra plot based on 
#' the function call returned from fRTM. 
#' @export
plot.rtm.spectra <- function(x, ...) {

    wl <- attr(x,"wavelength")
    main <- attr(x,"main")
    models <- attr(x,"models") 
    specnums <- sapply(c("rho","tau"), grepl, x=main)
    if(!is.null(ncol(specnums))) {
        specnums <- apply(specnums,2,which)
    }
    
    specnames <- c("reflectance", "transmission")[specnums]

    if (is.null(ncol(x))) {
        plot(x = wl, y = as.numeric(x),
             ylab=specnames,type="l",lwd=2,
             xlab="wavelength (nm)", ...)
        
    } else{

        if(ncol(x)>2) specnames <- colnames(x)
        
        matplot(x = wl, y = x,
                ylab="Spectra",
                xlab="wavelength (nm)",
                col=rainbow(ncol(x)),
                type="l",lwd=2,
                ...)
        pos <- ifelse(max(wl)<1500,"topleft", "topright")
        legend(pos, legend=specnames,
               col=rainbow(ncol(x)),bty="n",
               lty=1,lwd=1.5)

      }
}


#' RTM generic print function
#'
#' @param x predictions from an RTM
#' @param \dots additional plot arguments
#' @return prints the standard information from a simulated 
#' ccrtm spectra plot  
#' @export
print.rtm.spectra <- function(x, ...) {

    wl <- attr(x,"wavelength")
    main <- attr(x,"main")
    models <- attr(x,"models") 
    specnums <- sapply(c("rho","tau"), grepl, x=main)
    if(!is.null(ncol(specnums))) {
        specnums <- apply(specnums,2,which)
    }

    specnames <- c("reflectance", "transmission")[specnums]

    
    if(length(models)>1) models <- paste(models,collapse=", ")
    
    if(length(specnums)>1){
        
        specnames <- paste(specnames, collapse=" and ")
    }
   
    cat("RTM predicted spectra: ",specnames, "\n")
    cat("Generating model(s): ",models, "\n")
    cat("Wavelength range ",paste(range(wl),collapse="-"), "(nm) \n")
           
        
}


#' RTM inversion
#'
#' @param x predictions from an RTM
#' @param \dots additional plot arguments
#' @return prints the inverted parameters
#' @export
print.rtm.inversion <- function(x, ...) {

  models <- attr(x,"models") 

  cat("RTM inversions: ",nrow(x$mu), "\n")
  cat("inverted  model(s): ",models, "\n")
}



#' RTM inversion summary
#'
#' @param x predictions from an RTM
#' @param \dots additional plot arguments
#' @return summarizes the inverted parameters
#' @export
summary.rtm.inversion <- function(x, ...) {

  models <- attr(x,"models")

  mu <- colMeans(x$mu)
  mu.l <- colMeans(x$lower.ci)
  mu.u <- colMeans(x$upper.ci)

  xformated <- paste0("Estimate ",colnames(x$mu)," : " ,signif(mu,3)," (CI: ",signif(mu.l,2),"-",signif(mu.u,2),")\n")
  
  xhead <- paste0(signif(head(x$mu),3)," (CI: ",
                  signif(head(x$lower.ci),2),"-",signif(head(x$upper.ci),2),")")
 # colnames(xhead) <- colnames(x$mu)

  cat("\n")
  cat("------------------------------------------------------\n")
  cat("Inverted  model(s): ",models, "\n")
  cat("------------------------------------------------------\n")
  cat("Inverted (mean) parameters and inversion uncertainty: \n",xformated, "\n\n")
  cat("Head of inverted values matrix: \n ")
  print(head(x$mu))
  cat("\n A total of",nrow(x$mu), "RTM inversions performed. \n\n")

  invisible(data.frame(mean=mu,lower=mu.l,upper=mu.u))
}




