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


#' Plot RTM return spectra vs. wavelength
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
