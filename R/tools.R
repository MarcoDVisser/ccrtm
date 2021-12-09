##' Kullback-Lieber divergence function
##' D(spec1 || spec2) = sum(spec1 * log(spec1 / spec2))
##' @param spec1 spectral signal 1
##' @param spec2 spectral signal 2 at identical wavelengths
##' 
##' @return the KL divergence between the vector inputs 
##' @export
KLd <- function(spec1,spec2){

    if(length(spec1)==length(spec2)){
        return(sum(spec1 * log(spec1 / spec2)))
    } else {

        stop("spectra not of equal length")
    }
}





