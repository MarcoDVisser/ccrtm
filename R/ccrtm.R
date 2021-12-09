#' ccrtm: Coupled Chain Radiative Transfer Models.
#' 
#' A collection of radiative transfer models that can form
#' a coupled chain to model radiative transfer across multiple
#' spatial scales from leaf to canopy to stand.
#'
#' Currently implemented models:
#'
#'  \itemize{
#'
#' \item [1] = PROSPECT 5, 5B and D
#' \item [2] = FOURSAIL, and FOURSAIL2 
#' \item [3] = FLIM 
#'
#' }
#'
#' Currently being tested or to be implemented models
#'
#'  \itemize{
#'
#' \item [1] = LIBERTY, PROCOSINE
#' \item [2] = INFORM
#'
#' }
#'
#' 
#' @reference Verhoef, W. and Bach, H., 2003. Remote sensing data assimilation
#' using coupled radiative transfer models.
#' Physics and Chemistry of the Earth, Parts A/B/C, 28(1-3), pp.3-13.
#'
#' @reference Feret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
#' Properties Model Separating Photosynthetic Pigments, 
#' Remote Sensing of Environment.
#' 
#' @docType package
#' @author Marco D. Visser 
#' @import Rcpp graphics grDevices stats testthat
#' @importFrom Rcpp sourceCpp
#' @useDynLib ccrtm
#' @name ccrtm
NULL 

#' refractive index and specific absorption coefficient for
#' PROSPECT 5
#'
#' see http://teledetection.ipgp.jussieu.fr/prosail/ for more details
#' on the data.
#' 
#' @section details:
#' ***********************************************************************
#' data_prospect5 (february, 25th 2008)
#'  The dataset contains the following labels (columns):
#' ***********************************************************************
#'  \itemize{
#'
#' \item [1] = wavelength (nm)
#' \item [2] = refractive index of leaf material (
#'  or the ratio of the velocity of light in a vacuum to its velocity in "leaf
#' medium").
#' \item [3] = specific absorption coefficient of chlorophyll (a+b) (cm2.microg-1)
#' \item [4] = specific absorption coefficient of carotenoids (cm2.microg-1)
#' \item [5] = specific absorption coefficient of brown pigments (arbitrary units)
#' \item [6] = specific absorption coefficient of water (cm-1)
#' \item [7] = specific absorption coefficient of dry matter (g.cm-1)
#' \item [8] = direct light
#' \item [9] = diffuse light
#' \item [10] = dry soil
#' \item [11] = wet soil
#'
#' }
#'
#' @section references:
#' Feret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
#' Properties Model Separating Photosynthetic Pigments, Remote Sensing of
#' Environment
#' 
#' @docType data
#' @keywords datasets
#' @name data_prospect5
#' @usage data(prospect5)
NULL

#' refractive index and specific absorption coefficient for PROSPECT D
#'
#' see http://teledetection.ipgp.jussieu.fr/prosail/ for more details
#' on the data.
#' 
#' @section details:
#' ***********************************************************************
#' data_prospect5 (february, 25th 2008)
#'  The dataset contains the following labels (columns):
#' ***********************************************************************
#'  \itemize{
#'
#' \item [1] = wavelength (nm)
#' \item [2] = refractive index of leaf material (
#'  or the ratio of the velocity of light in a vacuum to its velocity in "leaf
#' medium").
#' \item [3] = specific absorption coefficient of chlorophyll (a+b) (cm2.microg-1)
#' \item [4] = specific absorption coefficient of carotenoids (cm2.microg-1)
#' \item [5] = specific absorption coefficient of brown pigments (arbitrary units)
#' \item [6] = specific absorption coefficient of water (cm-1)
#' \item [7] = specific absorption coefficient of dry matter (g.cm-1)
#' \item [8] = direct light
#' \item [9] = diffuse light
#' \item [10] = dry soil
#' \item [11] = wet soil
#'
#' }
#'
#' @section references:
#' Feret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
#' Properties Model Separating Photosynthetic Pigments, Remote Sensing of
#' Environment
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name data_prospectd
#' @usage data(prospectd)
NULL


#' direct and diffuse light 
#'
#' @section details:
#' ***********************************************************************
#'  \itemize{
#' \item [1] = direct light
#' \item [2] = diffuse light
#' }
#'
#' @section references:
#' Feret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
#' Properties Model Separating Photosynthetic Pigments, Remote Sensing of
#' Environment
#' 
#' @docType data
#' @keywords datasets
#' @name solar
#' @usage data(solar)
NULL


#' soil reflectance
#'
#' @section details:
#' ***********************************************************************
#'  \itemize{
#' \item [1] = wet soil
#' \item [2] = dry soil 
#' }
#'
#' @section references:
#' Feret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
#' Properties Model Separating Photosynthetic Pigments, Remote Sensing of
#' Environment
#' 
#' @docType data
#' @keywords datasets
#' @name soil
#' @usage data(soil)
NULL
