#' ccrtm: Coupled Chain Radiative Transfer Models.
#'
#' A collection of radiative transfer models that can form
#' a coupled chain to model radiative transfer across multiple
#' spatial scales from leaf to canopy to stand.
#'
#' Currently implemented models that can be coupled:
#'
#'  \itemize{
#'
#' \item [1] = PROSPECT 5, 5B and D
#' \item [2] = FOURSAIL, and FOURSAIL2
#' \item [3] = FLIM
#'
#' }
#'
#' Currently being tested / or to be implemented models
#'
#'  \itemize{
#'
#' \item [1] = LIBERTY, PROCOSINE
#' \item [2] = INFORM*
#'
#' }
#'
#' *available as lower-level library (see ccrtm github page).
#'
#' @details
#' ## Generating predictions
#' To generate model prediction the typical approach is to use
#' the fRTM function and apply a formula that specifies the
#' both the expected output (left hand) and the different models
#' you would like to couple to generate the output (right hand).
#'
#' At current the following radiative transfer models,
#' and corresponding formula, are given in the next table
#'
#' | Example of a formula                     | Model     |
#' | :--------------------------------------  | :-------: |
#' | rho ~ prospect5                          | prospect5 |
#' | rho ~ prospectd                          | prospectd |
#' | rho ~ prospect5 + foursail               | PROSAIL   |
#' | rho ~ prospect5 + foursail               | PROSAIL   |
#' | rho ~ prospectd + foursail2              | PROSAIL   |
#' | rho ~ prospectd + prospect5 + foursail2  | PROSAIL2  |
#' | rho ~ prospectd + foursail2b             | PROSAIL2b |
#' | rho ~ prospectd + foursail + flim + sky  | INFORM*   |
#' | rho ~ prospectd + foursail + flim        | INFORM*   |
#'
#' *INFORM is currently restricted to a lower-level function
#' only. See the cctrm github readme page on how to use it.
#'
#' In the above examples, prospectd can be replaced with prospect5
#' - or vice versa - if so desired.
#'
#' Also, note that the formula "rho ~ prospectd + foursail2" is the
#' same as "rho ~ prospectd + prospectd + foursail2" and
#' both expect a names list of 3 parameter vectors.
#'
#' See the help files for details on each right hand
#' component. For instance, ?foursail provides more elaboration
#' on the 4SAIL model and gives an example for lower-level
#' implementations of each component model. See also ?flim,
#' ?foursail2, ?foursail2b, ?prospect5, and ?prospectd.
#'
#' Tranmission can also be returned if specified in
#' the left-hand component of the formula:
#'
#' | Formula                             | Model     |
#' | :---------------------------------  | :-------: |
#' | rho + tau ~ prospect5               | prospect5 |
#' | rho + tau ~ prospectd + foursail    | 4SAIL     |
#'
#' The examples above indicate that the users wishes to predict
#' transmission next to reflectance. More specifically, The first
#' returns leaf reflectance and transmission while
#' the second returns 4 components of canopy reflectance and
#' canopy transmission in the solar and viewing direction.
#'
#' ## PROSAIL and diffuse and direct light
#' In contrast to some FORTRAN and MATLAB implementation,
#' the sky light model is not implemented by default in
#' ccrtm. This is because it is not a standard component of 4SAIL.
#' In addition, this would affect limit the application of other
#' more realistic atmospheric models. You can apply it by using ?skyl
#' on model output obtained from fRTM (see example in ?skyl).
#' The sky light model _is_ implemented for the INFORM model as per March 2022.
#'
#' ## Parameters
#' Parameters are given as input to fRTM as a named list. See ?getDefaults for
#' examples on how to structure parameters. Individual models can
#' consulted on each parameter (e.g. see ?foursail2b or ?prospect5).
#'
#' ## Leaf inclination models.
#' Canopy models such as foursail use leaf inclination models.
#' In cctrm four inclination models are implemented.
#' see ?lidf for more details.
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

#' refractive index and specific absorption coefficients for
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

#' refractive index and specific absorption coefficients for PROSPECT D
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
