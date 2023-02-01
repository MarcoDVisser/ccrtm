##' Sky light model
##'
##' Simple skyl atmospheric model.
##'
##' The version implemented here can also include a dependence of the sun zenith angle
##' after Danner et al. (2019) who build on recommendations from Francois et al. (2002).
##'
##' @param rddt Bi-hemispherical reflectance
##' @param rsdt Directional-hemispherical reflectance for
##'             solar incident flux
##' @param rdot Hemispherical-directional reflectance in viewing direction
##' @param rsot Bi-directional reflectance factor
##' @param Es Solar flux
##' @param Ed Diffuse flux
##' @param tts solar angle
##' @param skyl diffuse fraction, if NULL skyl is estimated using the
##' tts (solar angle).
##'
##' @references Francois, C., Ottle, C., Olioso, A., Prevot, L., Bruguier, N.,
##'   Ducros, Y.(2002). Conversion of 400-1100 nm vegetation albedo measurements
##'   into total shortwave broadband albedo using a canopy radiative transfer model.
##'   Agronomie 22, 611-618.
##' @references Danner M, Berger K, Wocher M, Mauser W, Hank T.
##'   Fitted PROSAIL Parameterization of Leaf Inclinations, Water Content
##'   and Brown Pigment Content for Winter Wheat and Maize Canopies.
##'   Remote Sensing. 2019; 11(10):1150.
##' @examples
##'
##' data(solar)
##' rt<-fRTM(rho~prospect5+foursail)
##' skyl(rt[,"rddt"],rt[,"rsdt"],rt[,"rdot"],rt[,"rsot"],
##' Es=solar[,1],Ed=solar[,2],tts=45,skyl=NULL)
##'
##' @return a list with hemispherical and directional reflectance.
##' @export
skyl <- function(rddt,rsdt,rdot,rsot,Es,Ed,tts,skyl=NULL){

  if(is.null(skyl)){
    rd <- pi/180
    ## skyl following the recommendations from
    ## Francois et al. (2002) Conversion of 400-1000 nm vegetation albedo
    ## measurements into total shortwave broadband albedo using a canopy
    ## radiative transfer model, Agronomie
    skyl   <-   0.847- 1.61*sin((90-tts)*rd)+
      1.04*sin((90-tts)*rd)*sin((90-tts)*rd)
  }

  ## Direct / diffuse light
  Rdir <- (1-skyl)*Es
  Rdif <- (skyl)*Ed

  resh <-  (rsdt*Rdir+rddt*Rdif)/(Rdir+Rdif) ## hemispherical
  resv <-  (rsot*Rdir+rdot*Rdif)/(Rdir+Rdif) ## observer

  return(list(hemispherical=resh,directional=resv))
}
