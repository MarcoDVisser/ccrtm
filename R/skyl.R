##' Sky light model
##' 
##' Simple atmospherical model that builds on recommendations
##' from Francois et al. (2002).
##'
##' @param rddt Bi-hemispherical reflectance
##' @param rsdt Directional-hemispherical reflectance for
##'             solar incident flux
##' @param rdot Hemispherical-directional reflectance in viewing direction
##' @param rsot Bi-directional reflectance factor
##' @param Es Solar flux
##' @param Ed Diffuse flux
##' @param tts solar angle
##' @param skyl fraction diffuse
##' @references Francois, C., Ottle, C., Olioso, A., Prevot, L., Bruguier, N., 
##'   Ducros, Y.(2002). Conversion of 400-1100 nm vegetation albedo measurements 
##'   into total shortwave broadband albedo using a canopy radiative transfer model. 
##'   Agronomie 22, 611-618.
##' @examples
##' data(solar)
##' @return a list with hemispherical and directional reflectance.
##' rt<-fRTM(rho~prospect5+foursail)
##' skyl(rt[,"rddt"],rt[,"rsdt"],rt[,"rdot"],rt[,"rsot"],
##' Es=solar[,1],Ed=solar[,2],tts=45,skyl=NULL)
##'
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
    PARdiro <- (1-skyl)*Es
    PARdifo <- (skyl)*Ed  
    
    resh <-  (rddt*PARdifo+rsdt*PARdiro)/(PARdiro+PARdifo) ## hemispherical
    resv <-  (rdot*PARdifo+rsot*PARdiro)/(PARdiro+PARdifo) ## observer
    
    return(list(hemispherical=resh,directional=resv))
}
