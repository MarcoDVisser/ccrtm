#' S3- methods for Generate defaults settings and parameters
#' for all supported models. 
#'
#' @param model a ccrtm formula or character vector of modelnames
#' @param \\dots not used.
#' (e.g. "prospect5")
#' @return a data.frame with default model parameters
#' @export
getDefaults<-function(model=NULL, ...){

    rtmModels <- getModels()

    if(class(model)=="formula"){

        ## get order of RTM
        reqMods <- attr(terms(model),"term.labels")

        ## assert that repeated terms are not lost!
        test <-unlist(strsplit(as.character(model)[[3]]," \\+ "))
        if(length(test)>length(reqMods)) reqMods <- test
        

        ## check models
        if(any(!reqMods%in%rtmModels$model)){
            
            stop(reqMods[!reqMods%in%rtmModels$model],
                 " not implemented")
        }
    } else if(class(model)=="character"){

        test <- any(sapply(model,function(X) !X%in%rtmModels$model))
        if(test) stop("model not recognized")
        
        reqMods <- model
        }
        
    def <- vector("list",length(reqMods))
    
    for(i in seq_len(length(reqMods))){
        
        mod <- reqMods[i]
        class(mod) <- reqMods[i]
        def[[i]] <- defaults(mod)
    }
    
    names(def) <- reqMods
    if(length(def)==1) return(def[[1]])
    
    return(def)     
}

## defualts function
defaults <- function(x, ...){
    UseMethod("defaults",x)
}

## get some default values for the prospect model 
##      - N     = leaf structure parameter
##      - Cab   = chlorophyll a+b content 
##      - Car   = carotenoids content 
##      - Cbrown= brown pigments content in arbitrary units
##      - Cw    = equivalent water thickness
##      - Cm    = dry matter content 
defaults.prospect5<-function(x){
   
    ## typical values are the rice values
    typical = c(2.698,70.8,20, 0.0117,0.009327,0.001)
    ## N, Cab,Car, Cw, Cm, Cbrown, resid
    omega <- c(5,250,80,.4,0.3,10) 
    alpha <- c(1,10,1,1e-5,1e-5,0) 
    
    names(typical)<-c("N", "Cab","Car", "Cw", "Cm", "Cbrown")
    names(omega)<-c("N", "Cab","Car", "Cw", "Cm", "Cbrown")
    names(alpha)<-c("N", "Cab","Car", "Cw", "Cm", "Cbrown")
    
    def = data.frame(best = typical)
    def$lower = alpha
    def$upper = omega
    
    return(def)
}

## get some default values for the prospect model 
##       - N     = leaf structure parameter
##       - Cab   = chlorophyll a+b content 
##       - Car   = carotenoids content 
##       - Canth = Leaf anthocyanin content
##       - Cbrown= brown pigments content in arbitrary units
##       - Cw    = equivalent water thickness
##       - Cm    = dry matter content 
##
defaults.prospectd<-function(x){
    
    
    ## typical values are the rice values
    typical = c(2.698,70.8,20, 0.0117,0.009327,5,0.001)
    ## N, Cab,Car, Cw, Cm, Cbrown, resid
    omega <- c(5,250,80,.4,0.3,50,10) 
    alpha <- c(1,10,1,1e-5,1e-5,0,0) 
    
    names(typical)<-c("N", "Cab","Car", "Cw", "Cm","Canth", "Cbrown")
    names(omega)<-c("N", "Cab","Car", "Cw", "Cm", "Canth","Cbrown")
    names(alpha)<-c("N", "Cab","Car", "Cw", "Cm","Canth","Cbrown")
    
    def = data.frame(best = typical)
    def$lower = alpha
    def$upper = omega
    
     return(def)
}




## Generate default values for the foursail model
##       - psoil = linear mixing weight for 2 soil
##                 reflectance types
##       - LAI   =  Leaf area index
##       - TypeLidf= Type of leaf angle distribution. 
##       - lidfa = Leaf angle distribution parameter a. 
##       - lidfb = Leaf angle distribution parameter b.
##       - hspot = Hotspot parameter
##       - tts   = Solar zenith angle (degrees)
##       - tto   = Observer zenith angle (degrees)
##       - psi   = Relative azimuth angle (degrees)
defaults.foursail<-function(x){

    nms <- c("psoil","LAI","TypeLidf","lidfa","lidfb",
             "hspot","tts","tto","psi")
#             "sigma")
    
    typical = c(0,4.77,1,-0.35,-0.15,
                0.01,30,10,0)
#                0.1)
    
    omega <- c(1,50,1,1,1,
               1,90,90,180)
#               4) 
    alpha <- c(0,0.1,1,-1,-1,
               1e-6,0,0,0)
#               1e-8) 
    
    names(typical)<-nms
    names(omega)<-nms
    names(alpha)<-nms
    
    def = data.frame(best = typical)
    def$lower = alpha
    def$upper = omega
    
    return(def)

}   

## Generate default values for the foursail2 model 
##       - psoil = psoil: Dry/Wet soil factor
##       - LAI   =  Leaf area index
##       - TypeLidf= Type of leaf angle distribution. 
##       - lidfa = Leaf angle distribution parameter a. 
##       - lidfb = Leaf angle distribution parameter b.
##       - Cv    = vertical crown coverage fraction.
##      - Zeta0 = tree shape factor (D/H).
##      - fb    = fraction secondary particles (fraction "brown").
##       - diss  = Canopy dissociation factor for secondary particles.
##       - skyl  = linear mixing weight for diffuse & direct radiation.
##                 (ratio of diffuse to total incident radiation)
##       - hspot = Hotspot parameter.
##       - tts   = Solar zenith angle (degrees).
##       - tto   = Observer zenith angle (degrees).
##       - psi   = Relative azimuth angle (degrees).
defaults.foursail2<-function(x){

    nms <- c("psoil","LAI","TypeLidf","lidfa","lidfb",
             "hspot",
             "Cv","Zeta","fb","diss",
             "tts","tto","psi")
#             "sigma")
    
    typical  <-  c(0,8.77,1,-0.35,-0.15,
                   0.01,                   
                   0.95,1/3,0.5,0.5,
                   30,10,0)
#                0.1)
    
    omega <- c(1,50,1,1,1,
               1,
               1,1,1,1,
               90,90,180)
#               4) 
    alpha <- c(0,0.1,1,-1,-1,
               1e-6,
               0,0,0,0,
               0,0,0)
#               1e-8) 
    
    names(typical)<-nms
    names(omega)<-nms
    names(alpha)<-nms
    
    def = data.frame(best = typical)
    def$lower = alpha
    def$upper = omega

    return(def)
    
}


## Generate default values for the foursail2b model 
##       - psoil = psoil: Dry/Wet soil factor
##       - LAI   =  Leaf area index
##       - TypeLidf= Type of leaf angle distribution. 
##       - lidfa = Leaf angle distribution parameter a. 
##       - lidfb = Leaf angle distribution parameter b.
##      - Cv    = vertical crown coverage fraction.
##      - Zeta0 = tree shape factor (D/H).
##      - fb    = fraction secondary particles (fraction "brown").
##       - diss  = Canopy dissociation factor for secondary particles.
##       - skyl  = linear mixing weight for diffuse & direct radiation.
##                 (ratio of diffuse to total incident radiation)
##       - hspot = Hotspot parameter.
##       - tts   = Solar zenith angle (degrees).
##       - tto   = Observer zenith angle (degrees).
##       - psi   = Relative azimuth angle (degrees).
defaults.foursail2b<-function(x){

    nms <- c("psoil","LAI","TypeLidf","lidfa","lidfb",
             "hspot",
             "Cv","Zeta","fb","diss",
             "tts","tto","psi")
#             "sigma")
    
    typical  <-  c(0,8.77,2,35,45,
                   0.01,                   
                   0.95,1/3,0.5,0.5,
                   30,10,0)
#                0.1)
    
    omega <- c(1,50,2,0,0,
               1,
               1,1,1,1,
               90,90,180)
#               4) 
    alpha <- c(0,0.1,2,90,90,
               1e-6,
               0,0,0,0,
               0,0,0)
#               1e-8) 
    
    names(typical)<-nms
    names(omega)<-nms
    names(alpha)<-nms
    
    def = data.frame(best = typical)
    def$lower = alpha
    def$upper = omega

    return(def)
    
}

## get some default values for the prospect model 
##       - N     = leaf structure parameter
##       - Cab   = chlorophyll a+b content 
##       - Car   = carotenoids content 
##       - Cbrown= brown pigments content in arbitrary units
##       - Cw    = equivalent water thickness
##       - Cm    = dry matter content 
##
defaults.skyl<-function(x){
   
    ## typical values are the rice values
    typical = c(0.1)
    omega <- c(1) 
    alpha <- c(0) 
    
    names(alpha)<-names(omega)<-names(typical)<-c("skyl")
    
    def = data.frame(best = typical)
    def$lower = alpha
    def$upper = omega
    
    return(def)
}
