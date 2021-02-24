#' Forest Light Interaction Model (FLIM)
#'
#' The FLIM model was first described by Rosema et al (1992).
#' In FLIM forests are assumed a discontinous mix of tree crowns
#' and gaps. Reflectance is modelled in terms of the
#' probabilty to observe either a gap (background)
#' or a tree crown. Both gaps and crowns may be shaded. 
#'
#' @param Rc Canopy reflectance at infinite depth
#' @param Rg soil/background reflectance
#' @param To transmission in viewing direction
#' @param Ts transmission in sun direction
#' @param area area of stand
#' @param params a named vector of parameters: 
#'  \itemize{
#' \item [1] = D, stand density (confounded with cd)
#' \item [2] = cd, crown diameter (confounded with D)
#' \item [3] = h, mean crown height
#' \item [6] = Solar zenith angle (tts)
#' \item [7] = Observer zenith angle (tto)
#' \item [8] = Sun-sensor azimuth angle (psi)
#' }
#'
#' @details
#'
#' Confounded parameters pairs cannot be inversely
#' estimated, one of the pairs should be set to 1.
#'
#' @return a list with reflectance, and the fractions of shaded 
#' and sunexplosed crowns, shaded and sun exposed open space. 
#'
#' @examlpes
#' parvec<- c(alpha = 0.5,lai=5,cd=15,h=30,d=10,tto=10,tts=20,psi=30)
#' flim(0.1,0.1,params=parvec)
#'    
#' @references Rosema, A., Verhoef, W., Noorbergen, H., Borgesius, J.J. (1992). 
#' A new forest light interaction model in support of forest monitoring. 
#' Remote Sens. Environ. 42, 23-41. 
#' @importFrom pracma sec
#' 
flim <- function(Rc,Rg,To=NULL,Ts=NULL,params,area=10000){

    ## degrees to radians
    rd <- pi/180

    test1 <- all(c(is.null(To),is.null(Ts)))
    test2 <- all(c("alpha","lai")%in%names(params))

    if(test1) {
        if(test2){
        alpha <- params["alpha"] # extinction coef
        lai <- params["lai"] # leaf area index
        } else {
            stop("Provide either To,Ts or alpha and LAI")
        }
    } 
    
    cd <- params["cd"] # mean crown diameter
    h <- params["h"] # mean crown height
    D <- params["d"] # stand density
    tto <- params["tto"]*rd  # observer zenith
    tts <- params["tts"]*rd # sun zenith
    psi <- params["psi"]*rd # relative azimuth
    

    ## Some house keeping
    k <- (pi*(cd/2)^2)/area
    g <- sqrt(tan(tto)^2 + tan(tts)^ 2 - 2*tan(tto)*tan(tts)*cos(psi)) ## geometric factor (eqn 4)

    ## Coverage and Shadowing 
    gmexp <- pracma::sec(tto*rd)/pracma::sec(tts) # geometric exponent
    co <- 1-exp(-k*D*pracma::sec(tto)) # eq 1 from Rosema et al 1992
    cs <- 1-(1-co)^gmexp #eqn 3

    p <- exp(-g*h/cd)

    ## ground surface fractions
    open <- 1-co
    sunlit <- 1-cs
    corf <- p*sqrt(co*open*cs*sunlit) # correlation factor
    
    Fcd <- co*cs + corf # eqn 6
    Fcs <- co*sunlit - corf # eqn 7
    Fod <- open*cs - corf # eqn 8
    Fos <- open*sunlit + corf #eqn 9

    if(test1){
        ## transmittance in sun and observer direction
        To <- exp(-alpha*lai*pracma::sec(tto)) # eqn 10
        Ts <- To^gmexp #eqn 12
    }
    
    G <- Fcd*Ts*To+Fcs*To+Fod*Ts+Fos # (eqn 13)
    C <- (1-Ts*To)*cs*co  # eqn 16
    R <- Rc*C + Rg*G # eqn 15

    names(R)<-"rho"
    return(list(R,crowndark=Fcd,crownsun=Fcs,opendark=Fod,opensun=Fos,open=open,sunlit=sunlit))
}

    
