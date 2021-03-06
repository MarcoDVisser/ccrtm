% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prospect.R
\name{prospectd}
\alias{prospectd}
\title{PROSPECT model version D}
\usage{
prospectd(param)
}
\arguments{
\item{param}{A named vector of PROSPECT parameters (note: program ignores case):
 \itemize{
\item [1] = leaf structure parameter (N)
\item [2] = chlorophyll a+b content in ug/cm2 (Cab)
\item [3] = carotenoids content in ug/cm2 (Car)
\item [4] = Leaf anthocyanin content (ug/cm2) (Canth)
\item [5] = brown pigments content in arbitrary units (Cbrown)
\item [6] = equivalent water thickness in g/cm2  (Cw)
\item [7] = leaf dry matter content in g/cm2 - lma - (Cm)
}}
}
\value{
spectra matrix with leaf reflectance and transmission 
for wavelengths 400 to 2500nm:
 \itemize{
\item [1] = leaf reflectance (rho)
\item [2] = leaf transmission (tau)
}
}
\description{
The PROSPECTD leaf reflectance model. The model was implemented based on
Jacquemoud and Ustin (2019), and is further described in detail in Feret et al (2017).
PROSPECT models use the plate models developed in 
Allen (1969) and Stokes (1862).
}
\references{
Jacquemoud, S., and Ustin, S. (2019). Leaf optical properties. 
  Cambridge University Press.

Feret, J.B., Gitelson, A.A., Noble, S.D., Jacquemoud, S. (2017).
  PROSPECT-D: Towards modeling leaf optical properties through a complete lifecycle. 
  Remote Sens. Environ. 193, 204-215.

Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R.
  (1969), Interaction of isotropic ligth with a compact plant leaf,
  Journal of the Optical Society of American, 59:1376-1379.

Stokes G.G. (1862), On the intensity of the light
  reflected from or transmitted through a pile of plates,
  Proceedings of the Royal Society of London, 11:545-556.
}
