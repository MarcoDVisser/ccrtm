% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lim.R
\name{lidf}
\alias{lidf}
\title{Leaf inclination distribution function models
s3 method for calling leaf models.}
\usage{
lidf(pars)
}
\arguments{
\item{pars}{a parameter vector c(angles, LIDFa, LIDFb) with a class
lidf.\link{modelnumber}. Models include:
\itemize{
\item \link{1} = Dlagden distribution (1, lidf.1)
\item \link{2} = Ellipsoid (Campebll) distribution (2, lidf.2)
\item \link{3} = Beta distribution (3, lidf.3)
\item \link{4} = One parameter beta distribution (4, lidf.4)
}
Models 1 and 2 are the standard models from the SAIL model.
Two parameter models use parameters LIDFa and LIDFb, while
single parameter models use only LIDFa (ignoring any supplied LIDFb).

More information on the Dlagden and Ellipsoid parameter is
given in Verhoef, W. (1998),theory of radiative transfer models applied
in optical remote sensing of vegetation canopies (PhD thesis).

The beta distribution is the typical beta distribution as often implemented
(as in dbeta(x,LIDFa, LIDFb)). Where x is a value between 0 and 90,
that gives the angular density over 0 and 90 degrees (rescaled to 0 and 1).

The one parameter beta distribution is given by LIDFa*x^(LIDFa-1).
Where x is a value between 0 and 90, that given the angular density
over 0 and 90 degrees (rescaled to 0 and 1).}
}
\value{
a vector of proportions for each leaf angle calculated
from each leaf inclination model.
}
\description{
Leaf inclination distribution function models
s3 method for calling leaf models.
}
