% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sod.R
\name{SOD}
\alias{SOD}
\title{Sequential Optimal Design}
\usage{
SOD(geodata, add.pts, n = ceiling(sqrt(10 * nrow(geodata$coords))), util,
  kcontrol, parallel = T, qtl = 0.75, p = 0.5, shape = NULL)
}
\arguments{
\item{geodata}{A geodata object containing the initial sample}

\item{add.pts}{Number of points to be added to the initial sample}

\item{n}{Number of points in the sides of the candidate grid. If a vector
of length 1, both sides of the grid will have the same number of points. If
a vector of length 2, the first value indicates the number of points along
the x axis, and the second value indicates the number of points along the y
axis. Defaults to the square root of ten times the number of observations in
the geodata object}

\item{util}{Utility function to be used. Possibilities are `predvar`,
`extrprob` or `mixed`. See the Details section for further details}

\item{kcontrol}{Parameters for kriging as a krige.geoR object. See the help
for krige.control for further details}

\item{parallel}{Indicates if the code should run in parallel. Defaults to
TRUE}

\item{qtl}{Reference quantile for the extreme probability utility function.
Defaults to 0.75}

\item{p}{Weight of the predictive variance function for the calculation of
the mixed utility function. Defaults to 0.5}

\item{shape}{A SpatialPointsDataFrame object, read with function readOGR()
from rgdal, which contains the area of interest. The candidate grid will be
located inside the limits indicated by this shapefile. Defaults to NULL, and
in this case, uses the bounding box of the observed data to create the
candidate grid}
}
\value{
Coordinates of the new sampling locations
}
\description{
Given a set of georeferenced measurements, this function finds the `add.pts`
optimal locations for sampling.
}
\details{
The value `predvar` for util refers to the reduction in predictive variance
utility function. `extrprob` refers to the utility function that favours the
locations that will have a higher probability of observing extreme values.
`mixed` refers to a mixed utility function which uses weight p for the
reduction in predictive variance utility function and weight 1-p for the
extreme observations utility function.
}
\examples{
library(geoR)

add.pts <- 10 # Number of points to be added
n <- 10 # Number of points to define the candidate grid
N <- 15 # Number of points for the simulation (only for testing)
qtl <- 0.75

# Model parameters: 20, 0.45, 1
set.seed(300)
simdata <- grf(N, cov.pars = c(20, 0.45), nug  =  1)

# Visualization of simulated data:
# points(simdata, cex.min = 1, cex.max = 2.5, col = gray(seq(1, 0, l = 4)))

beta1 <- mean(simdata$data)
m1 <- as.matrix(simdata$coords)
emv <- ajemv(simdata, ini.phi = 0.4, plot = F,
             ini.sigma2 = 10, pepita = 1, modelo = 'exponential')

new.pts <- SOD(simdata, add.pts, n, util = 'predvar',
               kcontrol = krige.control(type.krige = "SK",
                                        trend.d = "cte",
                                        nugget = emv$tausq,
                                        beta = mean(simdata$data),
                                        cov.pars = emv$cov.pars),
               parallel = F)
new.pts.bayes <- SOD(simdata, add.pts, n, util = 'predvar',
                     kcontrol = prior.control(beta.prior = "flat",
                                              sigmasq.prior = "reciprocal",
                                              phi.prior="uniform",
                                              phi.discrete=seq(0,2,l=20),
                                              tausq.rel.prior = "uniform",
                                              tausq.rel.discrete = seq(0, 1, l=20)),
                     parallel = F)

# Old points and new points
par(mfrow = c(1, 2))
plot(simdata$coords, pch = 16, main = 'Classical kriging')
points(new.pts, pch = '+', col = 'red')
plot(simdata$coords, pch = 16, main = 'Bayesian kriging')
points(new.pts.bayes, pch = '+', col = 'red')
par(mfrow = c(1, 1))

}
