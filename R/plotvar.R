
#' Plots the robust sample variogram
#'
#' Plots the robust sample variogram
#'
#' @param dados Georeferenced data, in the form of a geodata object
#' @param ptos Number of points in the variogram
#'
#' @return Robust variogram
#'
#' @importFrom geoR variog
#'
#' @export
plotvar <- function(dados, ptos) {
    variograma <- variog(dados, uvec = ptos, bin.cloud = T,
                         estimator.type = 'modulus')
    plot(variograma)
    variograma
}