
#' Minimum squares fit
#'
#' Fits the spatial covariance model through minimum squares
#'
#' @param dados *
#' @param est *
#' @param pepita *
#' @param ini.phi *
#' @param ini.sigma2 *
#' @param modelo *
#' @param pesos *
#' @param ptos *
#' @param plot *
#' @param col *
#'
#' @return The fitted model
#'
#' @importFrom geoR variog
#' @importFrom geoR variofit
#'
#' @export
ajmqp <- function(dados, est =  T, pepita = 0, ini.phi,
                  ini.sigma2, modelo = 'exponential',
                  pesos = 'npairs', ptos = 15, plot = T, col = 2) {
    #Estimando

    variograma.robusto <- variog(dados, uvec = ptos, bin.cloud = T,
                                 estimator.type ='modulus')

    if (est) {
        mqp <- variofit(variograma.robusto, ini = c(ini.sigma2,ini.phi),
                        cov.model = modelo, nugget = pepita, weights = pesos)}
    else {
        mqp <- variofit(variograma.robusto, ini = c(ini.sigma2, ini.phi),
                        cov.model = modelo, fix.nugget = T, nugget = pepita) }

    # Vendo se o resultado faz sentido -> isso ? s? para o exponencial e
    # gaussiano n?? E os outros?

    m1 <- dados$coords

    dmax <- max(dist(m1))

    if (modelo == 'spherical') {
        if (dmax < mqp$cov.pars[2]) {
            message('O par?metro phi sugere que h? correla??o a uma dist?ncia maior que dmax')
        }
    }

    if (modelo == 'gaussian') {
        if (dmax < 3*(mqp$cov.pars[2])) {
            message('O par?metro phi sugere que h? correla??o a uma dist?ncia maior que dmax')
        }
    }

    #Plotando

    if (plot) {
        plot(variograma.robusto)
        lines(mqp, lty = 1, col = col)
    }

    mqp
}
