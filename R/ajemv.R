
#' Maximum likelihood fit
#'
#' Fits the spatial covariance model through maximum likelihood
#'
#' @param dados *
#' @param est *
#' @param pepita *
#' @param ini.phi *
#' @param ini.sigma2 *
#' @param modelo *
#' @param restrita *
#' @param plot *
#' @param col *
#' @param boxcox *
#' @param ini.lambda *
#' @param tendencia *
#' @param tend *
#'
#' @return The fitted model
#'
#' @importFrom geoR likfit
#'
#' @export
ajemv <- function(dados, est = T, pepita = 0, ini.phi, ini.sigma2,
                  modelo = 'exponential', restrita = F, plot = T, col = 3,
                  boxcox = F, ini.lambda = 0.5, tendencia = F, tend = '1st'){


    if (restrita) {
        if (est)
            emv <- likfit(dados, ini = c(ini.sigma2,ini.phi),
                          nug = pepita, method = 'RML')
        else
            emv <- likfit(dados, ini = c(ini.sigma2,ini.phi),
                          fix.nugget = T, nug = pepita, method = 'RML') }

    if (boxcox) {
        if (est)
            emv <- likfit(dados, ini = c(ini.sigma2,ini.phi), nug = pepita)
        else
            emv <- likfit(dados, ini = c(ini.sigma2,ini.phi),
                          fix.nugget = T, nug = pepita)}

    if (tendencia) {
        if (est)
            emv <- likfit(dados, ini = c(ini.sigma2,ini.phi),
                          nug = pepita, trend = tend)
        else
            emv <- likfit(dados, ini = c(ini.sigma2,ini.phi), fix.nugget = T,
                          nug = pepita, trend = tend)}

    #    if (restrita == T & boxcox == T)
    #    if (restrita == T & tendencia == T)
    #    if (boxcox == T & tendencia == T)

    if (restrita == F & boxcox == F & tendencia == F) {
        if (est)
            emv <- likfit(dados, ini = c(ini.sigma2,ini.phi), nug = pepita)
        else
            emv <- likfit(dados, ini = c(ini.sigma2,ini.phi),
                          fix.nugget = T, nug = pepita)}


    #Plotando - criar variograma tb

    if(plot) lines(emv, lty = 1, col = col)

    emv
}