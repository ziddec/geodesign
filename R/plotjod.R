#' @export
plot.geodesign.JOD <- function(obj) {
    par(mfrow = c(2,1))
    plot(obj$init.pts, pch = 16, main = 'Final design')
    points(obj$best.sample[(nrow(obj$init.pts)+1):nrow(obj$best.sample),], pch = '+', col = 'red')
    plot(obj$util.evolution, type = 'l', col = 'blue', main = 'Utility evolution',
         xlab = 'Number of iterations', ylab = 'Increase in utility function')
    par(mfrow = c(1,1))
}