#' @export
print.geodesign.JOD <- function(obj) {
    cat('Joint Optimal Design\n')
    cat('--------------------\n')
    cat(sprintf('Initial dataset cointaining %d points\n', nrow(obj$init.pts)))
    cat(sprintf('%d points added considering %d iterations\n', obj$add.pts, obj$num.it))
    cat(sprintf('Increase in utility function: %.1f%%\n', 100*obj$util.evolution[obj$num.it]))
}