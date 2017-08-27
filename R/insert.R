
#' Inserts an element in a vector
#'
#' Inserts an element `el` into the vector `vec` in the position `pos`
#'
#' @param vec Vector of interest
#' @param el Element to be inserted
#' @param pos Position to be inserted. Must be lesser or equal to the length of
#' `vec`
#'
#' @return A new vector based on `vec` with the new element
#'
#' @export
insert <- function(vec, el, pos){
    n <- length(vec)

    stopifnot(pos <= n)

    if (pos == 1) {
        c(el, vec[1:n])
    } else {
        c(vec[1:(pos-1)], el, vec[pos:n])
    }
}
