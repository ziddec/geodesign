
#' Sequential Optimal Design
#'
#' Given a set of georeferenced measurements, this function finds the `add.pts`
#' optimal locations for sampling.
#'
#' @param geodata A geodata object containing the initial sample
#' @param add.pts Number of points to be added to the initial sample
#' @param n Number of points in the sides of the candidate grid. If a vector
#' of length 1, both sides of the grid will have the same number of points. If
#' a vector of length 2, the first value indicates the number of points along
#' the x axis, and the second value indicates the number of points along the y
#' axis. Defaults to the square root of ten times the number of observations in
#' the geodata object
#' @param util Utility function to be used. Possibilities are `predvar`,
#' `extrprob` or `mixed`. See the Details section for further details
#' @param kcontrol Parameters for kriging as a krige.geoR object. See the help
#' for krige.control for further details
#' @param parallel Indicates if the code should run in parallel. Defaults to
#' TRUE
#' @param qtl Reference quantile for the extreme probability utility function.
#' Defaults to 0.75
#' @param p Weight of the predictive variance function for the calculation of
#' the mixed utility function. Defaults to 0.5
#' @param shape A SpatialPointsDataFrame object, read with function readOGR()
#' from rgdal, which contains the area of interest. The candidate grid will be
#' located inside the limits indicated by this shapefile. Defaults to NULL, and
#' in this case, uses the bounding box of the observed data to create the
#' candidate grid
#'
#' @return Coordinates of the new sampling locations
#'
#' @details
#' The value `predvar` for util refers to the reduction in predictive variance
#' utility function. `extrprob` refers to the utility function that favours the
#' locations that will have a higher probability of observing extreme values.
#' `mixed` refers to a mixed utility function which uses weight p for the
#' reduction in predictive variance utility function and weight 1-p for the
#' extreme observations utility function.
#'
#' @examples
#' library(geoR)
#'
#' add.pts <- 10 # Number of points to be added
#' n <- 10 # Number of points to define the candidate grid
#' N <- 15 # Number of points for the simulation (only for testing)
#' qtl <- 0.75
#'
#' # Model parameters: 20, 0.45, 1
#' set.seed(300)
#' simdata <- grf(N, cov.pars = c(20, 0.45), nug  =  1)
#'
#' # Visualization of simulated data:
#' # points(simdata, cex.min = 1, cex.max = 2.5, col = gray(seq(1, 0, l = 4)))
#'
#' beta1 <- mean(simdata$data)
#' m1 <- as.matrix(simdata$coords)
#' emv <- ajemv(simdata, ini.phi = 0.4, plot = F,
#'              ini.sigma2 = 10, pepita = 1, modelo = 'exponential')
#'
#' new.pts <- SOD(simdata, add.pts, n, util = 'predvar',
#'                kcontrol = krige.control(type.krige = "SK",
#'                                         trend.d = "cte",
#'                                         nugget = emv$tausq,
#'                                         beta = mean(simdata$data),
#'                                         cov.pars = emv$cov.pars),
#'                parallel = F)
#'
#' # Old points and new points
#' plot(simdata$coords, pch = 16)
#' points(new.pts, pch = '+', col = 'red')
#'
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom geoR krige.conv
#' @importFrom geoR as.geodata
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom sp bbox
#' @importFrom sp SpatialPoints
#' @importFrom sp proj4string
#' @importFrom spatstat as.owin
#' @importFrom spatstat gridcentres
#'
#'
#' @export
SOD <- function(geodata, add.pts, n = ceiling(sqrt(10*nrow(geodata$coords))),
                util, kcontrol, parallel = T, qtl = 0.75, p = 0.5, shape = NULL) {

    if (!("geodata" %in% class(geodata)))
        stop("Expected first argument to be a geodata")

    if (mode(add.pts) != "numeric" || add.pts <= 0)
        stop("Invalid value for `add.pts`")

    if (mode(n) != "numeric")
        stop("Expected `n` to be numeric")

    if (length(n) == 1) {
        nx <- ny <- n
    } else if (length(n) == 2) {
        nx <- n[1]
        ny <- n[2]
    } else {
        stop("Invalid length for n")
    }

    if (!(util %in% c('predvar', 'extrprob', 'mixed')))
        stop("Invalid value for util")

    if (!("krige.geoR" %in% class(kcontrol)))
        stop("Expected `kcontrol` to be a `krige.geoR` object")

    if (typeof(parallel) != "logical")
        stop("Expected `parallel` to be a logical value")

    if (mode(qtl) != "numeric" || (qtl < 0 | qtl > 1))
        stop("Expected `qtl` to be a value between 0 and 1")

    if (mode(p) != "numeric" || (p < 0 | p > 1))
        stop("Expected `p` to be a value between 0 and 1")

    if (parallel) {
        cl <- makeCluster(c("localhost", "localhost"), "SOCK")
        registerDoSNOW(cl)
    }

    if (!is.null(shape)) {
        if (typeof(shape) != "S4") #?
            stop("Expected `shape` to be a SpatialPolygonsDataFrame value")
    }

    if (!is.null(shape)) {
        # 'shape' should be a SpatialPolygonsDataFrame object
        shape.window <- as.owin(shape) # Owin
        initgrid <- gridcentres(shape.window, nx, ny) # List
        initgridP <- SpatialPoints(initgrid) # SpatialPoints
        proj4string(initgridP) <- proj4string(shape)
        finalgrid <- initgridP[shape,] # SpatialPoints
        cgrid <- kgrid <- finalgrid@coords # matrix
    } else {
        box <- bbox(geodata$coords)
        kgrid <- cgrid <- expand.grid(seq(box[1,1], box[1,2], l = nx),
                             seq(box[2,1], box[2,2], l = ny))
    }



    #botar um if aqui pra criar malha preditiva, dependendo de se há shapefile
    #botar mensagens para quando está criando malha preditiva, pq demora
    #botar 2 mensagens: criando malha, otimizando

    cat("Optimizing...\n")



    if (util == 'predvar') {

        it.predvar.util <- list() # predictive variance utility f.
        best.pt <- NULL

        for (g in 1:(add.pts)) {
            cat(sprintf("Iteration %d out of %d%s\n",
                        g, add.pts, if (g == add.pts) "! Yay!" else ""))

            # "Original" kriging (to be compared with krigings at each point)
            capture.output(
                or.krig <- krige.conv(geodata, locations = kgrid,
                                     krige = kcontrol)
            )

            # Checks for common points between existing coords and 'cgrid'
            m1 <- as.matrix(geodata$coords)
            m2 <- as.matrix(cgrid)
            ptscommon <- NULL
            cont <-  0
            cont2 <-  0
            ptnr <- NULL
            for (i in 1:nrow(m2)) {
                cont2 <- cont2 + 1
                for (j in 1:nrow(m1)) {
                    if (all(m1[j,] == m2[i,])) {
                        cont <- cont + 1 # how many common points exist
                        ptscommon <- rbind(ptscommon,m1[j,]) # coordinates of common points
                        # Pode apenas fazer m2[ptnr, ] no final para obter mesma matriz
                        ptnr <- c(ptnr,cont2) # (sequence) number of the common points
                    }
                }
            }

            # Calculates decrease in predictive variance ("predvar")

            names(cgrid) <- names(as.data.frame(geodata$coords))
            geodata.list <- list()
            values <- c(geodata$data,0)
            df.list <- list()
            krig.list <- list()
            i <- 0
            predvar.util <- NULL

            if (cont == 0) { #if there are no points in common
                predvar.util <- foreach(i = 1:length(cgrid[,1]), .packages = 'geoR', .combine = "c") %dopar% {
                    df.list[[i]] <- rbind(geodata$coords,cgrid[i,])
                    geodata.list[[i]] <- as.geodata(data.frame(cbind(df.list[[i]],values)))
                    capture.output(
                        krig.list[[i]] <- krige.conv(geodata.list[[i]], locations = kgrid,
                                                     krige = kcontrol)
                    )
                    predvar.util <- mean(or.krig$krige.var - krig.list[[i]]$krige.var)
                    df.list[[i]] <- 0
                    krig.list[[i]] <- 0
                    geodata.list[[i]] <- 0
                    predvar.util

                }
            }

            if (cont > 0) { #if there are points in common
                predvar.util <- foreach(i = 1:length(cgrid[,1]), .packages = 'geoR', .combine = "c") %dopar% {
                    if (!(i %in% ptnr)) { # if point 'i' is NOT one of the common points
                        df.list[[i]] <- rbind(geodata$coords,cgrid[i,])
                        geodata.list[[i]] <- as.geodata(data.frame(cbind(df.list[[i]],values)))
                        capture.output(
                            krig.list[[i]] <- krige.conv(geodata.list[[i]], locations = kgrid,
                                                         krige = kcontrol)
                        )
                        predvar.util <- mean(or.krig$krige.var - krig.list[[i]]$krige.var)
                        df.list[[i]] <- 0
                        krig.list[[i]] <- 0
                        geodata.list[[i]] <- 0
                        predvar.util
                    }
                }
            }


            # Erases utility function value if location is already in sample

            if (cont > 0) {
                for (i in 1:(cont)) {
                    predvar.util <- insert(predvar.util, 0, ptnr[i])
                }
            }


            # Linear transformation in utility function ( R+ -> (0,1) )

            coef1 <- range(predvar.util)[1]
            coef2 <- range(predvar.util)[2]
            co <- matrix(c(coef1, coef2, 1, 1), 2)
            ld <- c(0, 1)
            sol <- solve(co, ld)
            a <- sol[1]
            b <- sol[2]
            tr.predvar.util <- a * predvar.util + b # transformed utility function

            # Saves corrected utility function values for iteration "g"
            it.predvar.util[[g]] <- tr.predvar.util

            # Optimal sampling location for utility function #1
            best.pt <- c(best.pt, which(tr.predvar.util == max(tr.predvar.util)))

            # Point coordinates and value
            best.coord <- cgrid[best.pt,]
            best.value <- or.krig$predict[best.pt]

            # Merges data with optimal point
            geodata$data <- c(geodata$data, best.value[g])
            geodata$coords <- rbind(geodata$coords, best.coord[g,])
            rownames(geodata$coords) <- NULL
            geodata <- as.geodata(cbind(geodata$coords, geodata$data))
        }
    } # if util == 'predvar' ends here

    if (util == 'extrprob') {

        it.extrprob.util <- list() # extreme probabilities utility f
        best.pt <- NULL

        for (g in 1:(add.pts)) {
            cat(sprintf("Iteration %d out of %d%s\n",
                        g, add.pts, if (g == add.pts) "! Yay!" else ""))

            # "Original" kriging (to be compared with krigings at each point)

            capture.output(
                or.krig <- krige.conv(geodata, locations = kgrid,
                                     krige = kcontrol)
            )

            # Checks for common points between 'coords' and 'cgrid'
            m1 <- as.matrix(geodata$coords)
            m2 <- as.matrix(cgrid)
            ptscommon <- NULL
            cont <-  0
            cont2 <-  0
            ptnr <- NULL

            for (i in 1:length(m2[,1])) {
                cont2 = cont2 + 1
                for (j in 1:length(m1[,1])) {
                    if (all(m1[j,] == m2[i,])) {
                        cont <- cont + 1 # how many common points exist
                        ptscommon <- rbind(ptscommon,m1[j,]) # coordinates of common points
                        ptnr <- c(ptnr,cont2) # number of the points that are in common
                    }
                }
            }
            # Utility function #2:
            # Increase in the probability of observing extreme events ("extrprob")
            # (reference quantile: "qtl")

            extrprob.util <- NULL
            tolerance <- quantile(geodata$data, qtl)

            for (i in 1:(nx*ny)) {
                extrprob.util[i] <- (1 - pnorm(tolerance, sd = sqrt(or.krig$krige.var[i]),
                                               mean = or.krig$predict[i]))^2
                # "probability of value in point i being greater than the tolerance"
            }

            # Erases utility function values if location is already in sample
            if (cont != 0) {
                extrprob.util[ptnr] <- 0
            }

            # Saves corrected utility function values for iteration "g"
            it.extrprob.util[[g]] <- extrprob.util

            # Optimal sampling location for utility function #2
            best.pt <- c(best.pt, which(extrprob.util == max(extrprob.util)))

            # Point coordinates and value
            best.coord <- cgrid[best.pt,]
            best.value <- or.krig$predict[best.pt]
            colnames(best.coord) <- c("x", "y")

            # Merges geodata with optimal point
            geodata$data <- c(geodata$data,best.value[g])
            geodata$coords <- rbind(geodata$coords,best.coord[g,])
            rownames(geodata$coords) <- NULL
            geodata <- as.geodata(cbind(geodata$coords,geodata$data))
        }
    } # if util == 'extrprob' ends here

    if (util == 'mixed') {
        it.extrprob.util <- list()
        it.predvar.util <- list()
        it.mixed.util <- list() # mixed utility f.
        best.pt <- NULL

        for (g in 1:(add.pts)) {
            cat(sprintf("Iteration %d out of %d%s\n",
                        g, add.pts, if (g == add.pts) "! Yay!" else ""))

            # "Original" kriging (to be compared with krigings at each point)
            capture.output(
                or.krig <- krige.conv(geodata, locations = kgrid,
                                     krige = kcontrol)
            )

            #       Checks for common points between existing coords and 'cgrid'
            m1 <- as.matrix(geodata$coords)
            m2 <- as.matrix(cgrid)
            ptscommon <- NULL
            cont <-  0
            cont2 <-  0
            ptnr <- NULL
            for (i in 1:length(m2[,1])) {
                cont2 = cont2 + 1
                for (j in 1:length(m1[,1])) {
                    if (all(m1[j,] == m2[i,])) {
                        cont <- cont + 1 # how many common points exist
                        ptscommon <- rbind(ptscommon,m1[j,]) # coordinates of common points
                        ptnr <- c(ptnr,cont2) # (sequence) number of the common points
                    }
                }
            }

            #       Calculates decrease in predictive variance ("predvar")

            names(cgrid) <- names(as.data.frame(geodata$coords))
            geodata.list <- list()
            values <- c(geodata$data,0)
            df.list <- list()
            krig.list <- list()
            i <- 0
            predvar.util <- NULL

            if (cont == 0) { #if there are no points in common

                predvar.util <- foreach(i = 1:length(cgrid[,1]), .packages = 'geoR', .combine = "c") %dopar% {
                    df.list[[i]] <- rbind(geodata$coords,cgrid[i,])
                    geodata.list[[i]] <- as.geodata(data.frame(cbind(df.list[[i]],values)))
                    capture.output(
                        krig.list[[i]] <- krige.conv(geodata.list[[i]], locations = kgrid,
                                                     krige = kcontrol)
                    )
                    predvar.util <- mean(or.krig$krige.var - krig.list[[i]]$krige.var)
                    df.list[[i]] <- 0
                    krig.list[[i]] <- 0
                    geodata.list[[i]] <- 0
                    predvar.util

                }
            }

            if (cont > 0) { #if there are points in common
                predvar.util <- foreach(i = 1:length(cgrid[,1]), .packages = 'geoR', .combine = "c") %dopar% {
                    if (!(i %in% ptnr)) { # if point 'i' is NOT one of the common points
                        df.list[[i]] <- rbind(geodata$coords, cgrid[i,])
                        geodata.list[[i]] <- as.geodata(data.frame(cbind(df.list[[i]], values)))
                        capture.output(
                            krig.list[[i]] <- krige.conv(geodata.list[[i]], locations = kgrid,
                                                         krige = kcontrol)
                        )
                        predvar.util <- mean(or.krig$krige.var - krig.list[[i]]$krige.var)
                        df.list[[i]] <- 0
                        krig.list[[i]] <- 0
                        geodata.list[[i]] <- 0
                        predvar.util
                    }
                }
            }


            # Erases utility function value if location is already in sample

            if (cont > 0) {
                for (i in 1:(cont)) {
                    predvar.util <- insert(predvar.util, 0, ptnr[i])
                }
            }


            # Linear transformation in utility function ( R+ -> (0,1) )

            coef1 <- range(predvar.util)[1]
            coef2 <- range(predvar.util)[2]
            co <- matrix(c(coef1, coef2, 1, 1), 2)
            ld <- c(0, 1)
            sol <- solve(co,ld)
            a <- sol[1]
            b <- sol[2]
            tr.predvar.util <- a * predvar.util + b # transformed utility function

            # Utility function #2:
            # Increase in the probability of observing extreme events ("extrprob")
            # (reference quantile: "qtl")

            extrprob.util <- NULL
            tolerance <- quantile(geodata$data, qtl)

            for (i in 1:(nx * ny)) {
                extrprob.util[i] <- (1 - pnorm(tolerance, sd = sqrt(or.krig$krige.var[i]),
                                               mean = or.krig$predict[i]))^2
                # "probability of value in point i being greater than the tolerance"
            }

            # Erases utility function values if location is already in sample
            if (cont != 0) {
                extrprob.util[ptnr] <- 0
            }

            # Utility function #3:
            # Mixed utility function
            # p (weight of predvar.util) defaults to 0.5

            mixed.util <- p * tr.predvar.util + (1 - p) * extrprob.util
            it.mixed.util[[g]] <- mixed.util

            # Saves (corrected) utility function values for iteration "g"
            it.predvar.util[[g]] <- tr.predvar.util
            it.extrprob.util[[g]] <- extrprob.util
            it.mixed.util[[g]] <- extrprob.util

            # Optimal sampling location for utility function #3
            best.pt <- c(best.pt, which(mixed.util == max(mixed.util)))

            # Point coordinates and value
            best.coord <- cgrid[best.pt,]
            best.value <- or.krig$predict[best.pt]

            # Merges data with optimal point
            geodata$data <- c(geodata$data, best.value[g])
            geodata$coords <- rbind(geodata$coords, best.coord[g,])
            rownames(geodata$coords) <- NULL
            geodata <- as.geodata(cbind(geodata$coords, geodata$data))
        }

    } # if util == 'mixed' ends here

    util.evolution <- NULL # Utility function evolution through iterations
    # (plot images or make line plot with mean at each iteration)
    if (util == 'predvar') util.evolution <- it.predvar.util else
        if (util == 'extrprob') util.evolution <- it.extrprob.util else
            util.evolution <- it.mixed.util

    nc <- nrow(geodata$coords)
    geodata$coords[(nc - add.pts + 1):nc,]
}
