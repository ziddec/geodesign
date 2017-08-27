
#' Joint Optimal Design
#'
#' Given a set of georeferenced measurements, this function finds 'add.pts' new
#' locations for sampling, calculated jointly using a GA-based heuristic.
#'
#' @param geodata A geodata object containing the initial sample
#' @param add.pts Number of points to be added to the initial sample
#' @param kcontrol Parameters for kriging as a krige.geoR object.
#' See the help for krige.control for further details
#' @param kgrid malha de krigagem. TODO mudar isso! botar
#' pra ser definido por 'n', que nem no SOD
#' @param num.gen Number of generations (iterations)
#' @param gen.size Number of samples in each iteration
#' @param elite.size Number of samples in elite. Should be lesser than or
#' equal to gen.size \%/\% 2 TODO escrever isso direito
#'
#' @return A list with all of the generations and values of the utility function
#' which can be plotted as shown in the examples. TODO RETORNAR SÓ A ÚLTIMA
#'
#' @examples
#' add.pts <- 4 # Number of points to be added
#' n <- 10 # Number of points to define the candidate grid
#' N <- 15 # Number of points for the simulation
#' gen.size <- 10 # Number of samples in each iteration
#' elite.size <- 4 # Number of samples in elite
#' num.gen <- 50 # Number of generations (iterations)
#'
#' library(geoR)
#'
#' set.seed(400)
#' geodata <- grf(N, cov.pars = c(20, 0.45), nugget = 1)
#' beta1 <- mean(geodata$data)
#' emv <- ajemv(geodata,ini.phi = 0.4,ini.sigma2 = 10, pepita = 1,
#'              modelo = 'exponential', plot = F)
#' kcontrol <- krige.control(type.krige = "SK",
#'                           trend.d = "cte",
#'                           nugget = emv$tausq,
#'                           beta = beta1,
#'                           cov.pars = emv$cov.pars)
#' kgrid <- expand.grid(seq(0, 1, l = n), seq(0, 1, l = n))
#'
#' newsample <- JOD(geodata, add.pts, kcontrol, kgrid, num.gen, gen.size,
#' elite.size)
#'
#' utili <- NULL
#' for(i in 1:num.gen) utili[i] <- mean(newsample$utilidades[[i]][1])
#' plot(1:num.gen, utili[1:num.gen], type='l')
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom snow makeCluster
#' @importFrom snow stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom geoR krige.conv
#'
#' @export
JOD <- function(geodata, add.pts, kcontrol, kgrid, num.gen,
                gen.size, elite.size, parallel = T) {

    if (parallel) {
        cl <- makeCluster(c("localhost", "localhost"), "SOCK")
        registerDoSNOW(cl)
    }

    stopifnot(elite.size <= gen.size%/%2)

    N <- length(geodata$data)

    xmin <- min(geodata$coords[,1])
    ymin <- min(geodata$coords[,2])
    xmax <- max(geodata$coords[,1])
    ymax <- max(geodata$coords[,2])

    or.krig <- krige.conv(geodata, locations = kgrid,
                            krige = kcontrol)

    children <- vector("list", gen.size)
    for (i in 1:gen.size) {
        pts.x <- runif(add.pts, min = xmin, max = xmax)
        pts.y <- runif(add.pts, min = ymin, max = ymax)
        new.pt.coord <- matrix(c(pts.x, pts.y), nrow = add.pts)
        sample.data <- append(geodata$data, numeric(add.pts))
        names(new.pt.coord) <- c("x", "y")
        sample.coord <- rbind(geodata$coords, new.pt.coord)
        rownames(sample.coord) <- NULL
        rownames(new.pt.coord) <- NULL
        children[[i]] <- as.geodata(cbind(sample.coord, sample.data))
    }

    #lista.amostras é como se fosse o 'children' da 1a geração
    #Pelo amor de deus não esquecer de mudar o '0' ali (numeric(add.pts)) na hr de
    #fazer a outra fç utilidade!!

    children.it <- util.it <- vector('list', num.gen)

    for (kk in 1:num.gen){
        krig.list <- vector("list", gen.size)
        predvar.util <- foreach(i = 1:gen.size, .packages = 'geoR', .combine = "c") %dopar% {
            krig.list[[i]] <- krige.conv(children[[i]], locations = kgrid,
                                               krige = kcontrol)
            mean((or.krig$krige.var -
                      krig.list[[i]]$krige.var)/or.krig$krige.var)
        }

        util.order <- order(-predvar.util)
        ranked.list <- list()
        for (i in 1:gen.size) {
            ranked.list[[i]] <- children[[util.order[i]]]
        }
        util.it[[kk]] <- predvar.util[util.order]

        ##### P5, P6, P7.1 e P7.2 #####

        children <- vector('list', gen.size)

        for (i in 1:elite.size) {
            children[[i]] <- ranked.list[[i]]

            chosen <- sample((N+1):(N+add.pts), add.pts%/%2 + if (add.pts%%2 != 0) 1 else 0)
            chosen.pts <- children[[i]]$coords[chosen,]

            pts.x <- runif(add.pts%/%2, min = xmin, max = xmax) #2 pontos sao aleatorios
            pts.y <- runif(add.pts%/%2, min = ymin, max = ymax)
            random.pts <- cbind(pts.x, pts.y)

            names(chosen.pts) <- names(random.pts) <- c("x", "y")
            pts <- rbind(chosen.pts, random.pts)

            children.data <- c(geodata$data, numeric(add.pts))
            children.coord <- rbind(geodata$coords, pts)
            rownames(children.coord) <- rownames(children.data) <- NULL
            children[[elite.size+i]] <- as.geodata(cbind(children.coord, children.data))
        }

        for (i in (2*elite.size):gen.size) {
            pts.x <- runif(add.pts, min = xmin, max = xmax)
            pts.y <- runif(add.pts, min = ymin, max = ymax)
            random.pts <- cbind(pts.x, pts.y)
            names(random.pts) <- c("x", "y")

            children.data <- c(geodata$data, numeric(add.pts))
            children.coord <- rbind(geodata$coords, random.pts)
            rownames(children.coord) <- rownames(children.data) <- NULL
            children[[i]] <- as.geodata(cbind(children.coord, children.data))
        }

        children.it[[kk]] <- children
    }

    if(parallel) stopCluster(cl)

    list(children = children.it, utilidades = util.it)
}
