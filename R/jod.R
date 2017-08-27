#' Joint Optimal Design
#'
#' @param dados A geodata object containing the initial sample
#' @param k número de pontos que quero adicionar à minha amostra
#' @param kcontrol controle pro kriging
#' @param malha malha de krigagem, é o kgrid
#' @param J n de gerações
#' @param G tamanho de cada geração
#' @param m n de amostras na elite. tem q ser <= G div 2 e >0 (botar p todos)
#'
#' @examples
#' k <- 4 # número de pontos que quero adicionar à minha amostra
#' NN <- 10 # n de pontos em um lado da malha de krigagem
#' N <- 15 # n de dados iniciais - para a simulação
#' G <- 10 # n de elementos em cada geração do AG
#' m <- 4 # n de amostras na elite
#' J <- 500 # n de gerações
#'
#' library(geoR)
#'
#' set.seed(400)
#' dados <- grf(N, cov.pars = c(20, 0.45), nugget = 1)
#' beta1 <- mean(dados$data)
#' emv <- ajemv(dados,ini.phi = 0.4,ini.sigma2 = 10, pepita = 1,
#'              modelo = 'exponential', plot = F)
#' kcontrol <- krige.control(type.krige = "SK",
#'                           trend.d = "cte",
#'                           nugget = emv$tausq,
#'                           beta = beta1,
#'                           cov.pars = emv$cov.pars)
#' malha <- expand.grid(seq(0, 1, l = NN), seq(0, 1, l = NN))
#'
#' novosfilhos <- JOD(dados, 5, kcontrol, malha, 500, 10, 4)
#'
#' utili <- NULL
#' for(i in 1:J) utili[i] = mean(novosfilhos$utilidades[[i]][1])
#' plot(1:J, utili[1:J], type='l')
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom snow makeCluster
#' @importFrom snow stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom geoR krige.conv
#'
#' @export
JOD <- function(dados, k, kcontrol, malha, J, G, m, parallel = T) {

    if (parallel) {
        cl <- makeCluster(c("localhost", "localhost"), "SOCK")
        registerDoSNOW(cl)
    }

    stopifnot(m <= G%/%2)

    N <- length(dados$data)

    xmin <- min(dados$coords[,1])
    ymin <- min(dados$coords[,2])
    xmax <- max(dados$coords[,1])
    ymax <- max(dados$coords[,2])

    krigeagem <- krige.conv(dados, locations = malha,
                            krige = kcontrol)

    filhotes <- vector("list", G)
    for (i in 1:G) {
        ptos.x <- runif(k, min = xmin, max = xmax)
        ptos.y <- runif(k, min = ymin, max = ymax)
        coords.ptos.novos <- matrix(c(ptos.x, ptos.y), nrow = k)
        amostra.dados <- append(dados$data, numeric(k))
        names(coords.ptos.novos) <- c("x", "y")
        amostra.coords <- rbind(dados$coords, coords.ptos.novos)
        rownames(amostra.coords) <- NULL
        rownames(coords.ptos.novos) <- NULL
        filhotes[[i]] <- as.geodata(cbind(amostra.coords, amostra.dados))
    }

    #lista.amostras é como se fosse o 'filhotes' da 1a geração
    #Pelo amor de deus não esquecer de mudar o '0' ali (numeric(k)) na hr de
    #fazer a outra fç utilidade!!

    filhotes.it <- utilidades.it <- vector('list', J)

    for (kk in 1:J){
        lista.krigagens <- vector("list", G)
        u.varpred <- foreach(i = 1:G, .packages = 'geoR', .combine = "c") %dopar% {
            lista.krigagens[[i]] <- krige.conv(filhotes[[i]], locations = malha,
                                               krige = kcontrol)
            mean((krigeagem$krige.var -
                      lista.krigagens[[i]]$krige.var)/krigeagem$krige.var)
        }

        ordem <- order(-u.varpred)
        lista.rankeada <- list()
        for (i in 1:G) {
            lista.rankeada[[i]] <- filhotes[[ordem[i]]]
        }
        utilidades.it[[kk]] <- u.varpred[ordem]

        ##### P5, P6, P7.1 e P7.2 #####

        filhotes <- vector('list', G)

        for (i in 1:m) {
            filhotes[[i]] <- lista.rankeada[[i]]

            escolhidos <- sample((N+1):(N+k), k%/%2 + if (k%%2 != 0) 1 else 0)
            ptos.escolhidos <- filhotes[[i]]$coords[escolhidos,]

            ptos.x <- runif(k%/%2, min = xmin, max = xmax) #2 pontos sao aleatorios
            ptos.y <- runif(k%/%2, min = ymin, max = ymax)
            ptos.aleatorios <- cbind(ptos.x, ptos.y)

            names(ptos.escolhidos) <- names(ptos.aleatorios) <- c("x", "y")
            ptos <- rbind(ptos.escolhidos, ptos.aleatorios)

            filhote.dados <- c(dados$data, numeric(k))
            filhote.coords <- rbind(dados$coords, ptos)
            rownames(filhote.coords) <- rownames(filhote.dados) <- NULL
            filhotes[[m+i]] <- as.geodata(cbind(filhote.coords, filhote.dados))
        }

        for (i in (2*m):G) {
            ptos.x <- runif(k, min = xmin, max = xmax)
            ptos.y <- runif(k, min = ymin, max = ymax)
            ptos.aleatorios <- cbind(ptos.x,ptos.y)
            names(ptos.aleatorios) <- c("x", "y")

            filhote.dados <- c(dados$data, numeric(k))
            filhote.coords <- rbind(dados$coords, ptos.aleatorios)
            rownames(filhote.coords) <- rownames(filhote.dados) <- NULL
            filhotes[[i]] <- as.geodata(cbind(filhote.coords, filhote.dados))
        }

        filhotes.it[[kk]] <- filhotes
    }

    if(parallel) stopCluster(cl)

    list(filhotes = filhotes.it, utilidades = utilidades.it)
}
