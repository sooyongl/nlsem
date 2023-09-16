#' @description
#' This function comes from `nlsem` package
# Calculate weights and node points for mixture functions via Gauss-Hermite
# quadrature as defined in Klein & Moosbrugger (2000)
quadrature <- function (m, k) 
{
  one.dim <- gaussquad::hermite.h.quadrature.rules(m)[[m]]
  test <- as.matrix(expand.grid(lapply(vector("list", k), function(x) {
    x <- 1:m
    x
  })))
  final.nodes <- matrix(one.dim$x[test], ncol = k, byrow = FALSE)
  permute.weights <- matrix(one.dim$w[test], ncol = k, byrow = FALSE)
  final.weights <- apply(permute.weights, 1, prod)
  n <- final.nodes * sqrt(2)
  w <- final.weights * pi^(-k/2)
  out <- list(n = n, w = w, k = k, m = m)
  out
}

#' @description
#' This function mimics `nlsem` package
#' 
# Fill a model created with specify_sem with parameters given as a vector;
# mostly needed to simulate data from a prespecified model; exported function
fill_model <- function(model, parameters) {
  
  matrices = model$matrices
  
  # for (class in names(matrices)) {
  
  matrices[[1]][["Lambda.x"]][is.na(matrices[[1]][["Lambda.x"]])] <- parameters[1:2]
  matrices[[1]][["Gamma"]][is.na(matrices[[1]][["Gamma"]])] <- parameters[3:5]
  matrices[[1]][["Theta.d"]][is.na(matrices[[1]][["Theta.d"]])] <- parameters[6:8]
  
  matrices[[1]][["Psi"]][is.na(matrices[[1]][["Psi"]])] <- parameters[9]
  
  matrices[[1]][["Phi"]][1,] <- parameters[10:12]
  matrices[[1]][["Phi"]][,1] <- parameters[10:12]
  matrices[[1]][["nu.x"]][2:3] <- parameters[13:14]
  matrices[[1]][["alpha"]][1] <- parameters[15]
  matrices[[1]][["tau"]][1] <- parameters[16]
  matrices[[1]][["Omega"]][1,2] <- parameters[17]
  
  out <- list(matrices=matrices, info=model$info)
  
  out
}

#' #' @description
#' #' This function comes from `nlsem` package
#' #' 
#' # Fill a model created with specify_sem with parameters given as a vector;
#' # mostly needed to simulate data from a prespecified model; exported function
#' fill_model <- function(model, parameters) {
#'   
#'   # stopifnot(inherits(model, "singleClass") || inherits(model, "semm")
#'   #           || inherits(model, "nsemm"))
#'   # 
#'   # stopifnot(count_free_parameters(model) == length(parameters))
#'   
#'   matrices <- model$matrices
#'   
#'   # if (inherits(model, "singleClass")) {
#'   constraints <- "direct1"
#'   # } else {
#'   # constraints <- model$info$constraints
#'   # }
#'   
#'   switch(EXPR = constraints,
#'          
#'          indirect = {
#'            
#'            parnames <- c("Phi", "tau")
#'            ind.list <- list()
#'            num.list <- list()
#'            for (class in names(matrices)) {
#'              for (parname in parnames) {
#'                ind.list[[class]][[parname]] <- is.na(matrices[[class]][[parname]])
#'                num.list[[class]][[parname]] <-
#'                  length(matrices[[class]][[parname]][is.na(matrices[[class]][[parname]])])
#'              }
#'            }
#'            
#'            for (i in seq_along(matrices$class1)) {
#'              matrix.i <- matrices$class1[[i]]
#'              # number of NA's in matrix
#'              num.na <- length(matrix.i[is.na(matrix.i)])
#'              if (num.na > 0) {
#'                matrix.i[is.na(matrix.i)] <- parameters[1:num.na]
#'                parameters <- parameters[-(1:num.na)]
#'                for (class in names(matrices)) {
#'                  matrices[[class]][[i]] <- matrix.i
#'                }
#'              }
#'            }
#'            for (class in names(matrices)[-1]) {
#'              for (parname in parnames) {
#'                ind <- ind.list[[class]][[parname]]
#'                num <- num.list[[class]][[parname]]
#'                if (num > 0) {
#'                  matrices[[class]][[parname]][ind] <- parameters[1:num]
#'                  parameters <- parameters[-(1:num)]
#'                }
#'              }
#'            }
#'          },
#'          
#'          direct1 = {
#'            
#'            for (class in names(matrices)) {
#'              for (i in seq_along(matrices[[class]])) {
#'                matrix.i <- matrices[[class]][[i]]
#'                # number of NA's in matrix
#'                num.na <- length(matrix.i[is.na(matrix.i)])
#'                if (num.na > 0) {
#'                  matrix.i[is.na(matrix.i)] <- parameters[1:num.na]
#'                  parameters <- parameters[-(1:num.na)]
#'                  matrices[[class]][[i]] <- matrix.i
#'                }
#'              }
#'            }
#'            
#'          },
#'          
#'          direct2 = {
#'            
#'            if (inherits(model, "semm")) {
#'              parnames <- c("Gamma", "Beta", "Psi", "Phi", "alpha", "tau")
#'            } else {
#'              parnames <- c("Gamma", "Beta", "Psi", "Phi", "alpha", "tau", "Omega")
#'            }
#'            ind.list <- list()
#'            num.list <- list()
#'            for (class in names(matrices)) {
#'              for (parname in parnames) {
#'                ind.list[[class]][[parname]] <- is.na(matrices[[class]][[parname]])
#'                num.list[[class]][[parname]] <-
#'                  length(matrices[[class]][[parname]][is.na(matrices[[class]][[parname]])])
#'              }
#'            }
#'            
#'            for (i in seq_along(matrices$class1)) {
#'              matrix.i <- matrices$class1[[i]]
#'              # number of NA's in matrix
#'              num.na <- length(matrix.i[is.na(matrix.i)])
#'              if (num.na > 0) {
#'                matrix.i[is.na(matrix.i)] <- parameters[1:num.na]
#'                parameters <- parameters[-(1:num.na)]
#'                for (class in names(matrices)) {
#'                  matrices[[class]][[i]] <- matrix.i
#'                }
#'              }
#'            }
#'            for (class in names(matrices)[-1]){
#'              for (parname in parnames) {
#'                ind <- ind.list[[class]][[parname]]
#'                num <- num.list[[class]][[parname]]
#'                if (num > 0) {
#'                  matrices[[class]][[parname]][ind] <- parameters[1:num]
#'                  parameters <- parameters[-(1:num)]
#'                }
#'              }
#'            }
#'          }
#'   )     # end of switch
#'   
#'   # fill symmetric matrices
#'   for (class in names(matrices)) {
#'     tryCatch({matrices[[class]]$Phi <-
#'       fill_symmetric(matrices[[class]]$Phi)}, error=function(e) e,
#'       warning=function(w) w)
#'     tryCatch({matrices[[class]]$Psi <-
#'       fill_symmetric(matrices[[class]]$Psi)}, error=function(e) e,
#'       warning=function(w) w)
#'   }
#'   
#'   out <- list(matrices=matrices, info=model$info)
#'   class(out) <- class(model)
#'   out
#' }