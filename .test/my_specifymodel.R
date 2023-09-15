specify_sem

function (num.x, num.y, num.xi, num.eta, xi, eta, constraints = c("indirect", 
                                                                  "direct1", "direct2"), num.classes = 1, rel.lat = "default", 
          interaction = "none") 
{
  
  num.x = 5
  num.y = 1
  num.xi = 3
  num.eta = 1
  xi = "x1-x3,x4,x5"
  eta = "y1"
  num.classes = 1
  rel.lat = "default"
  constraints = 'indirect'
  interaction = "xi1:xi2"
  
  
  
  # if (!is.numeric(num.x) || !is.numeric(num.y) || !is.numeric(num.xi) || 
  #     !is.numeric(num.eta) || !is.numeric(num.classes)) {
  #   stop("Number of variables or classes must be numeric.")
  # } else if (num.x < num.xi || num.y < num.eta) {
  #   stop("The model contains not enough observed variables.")
  # }
  # stopifnot(num.x > 0, num.y >= 0, num.xi > 0, num.eta >= 0, num.classes > 0)
  
  
  # if (interaction != "none") {
  #   interact.matrix <- calc_interaction_matrix(unlist(strsplit(interaction, 
  #                                                              ",")))
  #   if (max(interact.matrix) > num.xi) {
  #     stop("Interaction effects contain more xi's than defined.")
  #   }
  #   if (num.eta > 1) {
  #     if (!grepl("eta", interaction)) {
  #       stop("For more than one eta, specification for interaction must be something like eta1~xi1:xi2.")
  #     }
  #   }
  # }
  
  model.class <- nlsem:::get_model_class(num.classes, interaction)
  xi.s <- unlist(strsplit(xi, ","))
  if (length(xi.s) != num.xi) {
    stop("Number of xi's and assignation of x's to xi's does not match. See ?specify_sem.")
  }
  xi.ind <- list()
  for (i in seq_len(num.xi)) {
    xi.ind[[i]] <- nlsem:::grep_ind(xi.s[i])
    if (max(xi.ind[[i]]) > num.x) {
      stop("Number of x's assinged to xi exceeds x's specified. See ?specify_sem.")
    }
  }
  eta.s <- unlist(strsplit(eta, ","))
  if (length(eta.s) != num.eta) {
    stop("Number of eta's and assignation of y's to eta's does not match. See ?specify_sem.")
  }
  eta.ind <- list()
  for (i in seq_len(num.eta)) {
    eta.ind[[i]] <- nlsem:::grep_ind(eta.s[i])
    if (max(eta.ind[[i]]) > num.y) {
      stop("Number of y's assinged to eta exceeds y's specified. See ?specify_sem.")
    }
  }
  Lambda.x <- matrix(0, nrow = num.x, ncol = num.xi)
  for (i in seq_len(num.xi)) {
    Lambda.x[xi.ind[[i]], i] <- c(1, rep(NA, length(xi.ind[[i]]) - 
                                           1))
  }
  Lambda.y <- matrix(0, nrow = num.y, ncol = num.eta)
  for (i in seq_len(num.eta)) {
    Lambda.y[eta.ind[[i]], i] <- c(1, rep(NA, length(eta.ind[[i]]) - 
                                            1))
  }
  if (rel.lat == "default") {
    Gamma <- matrix(nrow = num.eta, ncol = num.xi)
    Beta <- diag(num.eta)
  } else {
    GB <- rel_lat(rel.lat, num.eta = num.eta, num.xi = num.xi)
    Gamma <- tryCatch({
      GB[[grep("G", names(GB))]]
    }, error = function(e) matrix(nrow = num.eta, ncol = num.xi))
    Beta <- tryCatch({
      GB[[grep("B", names(GB))]]
    }, error = function(e) diag(num.eta))
  }
  Theta.d <- diag(NA, nrow = num.x)
  Theta.e <- diag(NA, nrow = num.y)
  
  Psi <- matrix(NA, nrow = num.eta, ncol = num.eta)
  Psi[upper.tri(Psi)] <- 0
  
  Phi <- matrix(NA, nrow = num.xi, num.xi)
  Phi[upper.tri(Phi)] <- 0
  
  nu.x <- matrix(NA, nrow = num.x, ncol = 1)
  for (i in seq_len(num.xi)) nu.x[xi.ind[[i]][1]] <- 0
  
  nu.y <- matrix(NA, nrow = num.y, ncol = 1)
  for (i in seq_len(num.eta)) nu.y[eta.ind[[i]][1]] <- 0
  
  alpha <- matrix(NA, nrow = num.eta, ncol = 1)
  tau <- matrix(NA, nrow = num.xi, ncol = 1)
  
  Omega <- matrix(0, nrow = num.xi, ncol = num.xi)
  if (interaction != "none") {
    interaction.s <- unlist(strsplit(interaction, ","))
    eta.logical <- matrix(nrow = num.eta, ncol = length(interaction.s))
    for (i in seq_len(num.eta)) {
      eta.logical[i, ] <- grepl(paste0("eta[", i, "]"), 
                                interaction.s)
    }
    which.eta <- apply(eta.logical, 1, which)
    if (nrow(eta.logical) == 1) {
      ind <- nlsem:::calc_interaction_matrix(interaction.s)
      Omega[ind] <- NA
      nlsem:::test_omega(Omega)
    } else {
      Omega <- array(0, dim = c(num.xi, num.xi, num.eta))
      for (i in seq_len(num.eta)) {
        eta.row <- which(eta.logical[i, ])
        ind <- nlsem:::calc_interaction_matrix(interaction.s[eta.row])
        Omega[, , i][ind] <- NA
        nlsem:::test_omega(Omega[, , i])
      }
    }
  }
  
  matrices <- list()
  for (c in seq_len(num.classes)) {
    if (model.class == "singleClass") {
      matrices[[c]] <- list(Lambda.x = Lambda.x, Lambda.y = Lambda.y, 
                            Gamma = Gamma, Beta = Beta, Theta.d = Theta.d, 
                            Theta.e = Theta.e, Psi = Psi, Phi = Phi, nu.x = nu.x, 
                            nu.y = nu.y, alpha = alpha, tau = tau, Omega = Omega)
    }
    # else if (model.class == "semm") {
    #   matrices[[c]] <- list(Lambda.x = Lambda.x, Lambda.y = Lambda.y, 
    #                         Gamma = Gamma, Beta = Beta, Theta.d = Theta.d, 
    #                         Theta.e = Theta.e, Psi = Psi, Phi = Phi, nu.x = nu.x, 
    #                         nu.y = nu.y, alpha = alpha, tau = tau)
    # }
    # else {
    #   matrices[[c]] <- list(Lambda.x = Lambda.x, Lambda.y = Lambda.y, 
    #                         Gamma = Gamma, Beta = Beta, Theta.d = Theta.d, 
    #                         Theta.e = Theta.e, Psi = Psi, Phi = Phi, nu.x = nu.x, 
    #                         nu.y = nu.y, alpha = alpha, tau = tau, Omega = Omega)
    # }
  }
  
  names(matrices) <- paste0("class", seq_len(num.classes))
  w <- matrix(1/num.classes, nrow = num.classes, ncol = 1)
  # constraints <- match.arg(constraints)
  model <- list(matrices = matrices, info = list(num.xi = num.xi, 
                                                 num.eta = num.eta, num.x = num.x, num.y = num.y, constraints = constraints, 
                                                 num.classes = num.classes, par.names = list(), w = w))
  class(model) <- model.class
  model$info$par.names <- nlsem:::get_parnames(model = model, constraints = constraints)
  model$info$bounds <- nlsem:::bounds(model)
  model
}

