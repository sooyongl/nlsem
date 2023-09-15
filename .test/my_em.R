# Compute EM algorithm for LMS

# Argument inputs ---------------------------------------------
model = my_model
data = data

start =  pars.start
qml = F
verbose = F
convergence = 0.1
max.iter = 200
max.mstep = 1
max.singleClass = 1 
neg.hessian = TRUE
m = 16
optimizer = c("nlminb")

# EM Estimation begins ----------------------------------------------------
## Validate starting values (e.g., negative value for variances)  
par.new <- nlsem:::convert_parameters_singleClass(model, start)

ll.ret <- NULL # Save likelihoods across iterations (from M step)
ll.new <- 0    # Likelihood every iteration (from M step)

num.iter <- 0
run <- TRUE

# Computing E and M steps until run == FALSE ------------------------------
while (run) {
  
  # Check likelihood trends ------
  if (num.iter > 3) {
    if (ll.new - ll.old > 0) {
      warning("Loglikelihood should be increasing.")
    }
  }
  
  if (TRUE) {
    cat(paste("Iteration", num.iter + 1, "\n"))
    cat("Doing expectation-step \n")
  }
  
  ll.old <- ll.new
  par.old <- par.new # par.new will keep updated
  
  # names(model$matrices$class1)[grep("Phi", names(model$matrices$class1))] <- "A"
  #------------------------------------------------------------------------
  # Compute E-step --------------------------------------------------------
  #------------------------------------------------------------------------
  
  model = model
  parameters = par.old
  dat = data
  m = m
  
  ##### Compute E-step         <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  stopifnot(count_free_parameters(model) == length(parameters))
  mod.filled <- fill_model(model = model, parameters = parameters)
  mod.filled$matrices
  
  k <- 1 # nlsem:::get_k(mod.filled$matrices$class1$Omega)
  # Calculate Hermite Quadrature <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  if (k != 0) {
    one.dim <- gaussquad::hermite.h.quadrature.rules(m)[[m]]
    test <- as.matrix(expand.grid(
      lapply(vector("list", k), function(x) {x <- 1:m; x })))
    final.nodes <- matrix(one.dim$x[test], ncol = k, byrow = FALSE)
    permute.weights <- matrix(one.dim$w[test], ncol = k, byrow = FALSE)
    final.weights <- apply(permute.weights, 1, prod)
    n <- final.nodes * sqrt(2)
    w <- final.weights * pi^(-k/2)
    out <- list(n = n, w = w, k = k, m = m)
    out
    
    quad <- out
    
    V <- quad$n
    w <- quad$w
  } 
  # Calculate Hermite Quadrature <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  P <- NULL
  for (i in seq_along(w)) {
    # i = 1
    
    p.ij <- w[i] * mvtnorm::dmvnorm(dat, 
                                    # mean = nlsem:::mu_lms(mod.filled, V[i, ]), 
                                    # sigma = nlsem:::sigma_lms(mod.filled, V[i, ])
                                    mean = MU_lms(mod.filled, V[i, ]),
                                    sigma = SIGMA_lms(mod.filled, V[i, ])
    )
    P <- cbind(P, p.ij, deparse.level = 0)
  }
  P <- P/rowSums(P)
  P
  # Compute E-step         <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  #------------------------------------------------------------------------
  # Compute M-setp --------------------------------------------------------
  #------------------------------------------------------------------------
  # m.step <- nlsem:::mstep_lms(
  model = model#, 
  # ls.str(model)
  # model$info
  P = P#, 
  dat = data#, 
  parameters = par.old#, 
  m = m#, 
  optimizer = optimizer#,
  max.mstep = max.mstep#)
  ## nlsem:::mstep_lms <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  # optimizer <- match.arg(optimizer)
  # if (optimizer == "nlminb") {
  #   if (is.null(control$iter.max)) {
  #     control$iter.max <- max.mstep
  #   }
  #   else warning("iter.max is set for nlminb. max.mstep will be ignored.")
  #   suppress_NaN_warnings(
  #     est <- nlminb(start = parameters,
  #                   objective = loglikelihood_lms, dat = dat, model = model,
  #                   P = P, 
  #                   upper = model$info$bounds$upper, 
  #                   lower = model$info$bounds$lower,
  #                   control = control, ...))
  # }
  # else {
  # if (is.null(control$maxit)) {
  # control <- list()
  # control$maxit <- max.mstep
  # } 
  #   else warning("maxit is set for optim. max.mstep will be ignored.")
  
  control <- list()
  control$maxit <- max.mstep
  
  ## Likelihood for LMS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  if(F) {
    # nlsem:::loglikelihood_lms(parameters, model, dat, P, m = 16)
    # function (parameters, model, dat, P, m = 16, ...) 
    # {
    mod.filled <- nlsem::fill_model(model = model, parameters = parameters)
    mod.filled$matrices
    k <- nlsem:::get_k(mod.filled$matrices$class1$Omega)
    quad <- nlsem:::quadrature(m, k)
    V <- quad$n
    # if (k == 0) {
    #   V <- as.data.frame(0)
    # }
    res0 <- sapply(seq_len(nrow(V)), function(i) { # i = 1
      lls <- sum(
        mvtnorm::dmvnorm(dat, 
                # mean = nlsem:::mu_lms(model = mod.filled, z = V[i, ]), 
                # sigma = nlsem:::sigma_lms(model = mod.filled, z = V[i, ]), 
                         mean = MU_lms(model = mod.filled, z = V[i, ]),
                         sigma = SIGMA_lms(model = mod.filled, z = V[i, ]),
                         log = TRUE) * P[, i])
      lls
    })
    res <- sum(res0)
    -res
    # }
  }
  ## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  # Compute parameters acheiving Maximum Likelihood ------------------------
  if(TRUE){
    cat("Doing maximization-step \n")
  }
  
  # est <- optim(
  #   par = parameters,
  #   fn = nlsem:::loglikelihood_lms,
  #   model = model,
  #   dat = dat,
  #   P = P,
  #   # upper = model$info$bounds$upper,
  #   # lower = model$info$bounds$lower,
  #   # method = "L-BFGS-B",
  #   control = control)
  
  
  est <- nlminb(
    start=parameters, 
    # objective=nlsem:::loglikelihood_lms, 
    objective=LL_lms, 
    dat=dat,
    model=model, 
    P=P, 
    upper=model$info$bounds$upper,
    lower=model$info$bounds$lower, 
    control=control)
  
  names(est) <- gsub("value", "objective", names(est))
  # }
  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  # Calculate Hessian matrix for SE.-------------------------------------
  if (neg.hessian == TRUE) {
    est$hessian <- nlme::fdHess(pars = est$par, 
                                fun = LL_lms, 
                                model = model, dat = dat, P = P)$Hessian
  }
  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  m.step <- est
  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  ll.new <- m.step$objective
  ll.ret <- c(ll.ret, ll.new)
  par.new <- unlist(m.step$par)
  num.iter <- num.iter + 1
  
  if (num.iter == max.iter) {
    warning("Maximum number of iterations was reached. EM algorithm might not have converged.")
    break
  }
  
  # Check Difference in likelihoods for convergence and stop iterations
  if (abs(ll.old - ll.new) < convergence) 
    run <- FALSE
}
# <<<<<<<<<<<<<<<<<< EM algorithm ends <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Finalize M estimation ----------------------------------------------------
final <- 
  nlsem:::mstep_lms(
    model = model, 
    P = P, 
    dat = data, 
    parameters = par.new, 
    neg.hessian = neg.hessian, 
    m = m, 
    optimizer = optimizer, 
    max.mstep = max.mstep)

coefficients <- final$par
names(coefficients) <- model$info$par.names
# Transform parameters back to Phi
A <- matrix(0, nrow = model$info$num.xi, ncol = model$info$num.xi)
A[lower.tri(A, diag = TRUE)] <-
  c(coefficients[grep("Phi", names(coefficients))],
    model$matrices$class1$A[lower.tri(model$matrices$class1$A, diag = TRUE)][!is.na(model$matrices$class1$A[lower.tri(model$matrices$class1$A, diag = TRUE)])])

A[upper.tri(A)] <- t(A)[upper.tri(t(A))]
Phi <- A %*% t(A)
coefficients[grep("Phi", names(coefficients))] <- 
  Phi[lower.tri(Phi, diag = F)]

# Clean result output
if (num.iter == max.iter) {
  em_convergence <- "no"
} else {
  em_convergence <- "yes"
}

info <- model$info[c("num.xi", "num.eta", "num.x", "num.y", 
                     "constraints", "num.classes")]
info$n <- nrow(data)
out <- list(model.class = class(model), coefficients = coefficients, 
            objective = -final$objective, em.convergence = em_convergence, 
            neg.hessian = final$hessian, loglikelihoods = -ll.ret, 
            info = info)


out
