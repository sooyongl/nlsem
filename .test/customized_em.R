# rm(list = ls())
for(i in c("utils.R","customized_sampleStat.R",
           "customized_Likelihood.R")) { source(i)}

dmvnorm <- mvtnorm::dmvnorm
fdHess <- nlme::fdHess

inp <- readRDS("data_and_modelSpec.RDS")

data <- inp$data
inp_model <- list(
  matrices = inp$model$matrices,
  info = inp$model$info)
pars.start <- inp$pars.start


# ------------------------------------------------------
# Compute EM algorithm for LMS
# -------------------------------------------------------

# Argument inputs --------------------------------------
# These arguments come from `em` function in `nlsem` package

model = inp_model     # Model matrix
data = data          # input data
start =  pars.start  # Starting values

max.iter = 1000        # Maximum iteration for EM; default is 200
max.mstep = 1        # Iteration for mstep; default is 1
convergence = 0.01   # 
m = 16               # number of quadrature points

neg.hessian = FALSE  # for getting SE (For Hessian matrix)
optimizer = "nlminb" # Obtain MLE

max.singleClass = 1  
qml = F              # Alternative to LMS

# EM Estimation begins ----------------------------------------------------
## Name starting values (e.g., negative value for variances)  
par.new <- start

ll.ret <- NULL # Save likelihoods across iterations (from M step)
ll.new <- 0    # Likelihood every iteration (from M step)

num.iter <- 0
run <- TRUE

# Computing E and M steps until run == FALSE ------------------------------
while (run) {
  
  if (TRUE) {
    cat(paste("Iteration", num.iter + 1, "\n"))
    cat("Doing expectation-step \n")
  }
  
  ll.old <- ll.new
  par.old <- par.new # par.new will keep updated
  
  #------------------------------------------------------------------------
  # Compute E-step <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #------------------------------------------------------------------------
  
  model = model
  parameters = par.old
  dat = data
  m = m
  
  ##### Compute E-step         <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  mod.filled <- fill_model(model = model, parameters = parameters)
  
  # Calculate Hermite Quadrature <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  k <- 1 # Number of interaction effects
  quad <- quadrature(m, k)
  
  V <- quad$n
  w <- quad$w
  # Hermite Quadrature Ends <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  P <- NULL
  for (i in seq_along(w)) {
    # i = 1
    p.ij <- w[i] * mvtnorm::dmvnorm(dat, 
                                    mean = MU_lms(mod.filled, V[i, ]),
                                    sigma = SIGMA_lms(mod.filled, V[i, ]))
    P <- cbind(P, p.ij, deparse.level = 0)
  }
  P <- P/rowSums(P)
  P
  # E-step Ends         <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  #------------------------------------------------------------------------
  # Compute M-setp <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #------------------------------------------------------------------------
  
  model = model
  P = P
  dat = data
  parameters = par.old
  m = m
  optimizer = optimizer
  max.mstep = max.mstep
  control <- list(maxit = max.mstep)
  
  ## Compute Likelihood for LMS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ## This is for testing
  if(F) {

    mod.filled <- fill_model(model = model, parameters = parameters)
    
    mod.filled$matrices
    
    k <- 1
    quad <- quadrature(m, k)
    V <- quad$n # quadrature points

    res0 <- sapply(seq_len(nrow(V)), function(i) { # i = 1
      lls <- sum(
        mvtnorm::dmvnorm(dat,
                         mean = MU_lms(model = mod.filled, z = V[i, ]),
                         sigma = SIGMA_lms(model = mod.filled, z = V[i, ]),
                         log = TRUE) * P[, i])
      lls
    })
    res <- sum(res0)
    -res
    # }
  }
  ## Compute Likelihood Ends <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  # Compute parameters achieving Maximum Likelihood <<<<<<<<<<<<<<<<<<<<<<<
  if(TRUE){
    cat("Doing maximization-step \n")
  }
  
  est <- nlminb(
    start=parameters, 
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
  # Calculate Hessian matrix Ends <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  m.step <- est
  # M-setp Ends <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
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
# In my case, stopped at 445th iteration

# Finalize M estimation ----------------------------------------------------
final <- nlminb(
  start=par.new, 
  objective=LL_lms, 
  dat=data,
  model=model, 
  P=P, 
  upper=model$info$bounds$upper,
  lower=model$info$bounds$lower, 
  control=control)
names(est) <- gsub("LL", "objective", names(est))

# Final parameter estimates ------------------------------------------------
final$par

saveRDS(final, "LMS_estimates.RDS")

# Benchmark ---------------------------------------------------------------
file.show("test_inp.out")





