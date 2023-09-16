#' @describeIn 
#' Calcuate Likelihood values given parameter values for LMS
#' Customized function of `nlsem:::loglikelihood_lms`
LL_lms <- function (parameters, model, dat, P, m = 16, ...) 
{
  mod.filled <- fill_model(model = model, parameters = parameters)
  
  k <- 1
  quad <- quadrature(m, k)
  V <- quad$n
  
  if (k == 0) 
    V <- as.data.frame(0)
  
  res0 <- sapply(seq_len(nrow(V)), function(i) {
    lls <- sum(mvtnorm::dmvnorm(dat, 
                       mean = MU_lms(model = mod.filled, z = V[i, ]), 
                       sigma = SIGMA_lms(model = mod.filled,z = V[i, ]), 
                       log = TRUE) * P[, i])
    lls
  })
  res <- sum(res0)
  -res
}