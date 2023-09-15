# nlsem:::loglikelihood_lms

LL_lms <- function (parameters, model, dat, P, m = 16, ...) 
{
  mod.filled <- nlsem:::fill_model(model = model, parameters = parameters)
  k <- nlsem:::get_k(mod.filled$matrices$class1$Omega)
  quad <- nlsem:::quadrature(m, k)
  V <- quad$n
  
  if (k == 0) 
    V <- as.data.frame(0)
  
  res0 <- sapply(seq_len(nrow(V)), function(i) {
    lls <- sum(dmvnorm(dat, 
                       mean = MU_lms(model = mod.filled, z = V[i, ]), 
                       sigma = SIGMA_lms(model = mod.filled,z = V[i, ]), 
                       log = TRUE) * P[, i])
    lls
  })
  res <- sum(res0)
  -res
}