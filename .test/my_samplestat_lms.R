# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# nlsem:::mu_lms

MU_lms <- function (model, z) {
  # z = V[8, ]
  # model = mod.filled
  
  matrices <- model$matrices$class1
  k <- nlsem:::get_k(matrices$Omega)
  
  matrices$A <- t(chol(matrices$Phi))
  n <- nrow(matrices$A)
  if (k < n & k > 0) {
    z.1 <- c(z, rep(0, n - k))
  }
  # else if (k == 0) {
  #   z.1 <- rep(0, n - k)
  # }
  # else z.1 <- z
  
  A.z <- matrices$A %*% z.1
  mu.x <- matrices$nu.x + matrices$Lambda.x %*% (matrices$tau + A.z)
  
  mu.y <- 
    matrices$nu.y + 
    matrices$Lambda.y %*% (matrices$alpha + matrices$Gamma %*% (matrices$tau + A.z) + 
    t(matrices$tau + A.z) %*% matrices$Omega %*% (matrices$tau + A.z))
  
  mu <- c(mu.x, mu.y)
  mu
}

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# nlsem:::sigma_lms
SIGMA_lms <- function (model, z) {
  # z = V[8, ]
  # model = mod.filled
  
  matrices <- model$matrices$class1
  k <- nlsem:::get_k(matrices$Omega)
  
  matrices$A <- t(chol(matrices$Phi))
  n <- nrow(matrices$A)
  if (k < n & k > 0) {
    z.1 <- c(z, rep(0, n - k))
  }
  # else if (k == 0) {
  #   z.1 <- rep(0, n - k)
  # }
  # else z.1 <- z
  
  
  
  A.z <- matrices$A %*% z.1
  d.mat <- nlsem:::get_d(n = n, k = k)
  
  Lx.A <- matrices$Lambda.x %*% matrices$A
  temp <- matrices$Gamma %*% matrices$A + t(matrices$tau + A.z) %*% matrices$Omega %*% matrices$A
  s11 <- Lx.A %*% d.mat %*% t(Lx.A) + matrices$Theta.d
  s12 <- Lx.A %*% d.mat %*% t(temp) %*% t(matrices$Lambda.y)
  s21 <- t(s12)
  s22 <- 
    matrices$Lambda.y %*% temp %*% d.mat %*% t(temp) %*% t(matrices$Lambda.y) + 
    matrices$Lambda.y %*% matrices$Psi %*% t(matrices$Lambda.y) + 
    matrices$Theta.e
  
  sigma <- rbind(cbind(s11, s12), cbind(s21, s22))
  
  if (!isSymmetric(sigma)) 
    stop("Sigma has to be symmetric.")
  
  sigma
}




