# rm(list = ls())
# install.packages("nlsem")
# install.packages(c("gaussquad","mvtnorm"))

# for(i in fs::dir_ls("../R")) { source(i)}

# for(i in c("my_likelihood.R","my_samplestat_lms.R")) { source(i)}

library(data.table)
library(tidyverse)

# library(nlsem)
# library(mvtnorm)
# library(gaussquad)
# fdHess <- nlme::fdHess

for(i in c("utils.R","customized_sampleStat.R",
           "customized_Likelihood.R")) { source(i)}


data <- data.table::fread(fs::dir_ls("../data")[1])
data <- dplyr::select(data, x1:x3, x4, x5, y1) #%>%
  # mutate(
  #   x4 = if_else(x4 > mean(x4), 1, 0)
  # )

fwrite(data, "test_dt.csv", col.names = F)
writeLines("
DATA: file is test_dt.csv;
VARIABLE: names are v1-v3 z x y;

ANALYSIS:
 type = random;
  ALGORITHM=INTEGRATION;
model:

F1 by v1-v3;
[v1@0];
[F1];
F1*;

!zz by z@1;
![z@0];
! z@0;

F1z | F1 XWITH z;

y on z x F1 F1z;


F1 with x;
F1 with z;

", "test_inp.inp")


MplusAutomation::runModels("test_inp.inp")
file.show("test_inp.out")


# data <- data %>% select(x1:x3, x4, x5, y1) %>%
#   set_names(c("v1","v2","v3","z","x","y")) %>%
#   mutate(
#     z = if_else(z > mean(z), 1, 0)
#   ) %>%
#   select(y, z, x, v1:v3)

# 
model <- nlsem::specify_sem(
  num.x = 5, 
  num.y = 1,
  num.xi = 3,
  num.eta = 1,
  xi = "x1-x3,x4,x5",
  eta = "y1", 
  num.classes = 1,
  interaction = "xi1:xi2")

class(model)
model$matrices# 
model$info


specs <- nlsem:::as.data.frame.singleClass(model)
head(specs)

# User constraints
specs[specs$label %in% paste0("Lambda.x4", 1), "class1"] <- 1
specs[specs$label %in% paste0("Lambda.x5", 1), "class1"] <- 1
specs[specs$label %in% paste0("Theta.d", c(19,25)), "class1"] <- c(0,0)
specs[specs$label %in% paste0("Theta.e", c("")), "class1"] <- c(0)

specs[specs$label %in% paste0("tau", c(2,3)), "class1"] <- 
  c(mean(data$x4),mean(data$x5))

# Phi is A : cholesky decomposition
# A.chol <- t(chol(A))
specs[specs$label %in% paste0("Phi", c(5,6,9)), "class1"] <- 
  c(var(data$x4),cov(data$x5, data$x4), var(data$x5))
specs[specs["label"]%in%paste0("Phi",4:9)] <- NULL

my_model <- nlsem::create_sem(specs)

my_model$matrices
my_model$info$par.names[my_model$info$par.names == "tau1"] <- "tau"
my_model$info$par.names[my_model$info$par.names == "Omega4"] <- "Omega"

pars.start <- runif(nlsem:::count_free_parameters(my_model))
names(pars.start) <- my_model$info$par.names
pars.start["Phi1"] <- 0.5
pars.start[c("Phi2","Phi3")] <- c(0.08,0.1)

pars.start[c("Lambda.x2","Lambda.x3")] <- c(1,1)

pars.start[c("Gamma1","Gamma2","Gamma3")] <- c(0.2,0.1,0.1)

pars.start[c("Theta.d1","Theta.d7","Theta.d13")] <- c(0.1)
pars.start[c("Psi")] <- c(0.1)

pars.start[c("tau")] <- c(2.5)
pars.start[c("alpha")] <- c(0.2)
pars.start[c("nu.x2","nu.x3")] <- c(-0.3, -0.1)

names(pars.start)[names(pars.start)=="Omega4"] <- "Omega"

pars.start[c("Omega")] <- c(-0.02)

my_model$matrices$class1$Phi[1,2] <- NA
my_model$matrices$class1$Phi[2,3] <- my_model$matrices$class1$Phi[3,2]
my_model$matrices$class1$Phi[1,3] <- NA
my_model$matrices
my_model$info

saveRDS(list(data = data, model = my_model, pars.start = pars.start), "data_and_modelSpec.RDS")

# model.filled <- nlsem::fill_model(model = model, parameters = pars.start)
# model.filled$matrices

args(em)
res <- em(model = my_model, 
          data = data,
          start =  pars.start, 
          qml = F,
          verbose = T,
          convergence = 0.1, 
          max.iter = 200,
          max.singleClass = 1, 
          neg.hessian = TRUE,
          m = 16, 
          optimizer = c("nlminb")
          )

summary(res)



####
data <- fread(fs::dir_ls("../data")[1])

model <- specify_sem(
  num.x = 6, 
  num.y = 3,
  num.xi = 2,
  num.eta = 1,
  xi = "x1-x3,x4-x6",
  eta = "y1-y3", 
  num.classes = 1,
  interaction = "eta1~xi1:xi2")

model$matrices# 

specs <- as.data.frame(model)

my_model <- create_sem(specs)
my_model$matrices

pars.start <- runif(count_free_parameters(my_model))

res <- em(model = my_model, 
          data = data,
          start =  pars.start, 
          qml = F,
          verbose = T,
          convergence = 0.1, 
          max.iter = 200,
          max.singleClass = 1, 
          neg.hessian = TRUE,
          m = 16, 
          optimizer = c("nlminb")
)

summary(res)

fwrite(data, "test_dt.csv", col.names = F)

writeLines("
DATA: file is test_dt.csv;
VARIABLE: names are x1-x6 y1-y3;

ANALYSIS:
 type = random;
 ALGORITHM=INTEGRATION;

model:

F1 by x1-x3;
F2 by x4-x6;

F3 by y1-y3;
[y1@0];
[x1@0];
[x4@0];

F12 | F1 XWITH F2;

F3 on F1 F2 F12;

", "test_inp.inp")

MplusAutomation::runModels("test_inp.inp")
file.show("test_inp.out")

# -------------------------------------------------------------------------


data("PoliticalDemocracy", package = "lavaan")
dat <- as.matrix(PoliticalDemocracy[ ,c(9:11,1:8)])

model <- specify_sem(
  num.x = 8, num.y = 3, 
  num.xi = 2, num.eta = 1,
  xi = "x1-x3,y1-y4", 
  eta = "y5-y8", 
  num.classes = 1,
  interaction = "eta1~xi1:xi2")

em(model, data, start, qml = FALSE, verbose = FALSE, 
   convergence = 1e-02,
   max.iter = 100, m = 16, optimizer = c("nlminb", "optim"),
   max.mstep = 1, max.singleClass = 1, neg.hessian = TRUE)

