# for(i in fs::dir_ls("../R")) { source(i)}
# install.packages("nlsem")
# install.packages(c("gaussquad","mvtnorm"))
library(nlsem)
library(data.table)

data <- fread(fs::dir_ls("../data")[1])
# 
model <- specify_sem(
  num.x = 6, num.y = 3,
  num.xi = 2, num.eta = 1,
  xi = "x1-x3,x4-x6",
  eta = "y1-y3", 
  num.classes = 1,
  interaction = "xi1:xi2")
class(model)
model$matrices# 

specs <- nlsem:::as.data.frame.singleClass(model)
head(specs)
# User constraints
# specs[specs$label %in% paste0("Lambda.x", c(2, 3, 11, 12)), "class1"] <- 1

my_model <- create_sem(specs)

pars.start <- runif(count_free_parameters(my_model))
args(em)
res <- em(model = my_model, 
          data = data,
          start =  pars.start, 
          qml = F,
          verbose = F,
          convergence = 0.1, 
          max.iter = 200,
          max.singleClass = 1, 
          neg.hessian = TRUE,
          m = 16, 
          optimizer = c("optim")
          )

summary(res)






# -------------------------------------------------------------------------


data("PoliticalDemocracy", package = "lavaan")
dat <- as.matrix(PoliticalDemocracy[ ,c(9:11,1:8)])

model <- specify_sem(
  num.x = 3, num.y = 8, num.xi = 1, num.eta = 2,
  xi = "x1-x3", eta = "y1-y4,y5-y8", 
  rel.lat = "eta1~xi1,eta2~xi1,eta2~eta1",
  # num.classes = 1, 
  interaction = "eta1~xi1:xi2")

em(model, data, start, qml = FALSE, verbose = FALSE, 
   convergence = 1e-02,
   max.iter = 100, m = 16, optimizer = c("nlminb", "optim"),
   max.mstep = 1, max.singleClass = 1, neg.hessian = TRUE)

