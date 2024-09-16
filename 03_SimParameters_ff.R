# This file specifies the simulation setup

#rm(list=ls())

# setwd("C:/Users/Antonio/Desktop/population_adjustment_simstudy") 
effects <- "cond" #"marg"
simulation <- "ff"

#-------------------------------------------------------------------------------
#Reduced replicates to 500 (instead of 1000)
replicates <- 1000 # Monte Carlo replicates
allocation <- 1/2 # proportion of patients in active treatment (1:1 allocation ratio, vs. placebo)

no.chars <- 4 # number of baseline characteristics, these are prognostic variables
no.ems <- 2 # number of effect modifiers (subset of baseline characteristics)
pvs <- 1:4 # indices of prognostic variables
ems <- 1:2 # indices of effect modifiers 
N_AC <- c(150,300,600) # number of subjects in the AC trial
N_BC <- 600 # number of subjects in the BC trial 


#-------------------------------------------------------------------------------
#Adjusted the treatment effect below
b_trtA <- log(0.25) # baseline effect (log HR) of active treatment vs. common comparator
b_trtB <- ifelse(effects=="marg",b_trtA,log(0.4))# we chose 0.25 and 0.4 because with really common or really rare events the error marginal and conditional effects are the same.
#Q: Difference marginal - conditional are most profound when the outcome is common
#affect how big the difference between two set ups

b_X <- c(-log(0.67), -log(0.5), -log(0.33)) # effect of each prognostic variable
b_EM <- c(-log(0.67), -log(0.5), -log(0.33)) # interaction effect of each effect modifier
cens_prob <- 0.35 # censoring probability of time-to-event outcome
meanX_AC <- c(0.45,0.3,0.15) # mean of each normally-distributed covariate in AC trial
meanX_BC <- 0.6 # mean of each normally-distributed covariate in BC
sdX <- 0.2 # standard deviation of each covariate (same for AC and BC)
corX <- c(0,0.35) # covariate correlation coefficient  

# Weibull distribution parameters
weib_inv_scale <- 8.5 
weib_shape <- 1.3  

# parameter combinations for each scenario
param.combinations <- expand.grid(N_AC = N_AC, b_X = b_X, b_EM = b_EM, corX = corX, meanX_AC = meanX_AC)
      #expand.grid creates a data frame with all the value that can be formed with the combinations of all the vectors passed as argument.
pc <- param.combinations
scenarios <- nrow(param.combinations) # number of scenarios

# the rate parameter of the exponential distribution (from which censoring times are drawn from) 
# is selected to achieve specific censoring rate under active treatment at baseline. we select it 
# by simulating survival times and use optim to minimize difference between observed and target censoring rate
optim.function <- function(param, inv_scale, shape, b_trt, cens_prob, N) {
  Tlat <- -log(runif(N))/(inv_scale*exp(b_trt))^(1/shape) # latent survival time
  C <- rexp(n = N, rate = param) # censoring time
  prop_cens <- sum(Tlat > C)/N # observed censoring rate
  fit <- sum((prop_cens - cens_prob)^2) # minimize difference
  return(fit)
} 

#-------------------------------------------------------------------------------
#Created 2 cens_rates, 1 for A and 1 for B based on differing b_trts (instead of 1)
cens_rateA_func <- optim(par = 1, fn = optim.function, inv_scale = weib_inv_scale, shape = weib_shape, b_trt = b_trtA,
                   cens_prob = cens_prob, N = 1000000, method = "Brent", lower = 0, upper = 10)
cens_rateA <- cens_rateA_func$par

cens_rateB_func <- optim(par = 1, fn = optim.function, inv_scale = weib_inv_scale, shape = weib_shape, b_trt = b_trtB,
                    cens_prob = cens_prob, N = 1000000, method = "Brent", lower = 0, upper = 10)
cens_rateB <- cens_rateB_func$par
save.image(file = paste0("survival_settings_",effects,simulation,".RData")) # save simulation settings
