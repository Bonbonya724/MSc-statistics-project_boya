#load(file="survival_settings_condff.RData")
load(file="survival_settings_margff.RData")
library('SimDesign')
set.seed(444)  # set random seed
scenario.id <- 1
replicates <- 2000

#I have used this file for scenario 
#55, 56, 58, 59, 61, 62, 64, 65, 67, 68, 70,
#71, 73, 74, 76, 77, 79, 80, 81, 82, 85, 88, 89,
#91, 94, 95, 97, 100, 103, 104, 106, 107, 108, 109

# This file generates the simulation study data



if (!dir.exists("Data2/")) {dir.create("Data2/")}

# package with copula function to simulate the covariates
if (!require(simstudy)) {install.packages("simstudy"); library(simstudy)}
# package for data manipulation
if (!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
# Cox regression to summarize the outcomes for the ALD trial in terms of log HR and its variance
if (!require(survival)) {install.packages("survival"); library(survival)}



#Datasets below varied in:
#1.Level correlation between covariates (corX 0 or 0.35)

#-------------------------------------------------------------------------------
#Replaced cens_rate with cens_rateA and censrateB
gen.data <- function(no.chars, no.ems, N_AC, N_BC, b_trt_A, b_trt_B, 
                     b_X, b_EM_A, b_EM_B, meanX_AC, meanX_BC, sdX, weib_shape, 
                     weib_inv_scale, cens_rateA, cens_rateB, corX, allocation) {
  R <- matrix(corX, nrow = no.chars, ncol = no.chars) # set correlation matrix
  diag(R) <- rep(1, no.chars)#Diagonal is of course 1
  states.data <- .Random.seed # track random-number-generator state
  N_AC_A <- floor(N_AC*allocation) # number of patients under A in AC
  N_AC_C <- N_AC - N_AC_A # number of patients under C in AC
  N_BC_B <- floor(N_BC*allocation) # number of patients under B in BC
  N_BC_C <- N_BC - N_BC_B # number of patients under C in BC
  # simulate correlated (or uncorrelated) continuous covariates from a Gaussian copula
  # with normally-distributed marginals and the specified correlation structure  
  X_AC_A <- as.matrix(genCorGen(n = N_AC_A, nvar = no.chars, params1 = as.numeric(rep(meanX_AC, no.chars)), 
                                params2 = as.numeric(rep(sdX^2, no.chars)), dist = "normal", corstr = "cs", 
                                corMatrix = R, wide = TRUE)[,-1]) # patients under A in trial AC 
  #N is the number of patients simulated, no.chars is 4
  X_AC_C <- as.matrix(genCorGen(n = N_AC_C, nvar = no.chars, params1 = as.numeric(rep(meanX_AC, no.chars)), 
                                params2 = as.numeric(rep(sdX^2, no.chars)), dist = "normal", corstr = "cs", 
                                corMatrix = R, wide = TRUE)[,-1]) # patients under C in trial AC 
  X_BC_B <- as.matrix(genCorGen(n = N_BC_B, nvar = no.chars, params1 = as.numeric(rep(meanX_BC, no.chars)), 
                                params2 = as.numeric(rep(sdX^2, no.chars)), dist = "normal", corstr = "cs", 
                                corMatrix = R, wide = TRUE)[,-1]) # patients under B in trial BC 
  X_BC_C <- as.matrix(genCorGen(n = N_BC_C, nvar = no.chars, params1 = as.numeric(rep(meanX_BC, no.chars)), 
                                params2 = as.numeric(rep(sdX^2, no.chars)), dist = "normal", corstr = "cs", 
                                corMatrix = R, wide = TRUE)[,-1]) # patients under C in trial BC 
  # set log hazards for each patient
  betaX_AC_A <- rep(b_trt_A, N_AC_A) # intercepts and treatment effects
  betaX_AC_C <- rep(0,N_AC_C)
  betaX_BC_B <- rep(b_trt_B, N_BC_B)
  betaX_BC_C <- rep(0,N_BC_C)
  col.names <- NULL
  for (k in 1:no.chars) {
    col.names <- c(col.names, paste0('X', k))
    betaX_AC_A <- betaX_AC_A + b_X*X_AC_A[,k] # prognostic variable effects
    betaX_AC_C <- betaX_AC_C + b_X*X_AC_C[,k]
    betaX_BC_B <- betaX_BC_B + b_X*X_BC_B[,k]
    betaX_BC_C <- betaX_BC_C + b_X*X_BC_C[,k]
    if (k <= no.ems) {
      betaX_AC_A <- betaX_AC_A + b_EM_A*X_AC_A[,k] # effect modification
      betaX_BC_B <- betaX_BC_B + b_EM_B*X_BC_B[,k]
    }
  }
  X_AC = as.data.frame(rbind(X_AC_A, X_AC_C))
  colnames(X_AC) <- col.names
  betaX_AC = c(betaX_AC_A, betaX_AC_C)
  # Generate AC Weibull-Cox-distributed latent event times according to Bender et al. (2005)
  Tlat = -log(runif(n = N_AC))/(weib_inv_scale*exp(betaX_AC))^(1/weib_shape) 
  
  #-------------------------------------------------------------------------------
  #Replaced cens_rate with cens_rateA 
  C = rexp(n = N_AC, rate = cens_rateA) # AC censoring times
  # AC follow-up times and event indicators
  time = pmin(Tlat, C)
  status = as.numeric(Tlat <= C)
  trt <- c(rep("A", N_AC_A), rep("C", N_AC_C)) # treatment assignment
  IPD.AC <- as.data.frame(cbind(trt, X_AC, time, status))
  X_BC = as.data.frame(rbind(X_BC_B, X_BC_C))
  colnames(X_BC) <- col.names
  betaX_BC = c(betaX_BC_B, betaX_BC_C)
  
  # Generate BC Weibull-Cox-distributed latent event times
  Tlat = -log(runif(n = N_BC))/(weib_inv_scale*exp(betaX_BC))^(1/weib_shape) 
  #-------------------------------------------------------------------------------
  #Replaced cens_rate with censrateB
  C = rexp(n = N_BC, rate = cens_rateB) # BC censoring times
  # BC follow-up times and event indicators
  time = pmin(Tlat, C)
  status = as.numeric(Tlat <= C)
  trt <- c(rep("B", N_BC_B), rep("C", N_BC_C)) # treatment assignment
  IPD.BC <- as.data.frame(cbind(trt, X_BC, time, status))
  IPD.BC$trt <- factor(IPD.BC$trt, levels = c("C","B")) 
  # aggregate the data for the BC trial 
  ALD.BC <- as.data.frame(cbind(
    # Trial mean covariate stats and summary outcomes (log HR and variance) for BC
    summarise(IPD.BC, mean.X1 = mean(X1), mean.X2 = mean(X2), mean.X3 = mean(X3), mean.X4 = mean(X4),
              logHR_B = summary(coxph(Surv(time, status)~trt, data  =  IPD.BC))$coef[1],
              var_logHR_B = vcov(coxph(Surv(time, status)~trt, data  =  IPD.BC))[[1]],
              HR_B = summary(coxph(Surv(time, status)~trt, data  =  IPD.BC))$coef[2])))    
  list(IPD.AC, IPD.BC, ALD.BC, states.data)
}

for (i in scenario.id:scenarios) {
  #print(i) removing print as it makes code slow
  start_time <- Sys.time()  # Capture start time
  cat("Running scenario", i, "of", scenarios, "\n")
  
  IPD.AC <- IPD.BC <- ALD.BC <- vector(mode = "list", replicates)
  states.data <- vector(mode = "list", replicates + 1) # random-number-generator states
  for (j in 1:replicates) {
    #-------------------------------------------------------------------------------
    #Replaced b_trt with b_trtA and b_trtB in the function below
    #Replaced cens_rate with cens_rateA
    #Added cens_rateB
    gen.datasets <- gen.data(no.chars = no.chars, no.ems = no.ems, N_AC = pc$N_AC[i], N_BC = N_BC,
                             b_trt_A = b_trtA, b_trt_B = b_trtB, b_X = pc$b_X[i], b_EM_A = pc$b_EM[i],
                             b_EM_B = pc$b_EM[i], meanX_AC = pc$meanX_AC[i], meanX_BC = meanX_BC, sdX = sdX,
                             weib_shape = weib_shape, weib_inv_scale = weib_inv_scale, cens_rateA = cens_rateA,
                             cens_rateB  =  cens_rateB, corX = pc$corX[i], allocation = allocation)
    IPD.AC[[j]] <- gen.datasets[[1]]
    IPD.BC[[j]] <- gen.datasets[[2]]
    ALD.BC[[j]] <- gen.datasets[[3]]
    states.data[[j]] <- gen.datasets[[4]]
  }
  states.data[[replicates + 1]] <- .Random.seed # final random-number-generator state
  file.id <- paste0("N_AC", pc$N_AC[i], "b_X", round(pc$b_X[i], digits = 2), 
                    "b_EM", round(pc$b_EM[i], digits = 2), 
                    "meanX_AC", pc$meanX_AC[i], "corX", pc$corX[i]) 
  save(IPD.AC, file = paste0("Data2/",effects,simulation,"IPD_AC_", file.id, ".RData"))
  save(IPD.BC, file = paste0("Data2/",effects,simulation,"IPD_BC_", file.id, ".RData"))
  save(ALD.BC, file = paste0("Data2/",effects,simulation,"ALD_BC_", file.id, ".RData"))
  save(states.data, file = paste0("Data2/States_", file.id, ".RData"))
  
  end_time <- Sys.time()  # Capture end time
  cat("Completed scenario", i, "in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes.\n")
}  


source('02_Functions.R') # load MAIC functions
if (!dir.exists("Results/")) {dir.create("Results/")}
if (!dir.exists("Results/MAIC")) {dir.create("Results/MAIC")}

IPD.AC.all <- IPD.BC.all <- ALD.BC.all <- vector(mode = "list", scenarios)

#scenarios
for (j in scenario.id:scenarios) {
  file.id <- paste0("N_AC", pc$N_AC[j], "b_X", round(pc$b_X[j], digits = 2), 
                    "b_EM", round(pc$b_EM[j], digits = 2),"meanX_AC", pc$meanX_AC[j], "corX", pc$corX[j]) 
  load(paste0("Data2/",effects,simulation,"IPD_AC_", file.id, ".RData"))
  load(paste0("Data2/",effects,simulation,"IPD_BC_", file.id, ".RData")) 
  load(paste0("Data2/",effects,simulation,"ALD_BC_", file.id, ".RData"))  

  
  
  
  maic.aess2 <- rep(NA,replicates)
  d.AB.maic2 <- rep(NA,replicates)
  var.d.AB.maic2 <- rep(NA,replicates)
  
  for(i in 1:replicates){
    
    data.AC = IPD.AC[[i]] # at replicate level
    data.BC = ALD.BC[[i]]
    
    data.AC$trt <- factor(data.AC$trt, levels = c("C","A"))
    AC.ems <- data.AC[,1 + ems] # column 1 of IPD is treatment indicator 
    maic.weights <- maic(A.X = AC.ems, B.summary = data.BC[ems]) # maic weights through method of moments
    #maic.weights.adj <- maic.weights + 1e-2
    #maic.aess <- approx.ess(maic.weights.added)
    #maic.weights.sort <- sort(maic.reweights, TRUE)[10]
    #maic.weights.adj <- (maic.weights)*(1/maic.weight.sort) + 1e-2
    #maic.reweights <- maic.weights*(length(maic.weights)/(sum(maic.weights)))
    maic.aess2[i] <- approx.ess(maic.weights)  # approximate effective sample size
    #maic.weights.sort <- sort(maic.reweights, TRUE)[2]
    #maic.weights.adj <- maic.reweights + maic.weights.sort/(length(maic.weights)) + 1e-8
    #maic.weights.adj <- maic.reweights +  1
    #maic.weights.adj <- maic.reweights + maic.weights.sort* runif(length(maic.weights), 1e-4, 1e-3)
    
    #fit weighted Cox proportional hazards model using robust=TRUE
    #maic.weights.adj <- c(rep(1e-10,53),5.428771e-06, rep(1e-10,58), 2.495719e-09,rep(1e-10,37))
    
    #============================================================================
    #Here in the future a function could be added for bootstrapping SEs.
    data.AC.design <- svydesign(~0, data = data.AC)
    # data.AC.pos.design <- subset(data.AC.design, maic.weights>1e-6)
    # maic.pos.weights <- subset(maic.weights, maic.weights>1e-6)
    
    
    tryCatch({
      outcome.fit.maic <- svycoxph(Surv(time, status)~trt, design = data.AC.design, weights = maic.weights, data = data.AC)
      # outcome.fit.maic <- svycoxph(Surv(time, status)~trt, design = data.AC.design, weights = maic.weights.added, data = data.AC)
      # outcome.fit.maic <- svycoxph(Surv(time, status)~trt, design = data.AC.pos.design, weights = maic.pos.weights, data = data.AC)
      d.AC.maic <- quiet(summary(outcome.fit.maic)$coef[1])
      var.d.AC.maic <- vcov(outcome.fit.maic)[[1]] # sandwich-like variance estimator for A vs. C
      d.BC.maic <- with(data.BC,logHR_B)
      var.d.BC.maic <- with(data.BC, var_logHR_B)
      d.AB.maic2[i] <- d.AC.maic - d.BC.maic # ATE for A vs. B
      var.d.AB.maic2[i] <- var.d.AC.maic + var.d.BC.maic # Variance of A vs. B treatment effect
      #list(d.AB.maic, var.d.AB.maic, maic.aess)
    }, error=function(e){message(paste("Scenario", j, "replicate", i, "failed, but it is ok!"),e)}
    )
  }
  
  maic.aess <- head(maic.aess2[!is.na(d.AB.maic2)],1000)
  d.AB.maic <- head(d.AB.maic2[!is.na(d.AB.maic2)],1000)
  var.d.AB.maic <- head(var.d.AB.maic2[!is.na(d.AB.maic2)],1000)
  
  
  
  file.id <- paste0(effects,simulation, "N_AC", pc$N_AC[j], "b_X", round(pc$b_X[j], digits = 2), 
                    "b_EM", round(pc$b_EM[j], digits = 2),"meanX_AC", pc$meanX_AC[j], "corX", pc$corX[j]) 
  
  means <- d.AB.maic
  variances <- var.d.AB.maic
  approx.ess.maic <- maic.aess
  save(means, file = paste0("Results/MAIC/means_", file.id, ".RData"))
  save(variances, file = paste0("Results/MAIC/variances_", file.id, ".RData"))  
  save(approx.ess.maic, file = paste0("Results/MAIC/aess_", file.id, ".RData")) 
  
  
}







