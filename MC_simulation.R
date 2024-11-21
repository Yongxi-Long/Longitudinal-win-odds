#setwd("C:/Users/ylong/OneDrive - LUMC/Projects/longitudinal-win-odds/sorted_codes")
rm(list = ls())
require(geepack)
require(simstudy)
require(GenOrd) # for generate ordinal samples
require(dplyr)
source("functions.R")
source("lwo.R")
load("estimands_win_odds.RData")
################ simulation set up
nsim <- 1e4
# for data generation
# total sample size
N <- 50
# time point for each visit
visits <- c(0,4,8) # baseline is visit = 0 
#--- Derived from SID trial
# baseline state occupancy probabilities
baseprobs <- rev(c(0.06,0.11,0.12,0.50,0.21)) # reversed so larger outcomes are better
# covariate effects 
covs_effects_baseline <- c("trt"=0,"age"=-0.005,"pre_diarrhea"= 0.23)
#--- Derived from SID trial
# correlation structure and correlation coefficient 
corstr <- "ar1"
rho <- 0.6
corMatrix <- generate_corMatrix(n_visits = length(visits),rho=rho,corstr = corstr)
# simulation scenario, one of null, SID, pos and con
scenario <- "pos"
time_effect <- switch(scenario,
                      "null" = 0.15,
                      "SID" = 0.15,
                      "pos" = 0.15,
                      "con" = 0
) # beta1, positive time trend towards higher score (self-healing)
trt_effect <- switch(scenario,
                     "null" = 0,
                     "SID" = 0,
                     "pos" = 0,
                     "con" = 1
)
time_trt_interaction <- switch(scenario,
                               "null" = 0,
                               "SID" = -0.02,# -0.02 comes from SID trial (after week 1) which is very neutral
                               "pos" = 0.1,
                               "con" = 0
)
true_wos <- switch (scenario,
                    "null" = true_wos_null,
                    "SID" = true_wos_SID,
                    "pos" = true_wos_pos,
                    "con" = true_wos_con
)
# for model fitting
set.seed(86)
# working correlation 
corstr_working <- "ar1"
# store the simulation results
results <- vector(mode = "list",length = nsim)

for(i in 1:nsim)
{

  dat <- gen_data(N,baseprobs,covs_effects_baseline,
                  time_effect,trt_effect,time_trt_interaction,
                  visits,corMatrix)
  dat_wide <- dat$wide_format
  dat_long <- dat$long_format
  
  # create indicator for week 4 and week 8
  dat_long <- dat_long |>
    mutate(w4 = 1*(time==4),w8 = 1*(time==8))
  # model time categorically
  # if scenario is constant treatment effect, then no interaction terms in the model
  if(scenario == "con")
  {
    mod_lwo <- lwo(GBS_DS ~ trt + 
                     # w4 + w8 + # cannot have main effect for time
                     # trt:w4 + trt:w8+
                     age + pre_diarrhea ,
                   data = dat_long,
                   id = "id",
                   visit = "visit_cat",
                 #  time.varname = c("w4","w8"),
                   corstr = corstr_working)
    
    trans_matrix <- matrix(c(1,0,0,
                             1,0,0,
                             1,0,0),nrow = 3,byrow = T)
    L <- matrix(c(1,0,0,
                  1,0,0,
                  1,0,0),nrow = 3,byrow = T)
    p_val_global_null <- 2*pnorm(abs(mod_lwo$coefficients["trt"]/sqrt(mod_lwo$var["trt","trt"])),
                               lower.tail = FALSE)
  } else
  {
    mod_lwo <- lwo(GBS_DS ~ trt + 
                     # w4 + w8 + # cannot have main effect for time
                     trt:w4 + trt:w8+
                     age + pre_diarrhea ,
                   data = dat_long,
                   id = "id",
                   visit = "visit_cat",
                   time.varname = c("w4","w8"),
                   corstr = corstr_working)
    
    trans_matrix <- matrix(c(1,0,0,0,0,
                             1,0,0,1,0,
                             1,0,0,0,1),nrow = 3,byrow = T)
    L <- matrix(c(1,0,0,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow = 3,byrow = T)
    # generalized Wald test for global null hypothesis
    W_2 <- t(L%*%mod_lwo$coefficients)%*%solve(L%*%mod_lwo$var%*%t(L))%*%(L%*%mod_lwo$coefficients)
    p_val_global_null <- pchisq(q=W_2,df=3,lower.tail = FALSE)
  }
  est_log_wos <- trans_matrix%*%c(mod_lwo$coefficients)
  var_log_wos <- diag(trans_matrix%*%mod_lwo$var%*%t(trans_matrix))
  # coverage on the log win odds scale
  CI_95 <- cbind(est_log_wos - qnorm(0.975)*sqrt(var_log_wos),est_log_wos + qnorm(0.975)*sqrt(var_log_wos))
  coverage_CI_95 <- log(true_wos) >= CI_95[,1] & log(true_wos) <= CI_95[,2]
  results[[i]] <- list(est_log_wos,var_log_wos,CI_95,coverage_CI_95,p_val_global_null)
  if(i%%10==0)
   {print(i)
  # save(results,file=paste0("results/results-comp-N",N,"-rho",rho,"-",scenario,".RData"))
 }
}

# # see bias
# rowMeans(sapply(results, function(i) as.vector(i[[1]])-log(true_wos)))
# # see coverage
# rowMeans(sapply(results, function(i) i[[4]]))
