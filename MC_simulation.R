#setwd("C:/Users/ylong/OneDrive - LUMC/Projects/longitudinal-win-odds/sorted_codes")
rm(list = ls())
require(geepack)
require(simstudy)
require(GenOrd) # for generate ordinal samples
require(calculus) # for derivative
source("functions.R")
source("lwo.R")
load("estimands_win_odds.RData")
################ simulation set up
nsim <- 1e4
N <- 100
baseprobs <- rev(c(0.06,0.11,0.12,0.50,0.21))
visits <- c(0,4,8) # baseline is visit = 0 
# covariate effects estimated from SID trial
covs_effects_baseline <- c("trt"=0,"age"=-0.005,"pre_diarrhea"= 0.23)
corstr <- "ar1"
# working correlation 
corstr_working <- "ar1"
rho <- 0.6
corMatrix <- generate_corMatrix(n_visits = length(visits),rho=rho,corstr = corstr)
scenario <- "con"
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
set.seed(86)

for(i in 1:nsim)
{
  # seed_list[[i]] <- .Random.seed
  # for checking a specific repetition
  # .Random.seed <- seed_list[[7]]
  # generate the data (genordcat function can generate correlated ordinal outcomes)
  dat <- gen_data(N,baseprobs,covs_effects_baseline,
                  time_effect,trt_effect,time_trt_interaction,
                  visits,corMatrix)
  dat_wide <- dat$wide_format
  dat_long <- dat$long_format
  # get data on the pair level
  dat_pairs_long <- make_pseudo_pairs(data = dat_long,
                                      ind_id = "id",
                                      outcome_var = "GBS_DS",
                                      trt_var = "trt",
                                      time_var = "time",
                                      pseudo_vars = c("age","pre_diarrhea"))
  # model time categorically
  dat_pairs_long$week4_pseudo <- 1*(dat_pairs_long$time_pseudo == 4)
  dat_pairs_long$week8_pseudo <- 1*(dat_pairs_long$time_pseudo == 8) 
  mod_lwo <- lwo(GBS_DS_pseudo ~ 
                # time_pseudo+
                   week4_pseudo+week8_pseudo+
                   age_pseudo+pre_diarrhea_pseudo,
                 data=dat_pairs_long,
                 group = trt_pseudo,
                 id_pair = pair_ID,
                 id1 = ID_subject1,
                 id2 = ID_subject2,
                 corstr = corstr_working)
 
  trans_matrix <- matrix(c(1,0,0,0,0,1,1,0,0,0,1,0,1,0,0),nrow = 3,byrow = T)
  var_log_wos <- diag(trans_matrix%*%mod_lwo$var%*%t(trans_matrix))
  est_log_wos <- trans_matrix%*%c(mod_lwo$coefficients)
  # coverage on the log win odds scale
  CI_95 <- cbind(est_log_wos - qnorm(0.975)*sqrt(var_log_wos),est_log_wos + qnorm(0.975)*sqrt(var_log_wos))
  coverage_CI_95 <- log(true_wos) >= CI_95[,1] & log(true_wos) <= CI_95[,2]
  # generalized Wald test for global null hypothesis
  L <- matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0),nrow = 3,byrow = TRUE)
  W_2 <- t(L%*%mod_lwo$coefficients)%*%solve(L%*%mod_lwo$var%*%t(L))%*%(L%*%mod_lwo$coefficients)
  p_val_global_null <- pchisq(q=W_2,df=3,lower.tail = FALSE)
  results[[i]] <- list(est_log_wos,var_log_wos,CI_95,coverage_CI_95,p_val_global_null)
 if(i%%10==0)
   {print(i)
   save(results,file=paste0("results/results-comp-N",N,"-rho",rho,"-",scenario,".RData"))
 }
}

