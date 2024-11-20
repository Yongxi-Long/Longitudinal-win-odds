##################################################
# data generating function
##################################################
# N: total sample size, 1:1 randomization
# baseprobs: marginal probabilities for ordinal categories if all prognostic effects are zero/reference group
# covs_effects_baseline: prognostic effects on the log odds ratio scale
# time_effect: effect of time per unit on the log odds ratio scale
# time_trt_interaction: time*treatment interaction effect in the log odds ratio scale
# visits: visit times
# corMatrix: correlation matrix for within-subject dependency
#' @import simstudy
gen_data <- function(N,baseprobs,covs_effects_baseline,
                     time_effect,trt_effect,time_trt_interaction,
                     visits,corMatrix)
{
  # calculate marginal probabilities for trt/control at each visit, without covariate effects
  probs_control_allvisits <- lapply(visits,function(visit)
  {
    # for N individuals at each visit
    shift_control <- time_effect*visit
    p0 <- calculate_multinomial_probabilities(shift_control,baseprobs)
    return(p0)
  })
  probs_trt_allvisits <- lapply(visits,function(visit)
  {
    shift_trt <- trt_effect*(visit > 0) + time_effect*visit + time_trt_interaction*visit
    p1 <- calculate_multinomial_probabilities(shift_trt,baseprobs)
    return(p1)
  })
  # define covariate distribution
  # 2 prognostic covariates
  # age, mean 60 sd 10, from SID-GBS trial
  def <- defData(varname = "age",dist = "normal",formula = 60,variance = 10^2)
  # 40% patients have preceding diarrhea, from SID-GBS trial
  def <- defData(def,"pre_diarrhea",dist = "binomial",formula = 0.4,variance = 1)
  # 1:1 randomization
  def <- defData(def,varname = "trt",dist = "trtAssign",formula = "1;1")
  # relationship specification, for continuous variables, mean-centered
  def <- defData(def,varname = "z",
                 formula = paste0(
                                  ifelse(covs_effects_baseline['age']>=0,"+",""),covs_effects_baseline['age'],'*age',
                                  ifelse(covs_effects_baseline['pre_diarrhea']>=0,"+",""),covs_effects_baseline['pre_diarrhea'],'*pre_diarrhea'),
                 dist="nonrandom")
  data_cov <- genData(N,def)
  data_cov$age <- round(data_cov$age,0)
  data_trt_cov <- data_cov[data_cov$trt==1,]
  data_trt <- genOrdCat(data_trt_cov,adjVar = "z",
                        baseprobs = matrix(unlist(probs_trt_allvisits),byrow = TRUE,ncol = length(baseprobs)),
                        corMatrix = corMatrix,
                        prefix = "visit",asFactor = FALSE)
  data_control_cov <- data_cov[data_cov$trt==0,]
  data_control <- genOrdCat(data_control_cov,adjVar = "z",
                            baseprobs = matrix(unlist(probs_control_allvisits),byrow = TRUE,ncol = length(baseprobs)),
                            corMatrix = corMatrix,
                            prefix = "visit",asFactor = FALSE)
  data_wide <- rbind(data_trt,data_control)
  data_wide <- data_wide[order(data_wide$id),]
  # combine and convert to long format
  data_long <- reshape2::melt(data_wide,id.vars = c("id",names(covs_effects_baseline),"z"),measure.vars = paste0("visit",1:length(visits)),
                              variable.name = "visit_cat",value.name = "GBS_DS")
  data_long <- data_long[order(data_long$id),]
  data_long$time <- rep(visits,N)
  data_sim <- list("wide_format"=data_wide,
                   "long_format"=data_long)
  return(data_sim)
}


##################################################
# calculate_multinomial_probabilities
##################################################
# function to calculate category probabilities based on provided category
# probabilities and log odds ratios for each cut-point
# shift_logit: a scalar/vector of log odds ratios for each cut-point
# p0: reference category probabilities
# common: Logical, default assumes proportionality, unless set to false and provide a vector of different log odds ratios
calculate_multinomial_probabilities <- function(shift_logit,p0,common=TRUE) {
  k <- length(p0)
  if(common==TRUE) {shift_logit <- rep(shift_logit[1],k)}
  p1 <- numeric(k)
  for (i in 1:(k-1)) {
    logit_0 <- -qlogis(cumsum(p0)[i])
    logit_1 <- logit_0 + shift_logit[i]
    p1[i] <- plogis(-logit_1) - ifelse(i==1,0,sum(p1[1:(i-1)]))
  }
  p1[k] <- 1-sum(p1,na.rm = T)
  return(p1)
}
#############################################
# function to calculate true win odds
#############################################
calculate_win_odds_one_visit <- function(p0,p1)
{
  k <- length(p0)
  # P(Y1 > Y0)
  Y1_gt_Y0 <- sum(sapply(1:(k-1), function(x) sum(p1[(x+1):k])*p0[x]),na.rm=T)
  # P(Y1 < Y0)
  Y1_lt_Y0 <- sum(sapply(2:k, function(x) sum(p1[1:(x-1)])*p0[x]),na.rm=T)
  # P(Y1 = Y0)
  Y1_eq_Y0 <- sum(sapply(1:k, function(x) p1[x]*p0[x]),na.rm=T)
  # win probability
  WP_true <- Y1_gt_Y0+0.5*Y1_eq_Y0
  # win odds
  WO_true <- (Y1_gt_Y0+0.5*Y1_eq_Y0)/(Y1_lt_Y0+0.5*Y1_eq_Y0)
  return(c("true win probability"=WP_true,"true win odds"=WO_true))
}
calculate_win_odds_multiple_visits <-  function(baseprobs,conditional_effect,
                                              time_effect,trt_effect,time_trt_interaction,visits)
{
  res <- sapply(visits,function(visit)
  {
    # for N individuals at each visit
    shift_control <- time_effect*visit + conditional_effect
    p0 <- calculate_multinomial_probabilities(shift_control,baseprobs)
    shift_trt <- trt_effect * (visit>0) + time_effect*visit + time_trt_interaction*visit +  conditional_effect
    p1 <- calculate_multinomial_probabilities(shift_trt,baseprobs)
    win_odds <- calculate_win_odds_one_visit(p0,p1)["true win odds"]
    return(win_odds)
  })
  return(res)
}
# use large simulated sample to marginal out the covariate distribution
calculate_adjusted_win_odds_multiple_visits <-  function(N_approx,baseprobs,covs_effects_baseline,
                                              time_effect,trt_effect,time_trt_interaction,visits)
{
  dat_full <- gen_data(N_approx,baseprobs,covs_effects_baseline,
                       time_effect,trt_effect,time_trt_interaction,
                       visits,corMatrix=diag(length(visits)))[["long_format"]]
  dat_full$z_rounded <- round(dat_full$z,1)
  table_z <- data.frame(table(dat_full$z_rounded))
  colnames(table_z) <- c("conditional_effect","count")
  table_z$conditional_effect <- as.numeric(levels(table_z$conditional_effect))
  table_z$freq <- table_z$count/sum(table_z$count)
  wos_conditional=sapply(table_z$conditional_effect, function(i)
    calculate_win_odds_multiple_visits(baseprobs,conditional_effect=i,
                                     time_effect,trt_effect,time_trt_interaction,visits)
  )
  # weigh each conditional win odds w.r.t the strata frequency to get marginal win odds
  wos_marginal <- apply(wos_conditional, 1, function(x) t(c(x))%*%c(table_z$freq))
  return(wos_marginal)
}

##########################################
# function to generate correlation matrix
##########################################
# n_visits: number of visits
# rho: correlation coefficient
# corstr: correlation structure, either "independent","exchangeable", or "ar1"
generate_corMatrix <- function(n_visits,rho,corstr)
{
  if(corstr %in% c("independence","ind"))
  {
    diag(n_visits)
  } else if (corstr == "ar1")
  {
    exponent <- abs(matrix(1:n_visits - 1, nrow = n_visits, ncol = n_visits, byrow = TRUE) - 
                      (1:n_visits - 1))
    rho^exponent
  } else if (corstr %in% c("exchangeable","cs","compound symmetry","exch"))
  {
    mat <- matrix(rho, nrow = n_visits, ncol = n_visits)
    diag(mat) <- 1
    mat
  } else
  {
    warning("unknown correlation structure")
    return(NULL)
  }
}

############################################
# make pseudo pairs
############################################
# data: the data (individual level to be converted to pseudo pair level), require long format
# ind_id: variable name indicating individual id in the data set (cluster id)
# outcome_var: variable name of the outcome (to be converted to pair level as win,loss,tie)
# trt_var: variable name for the treatment group indicator
# visit_var: variable name indicating the cluster pattern, usually the visit id if longitudinal data
# pseudo_vars: a vector of variable names that would be converted to pseudo pair level, usually the prognostic factors to be adjusted
# larger: logical TRUE/FALSE, indicating whether larger outcome category is better or the other way around
#' @import dplyr
make_pseudo_pairs <- function(data,ind_id,outcome_var,trt_var,time_var,pseudo_vars,
                              larger=TRUE,
                              weights=NULL)
{
  # if lower scores are better, reverse it
  if(!larger)
  {
    data[,outcome_var] <- min(data[,outcome_var])+max(data[,outcome_var])-data[,outcome_var]
  }
  # list of all possible pairs
  tab <- table(data[,ind_id],data[,time_var]) # this records whether jth visit of ith subject is missing
  all_possible_pairs <- data.frame(t(combn(as.numeric(rownames(tab)),2)))
  colnames(all_possible_pairs) <- c("ID_subject1","ID_subject2")
  # individuals who never appear at the same visit can not form a pair
  pair_actual <- apply(all_possible_pairs,1,function(i)
  {
    # add the appearance indicator of two individuals for all visits
    # if no element is 2, then they never meet
    meet <- any(tab[i["ID_subject1"],] + tab[i["ID_subject2"],]==2)
    return(meet)
  })
  all_pairs <- all_possible_pairs[pair_actual,]
  all_pairs$pair_ID <- 1:nrow(all_pairs)
  
  # now for each visit, make pseudo variable values on the pair level
  # see how many visits are there
  visits <- unique(data[,time_var])[order(unique(data[,time_var]))]
  temp <- sapply(visits,function(i)
  {
    dat_this_visit <- data[data[,time_var]==i,]
    if(!is.null(weights))
    {
      weights_this_visit <- data[data[,time_var]==i,weights]
      names(weights_this_visit) <- dat_this_visit[,ind_id]
    }
    # see which pairs can be formed at this visit
    pairs_at_this_visit <- intersect(which(all_pairs$ID_subject1 %in% dat_this_visit[,ind_id]),
                                     which(all_pairs$ID_subject2 %in% dat_this_visit[,ind_id]))
    dat_pairs_this_visit <- all_pairs[pairs_at_this_visit,]
    # get information from subject 1 
    dat1 <- dplyr::left_join(dat_pairs_this_visit,
                              dplyr::select(dat_this_visit,all_of(c(ind_id,outcome_var,trt_var,pseudo_vars))),
                              by=c("ID_subject1"=ind_id))
    # get information from subject 2
    dat2 <- dplyr::left_join(dat1,
                              dplyr::select(dat_this_visit,all_of(c(ind_id,outcome_var,trt_var,pseudo_vars))),
                              by=c("ID_subject2"=ind_id))
    # calculate pseudo variable values as the difference between subject 1 (left) and subject 2 (right)
    covs_pseudo = sapply(c(outcome_var,trt_var,pseudo_vars), function(covs_name)
    {
      covs_name_pseudo <- paste0(covs_name,"_pseudo")
      # we always compare treated subject to control subject,
      # so its X value of the treated subject - X value of the control subject
      left <- dplyr::pull(dat2,paste0(covs_name,".x"))
      right <- dplyr::pull(dat2,paste0(covs_name,".y"))
      # in case of factor dummy variables
      #left <- ifelse(is.factor(left),as.numeric(left),left)
      #right <- ifelse(is.factor(right),as.numeric(right),right)
      return(left-right)
    })
    
    colnames(covs_pseudo) <- paste0(c(outcome_var,trt_var,pseudo_vars),"_pseudo")
    dat_pairs_this_visit <- cbind(dat_pairs_this_visit,covs_pseudo)
    
    # pseudo outcome on the pair level, whether the treated subject wins over the control subject
    # tie is counted as 0.5
    dat_pairs_this_visit[,paste0(outcome_var,"_pseudo")] <- 0.5*(dat_pairs_this_visit[,paste0(outcome_var,"_pseudo")]==0)+
      1*(dat_pairs_this_visit[,paste0(outcome_var,"_pseudo")]>0)
    # for treatment effect (at all time points), comparison is only between groups
    # so visit_pseudo should multiply trt_pseudo to avoid forming within-group pairs at non-baseline visits
    # trt_pseudo = -1 means win or loss is assess by using control against treated, so y should be reversed
    # Update: it doesn't really matter as long as we keep consistent it is always the left - right
    control_to_trt_comparison <- dat_pairs_this_visit[,paste0(trt_var,"_pseudo")]== -1
    dat_pairs_this_visit[control_to_trt_comparison,paste0(outcome_var,"_pseudo")] <- 1 - dat_pairs_this_visit[control_to_trt_comparison,paste0(outcome_var,"_pseudo")]
    # corresponding pseudo covariates value should also change sign
    dat_pairs_this_visit[control_to_trt_comparison,paste0(c(trt_var,pseudo_vars),"_pseudo")] <- - dat_pairs_this_visit[control_to_trt_comparison,paste0(c(trt_var,pseudo_vars),"_pseudo")]

    dat_pairs_this_visit$visit <- which(visits==i)
    dat_pairs_this_visit$time <- i
    # time effect is converted to time*treatment interaction effect on the pair level, by multiplying the pseudo treatment indicator
    # thus only between-group pairs will have non-zero time entries
     dat_pairs_this_visit$time_pseudo <- (dat_pairs_this_visit[,paste0(trt_var,"_pseudo")]!=0)*dat_pairs_this_visit$time
    # dat_pairs_this_visit$time_pseudo <- dat_pairs_this_visit$time
    # if individual weights are supplied
    if (!is.null(weights)) 
    {
      dat_pairs_this_visit$weights <- weights_this_visit[as.character(dat_pairs_this_visit$ID_subject1)]*weights_this_visit[as.character(dat_pairs_this_visit$ID_subject2)]
    }
    return(dat_pairs_this_visit)
  },simplify=FALSE)
  dat_pairs_all_visits <- dplyr::bind_rows(temp)
  dat_pairs_all_visits <- dat_pairs_all_visits[order(dat_pairs_all_visits$pair_ID),]
  return(dat_pairs_all_visits)
}


assess_results <- function(type=NULL,N=NULL,rho=NULL,hypothesis=NULL,results=NULL)
{
  load("estimands_win_odds.RData")
  if(is.null(results))
  {
    load(paste0("C:/Users/ylong/OneDrive - LUMC/Projects/Longitudinal_win_odds/codes/results/results-",type,"-N",N,"-rho",rho,"-",hypothesis,".RData"))
  }
  if(hypothesis=="null")
  {
    true_wos <- true_wos_null
  } else if (hypothesis=="SID")
  {
    true_wos <- true_wos_SID
  } else if (hypothesis=="pos")
  {
    true_wos <- true_wos_pos
  } else
  {
    true_wos <- true_wos_con
  }
  # assess coverage of 95% CI for log win odds
  coverage <- rowMeans(sapply(results, function(i) 
  {
    i[[4]]
  }),na.rm=T)
  # power/type I error
  power_each_visit <- rowMeans(sapply(results, function(i) 
  {
    CI <- i[[3]]
    # CI does not contain 0
    CI[,1] > 0 | CI[,2] < 0
  }),na.rm=T)
  # assess bias of log win odds estimator
  est_log_wos <- sapply(results, function(i) return(i[[1]]))
  bias_log_wos <- rowMeans(est_log_wos,na.rm = TRUE) - log(true_wos)
  # convert to odds scale
  est_wos <- exp(est_log_wos)
  bias_wos <- rowMeans(est_wos,na.rm = TRUE) - true_wos
  # on the probability scale it should be unbiased
  bias_MWS <- rowMeans(est_wos/(1+est_wos),na.rm = TRUE) - true_wos/(1+true_wos)
  
  # assess estimated variance
  # "true" variance by directly calculating the variance of MC estimates
  var_log_wos_MC <- apply(est_log_wos, 1, function(x) var(x,na.rm = TRUE))
  # average estimated variance
  var_log_wos_avg <- rowMeans(sapply(results, function(i) return(i[[2]])),na.rm = TRUE)
  # global power
  power_global <- mean(sapply(results, function(i) i[[5]]<0.05))
  print(paste0("global power is: ",power_global))
  res_assessed <- data.frame(rbind(coverage,power_each_visit,true_wos,bias_log_wos,bias_wos,bias_MWS,
                                   var_log_wos_MC,var_log_wos_avg))
  colnames(res_assessed) <- paste0("visit",1:length(power_each_visit))
  return(res_assessed)
  
}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
