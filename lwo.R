
#' @title longitudinal win odds
#' @description
#' function for fitting a longitudinal probabilistic index model with logit link
#' to get longitudinal win odds
#' @param formula formula 
#' @param data data frame containing relevant variables
#' @param family family of the link function, only logit link gives win odds
#' @param id name of the column that contains cluster id. For example, patient id for repeated measurements on the same patient
#' Data are assumed to be sorted so that observations on each cluster appear as contiguous rows in data. If data is not sorted this way,
#'  the function will not identify the clusters correctly. 
#' @param visit name of the column that contains visit within each cluster. We need this because subjects forming a pair must come from the same visit
#' @param time.varname (scalar or vector of) variable names that represent time in the formula. This is IMPORTANT because
#' time variable will not be converted into pair level. 
#' @param corstr working correlation structure, one of "independence"', '"exchangeable"', '"ar1"', '"unstructured"' and '"userdefined"'
#' @param zcor Used for entering a user defined working correlation structure.
#' @param weights weight for each observation
#' @import geepack
#' @export lwo
lwo <- function (formula, 
                  data = parent.frame(), 
                  family = binomial(),
                  id,
                  visit,
                  time.varname = NULL,
                  corstr = "independence", 
                  weights,subset,na.action, start = NULL, etastart, mustart, 
                  offset, control = geese.control(...), method = "glm.fit", 
                  contrasts = NULL, waves = NULL, 
                  zcor = NULL, scale.fix = TRUE, 
                  scale.value = 1,std.err="san.se", ...) 
{
  if (missing(data)) 
    data <- environment(formula) #If data is not explicitly provided, R uses the environment in which the formula was created to find the variables (e.g., the global environment or a specific data frame).
  ##---- Prepare model matrix and response vector
  # returns a full call to the function. e.g., glm(formula = y~x data=df, family = binomial())
  # expand.dot: expands any ... arguments to include all their contents explicitly. This is useful for saving the full call for reproducing results or diagnostics.
  call <- match.call(expand.dots = TRUE)
  mf.call <- match.call(expand.dots = FALSE)
  # find position of the following terms in mf
  m <- match(c("formula", 
               "data",
               #   "id",
               #   "visit",
               #   "time",
               "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf.call), 0L) # returns 0 if not the term is not inputted
  # reduce mf to only contain the above terms
  mf.call <- mf.call[c(1L, m)] #c(1L, m) ensures the function name (first element) and the matched arguments are retained.
  # Other arguments, such as family or additional parameters, are excluded temporarily.
  mf.call$drop.unused.levels <- TRUE # If factors in the data frame have levels not present in the subset used for modeling, this avoids errors or inconsistencies during the fitting process.
  # The first element of mf (which was previously the function name, e.g., lwo) is replaced with stats::model.frame.
  # stats::model.frame(formula = y ~ x, data = df, ...)
  mf.call[[1L]] <- quote(stats::model.frame)
  
  # mf is the model frame
  # eval() executes the stats::model.frame() call in the appropriate context (parent.frame()), producing the model frame (mf).
  # The model frame includes the response variable, predictors, and any weights, offsets, or subsets specified in the arguments.
  mf <- eval(mf.call, parent.frame())
  # stats::model.frame() attaches a "terms" attribute to the resulting model frame. This "terms" object contains detailed information about the model formula, such as:
  # The response and predictor variables.
  # Their roles (response vs. predictors).
  # Any interactions or transformations specified in the formula.
  # mt is used later in the fitting process to handle predictors correctly.
  mt <- attr(mf, "terms")
  # return(list(mf,mt))
  
  ################################################################################
  # Create new data frame from input dataframe: convert outcome and covariate to pair level        
  ################################################################################
  # we need input id name and time name, because time variable will not be converted to pair level
  relevant_vars <- extract_variable_names(formula)
  outcome.label <- relevant_vars[1]
  covariate.labels <-  relevant_vars[2:length(relevant_vars)]
  
  # only form pairs at the same visit, not from different visits
  # list of all possible pairs
  tab <- table(data[,id],data[,visit])
  #table(mf[,"(id)"],mf[,"(visit)"]) # this records whether jth visit of ith subject is missing
  all.possible.pairs <- data.frame(t(combn(as.numeric(rownames(tab)),2)))
  colnames(all.possible.pairs) <- c("L","R")
  meet <- apply(all.possible.pairs,1,function(i)
  {
    # add the appearance indicator of two individuals for all visits
    # if no element is 2, then they never meet
    meet <- any(tab[i["L"],] + tab[i["R"],]==2)
    return(meet)
  })
  all.pairs <- all.possible.pairs[meet,]
  all.pairs$pair.id <- 1:nrow(all.pairs)
  
  # now for each visit, make pseudo variable values on the pair level
  # see how many visits are there
  visits <- unique(data[,visit])[order(unique(data[,visit]))]
  temp <- sapply(visits,function(i)
  {
    data.this.visit <- data[data[,visit]==i,]
    # see which pairs can be formed at this visit
    pairs.this.visit <- intersect(which(all.pairs$L %in% data.this.visit[,id]),
                                  which(all.pairs$R %in% data.this.visit[,id]))
    dat.pairs.this.visit <- all.pairs[pairs.this.visit,]
    # get information from subject 1 
    dat.L <- dplyr::left_join(dat.pairs.this.visit,
                              dplyr::select(data.this.visit,all_of(c(id,relevant_vars))),
                              by=c("L"=id))
    # get information from subject 2
    dat.R <- dplyr::left_join(dat.pairs.this.visit,
                              dplyr::select(data.this.visit,all_of(c(id,relevant_vars))),
                              by=c("R"=id))
    # calculate pseudo variable values as the difference between subject 1 (left) and subject 2 (right): right - left always
    dat.RminusL <- dat.R - dat.L 
    dat.RminusL <- dplyr::select(dat.RminusL,-c("L","R","pair.id"))
    
    # convert pair level outcome to either a win (1), a tie (0.5) or a loss (0)
    dat.RminusL[,outcome.label] <- 0.5*(dat.RminusL[,outcome.label]==0) + 1*(dat.RminusL[,outcome.label] > 0)
    
    
    # We should not convert time variable to pair difference, because we only form pairs with subjects
    # at the same time/visit, so pair level time variable is always zero
    # use the mean value of time between left and right subject
    dat.RminusL[,time.varname] <- (dat.L[,time.varname]+dat.R[,time.varname])/2
    dat.comb <- cbind(dat.pairs.this.visit,dat.RminusL)
    return(dat.comb)
  },simplify=FALSE)
  data.pair <- dplyr::bind_rows(temp)
  data.pair <- data.pair[order(data.pair$pair.id),]
  ##---- Modify model framework mf
  # since we have converted data frame from individual level to pair level
  # we have to modify the model framework accordingly
  mf.pair <- model.frame(formula, data.pair, drop.unused.levels = TRUE)
  mt.pair <- terms(mf.pair)
  
  ##---- Construct X and Y 
  Y  <- model.response(mf.pair, "numeric")
  N <- NROW(Y) # number of pairs
  X  <- if (!is.empty.model(mt.pair)) 
  {
    model.matrix(mt.pair, mf.pair, contrasts)
  }else 
  {matrix(, N, 0)}
  # remove intercept
  X <- X[,colnames(X)!="(Intercept)"]
  
  ##---- Clustering variable
  id.pair <- data.pair[,"pair.id"]
  if (is.null(id.pair)) stop("pair id variable not found.")
  
  ##---- Get family, for longitudinal win odds, only logit link is valid
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family not recognized")
  }
  if (family$family!="binomial")
  {
    warning("Longitudinal win odds is only valid with a logit link!")
  }
  
  ##---- Match method, can only return model data frame if user inputs "model.frame"
  if (identical(method, "model.frame")) 
    return(mf.pair)
  if (!is.character(method) && !is.function(method)) 
    stop("invalid 'method' argument")
  if (identical(method, "glm.fit")) 
    control <- do.call("geese.control", control)
  
  ##---- Working correlation structures
  if (corstr=="fixed" && is.null(zcor)){
    stop("When corstr is 'fixed' then 'zcor' must be given\n")
  }
  CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
               "userdefined", "fixed")
  # do partial match if users put ind instead of independence
  corstrv <- pmatch(corstr, CORSTRS, -1)
  corstr  <- CORSTRS[corstrv]
  
  ##---- Waves (not used)
  waves <- model.extract(mf, waves)
  if (!is.null(waves)) waves <- as.factor(waves)
  
  ##---- Weights for the observations (not used)
  weights <- model.weights(mf)
  if (is.null(weights)) weights <- rep(1, N)
  
  ##---- Offset for the mean and the scale (not used)
  offset <- model.offset(mf)
  if (is.null(offset)) offset <- rep(0, N)
  soffset <- rep(0, N)
  
  ## Check that factors in model do not have unused levels in data 
  ## (otherwise R crashes).
  vars <- all.vars(formula)
  stopIt <- FALSE
  for(ii in seq_along(vars)){
    vv <- vars[ii]
    if(!is.na(match(vv, names(mf))) && is.factor(mf[,vv])){
      if (length(unique(mf[,vv])) != length(levels(mf[,vv]))){
        cat("Factors not allowed to have unused levels...\n")
        cat(" Levels of factor", vv,":", paste(levels(mf[, vv]), sep=' '),"\n")  
        cat(" Used levels of factor", vv,":", paste(unique(mf[, vv]), sep=' '),"\n")
        stopIt <- TRUE
      }
    }
  }
  if (stopIt)
    stop("Can not continue...\n")
  
  
  ##---- Create dummy glm object in order for the model to fit in generic functions that have glm method
  glmFit <- glm.fit(x = X, y = Y, weights = rep(1,N), start = NULL, etastart = NULL, 
                    mustart = NULL, offset = rep(0,N), family = binomial(), 
                    control = list(), intercept = TRUE, singular.ok = TRUE)
  class(glmFit) <- "glm"
  glmFit$terms <- mt.pair
  glmFit$model <- mf.pair
  # modelmat <- model.matrix(glmFit)[-1] # no intercept
  # qqrr     <- qr(modelmat)
  # if (qqrr$rank < ncol(modelmat)){
  #   print(head(modelmat))
  #   stop("Model matrix is rank deficient; geeglm can not proceed\n")
  # }
  
  ##---- Use geese.fit to estimate parameters, the variances are wrong which will be 
  # fixed later
  ans <- geese.fit(x=X, y=Y, id=id.pair, 
                   family=family,
                   corstr=corstr,zcor = zcor,
                   offset=offset, soffset=soffset,
                   weights=weights, waves = waves, control = control, 
                   corp = NULL, b = NULL, 
                   alpha = NULL, gm = NULL,  mean.link = NULL, variance = NULL, 
                   cor.link = "identity", sca.link = "identity", link.same = TRUE, 
                   scale.fix = scale.fix, scale.value = 1)
  ans <- c(ans, list(call = call, formula = formula))
  class(ans)  <- "geese"
  ans$model.matrix <- X
  ans$id.pair <- id.pair
  # return(ans)
  ################################################################################
  # update to geese, modify the standard error vbeta
  ################################################################################
  # extract U evaluated at estimated regression coefs
  N_pairs <- length(unique(id.pair))
  betas_hat <- c(ans$beta)
  rho_hat <- ans$alpha
  x <- X
  y <- Y
  # fitted values
  mu <- plogis(c(x%*%betas_hat))
  num_visits_per_pair <- c(table(id.pair))
  # list of inverse of the correlation matrix, for every subject, number of visits may vary under missing data
  R_inv_list <- sapply(1:max(num_visits_per_pair), function(n)
  {
    mat <- generate_corMatrix(n_visits = n,rho=rho_hat,corstr = corstr)
    mat_inv <- solve(mat)
  })
  U_l <- x*c(mu*(1-mu))
  # U_l <- diag(mu*(1-mu))%*%x # this is also valid but super slow
  cum_num_visits_per_pair <- c(0,cumsum(num_visits_per_pair))
  U <- t(sapply(1:N_pairs, function(i) 
  {
    row_indexes <- (cum_num_visits_per_pair[i]+1):cum_num_visits_per_pair[i+1]
    R_inv <- R_inv_list[[length(row_indexes)]]
    V_sqrt_inv <- solve(diag(x=sqrt(mu[row_indexes]*(1-mu[row_indexes])),nrow =length(row_indexes)))
    sigma_inv <- V_sqrt_inv%*%R_inv%*%V_sqrt_inv
    res <- t(matrix(U_l[row_indexes,],ncol = ncol(x)))%*%sigma_inv%*%c(y[row_indexes]-mu[row_indexes])
    return(res)
  }
  ))
  # make U the right dimension (N_pairs * length(betas))
  U <- matrix(U,nrow = N_pairs)
  
  ########## sandwich estimator
  # the bread part which takes care of variance
  # analytically, faster than numerical 
  U.diff <- sapply(1:N_pairs, function(i)
  {
    row_indexes <- (cum_num_visits_per_pair[i]+1):cum_num_visits_per_pair[i+1]
    Z_ij <- matrix(x[row_indexes,],ncol = ncol(x))
    mu_ij <- matrix(mu[row_indexes],ncol=1)
    mu_ij_diff <- Z_ij*c(mu_ij*(1-mu_ij))
    Y_ij <- matrix(y[row_indexes],ncol=1)
    t <- nrow(Z_ij)
    R_inv <- R_inv_list[[length(row_indexes)]]
    V_sqrt <- diag(x=sqrt(mu[row_indexes]*(1-mu[row_indexes])),nrow = t)
    t_p <- expand.grid(1:t,1:length(betas_hat))
    V_sqrt_diff <- matrix(apply(t_p,1,function(i)
    {
      t_i <- i[1]; p_i <- i[2]
      mu_ti <- mu_ij[t_i]; Z_ti_pi <- Z_ij[t_i,p_i]
      return(1/2*sqrt(mu_ti*(1-mu_ti))*(1-2*mu_ti)*Z_ti_pi)
    }),byrow = FALSE,nrow = t,ncol = length(betas_hat))
    V_sqrt_inv <- solve(V_sqrt)
    V_sqrt_inv_diff <- matrix(apply(t_p,1,function(i)
    {
      t_i <- i[1]; p_i <- i[2]
      mu_ti <- mu_ij[t_i]; Z_ti_pi <- Z_ij[t_i,p_i]
      return(-1/2*(1/sqrt(mu_ti*(1-mu_ti)))*(1-2*mu_ti)*Z_ti_pi)
    }),byrow = FALSE,nrow = t,ncol = length(betas_hat))
    temp1 <- V_sqrt%*%R_inv%*%V_sqrt_inv_diff
    temp2 <- t(t(V_sqrt_diff)%*%R_inv%*%V_sqrt_inv)
    temp3 <- (temp1+temp2)*c(Y_ij-mu_ij)
    temp4 <- V_sqrt%*%R_inv%*%V_sqrt_inv%*%(-mu_ij_diff)
    temp5 <- t(Z_ij)%*%(temp3+temp4)
    return(temp5)
  },simplify = FALSE)
  U.diff.reduced <- Reduce("+",U.diff)
  U.diff.inv <- solve(U.diff.reduced)
  
  # the meat, use residuals to estimate correlation
  # use matrix multiplication, much faster
  # IDs for subject 1 and subject 2 in each pair
  # make pair identifier data frame
  pair_df <- dplyr::select(data.pair,c("L","R","pair.id"))
  IDs_subjectL <- unique(pair_df[,c("L","pair.id")])[,"L"]
  IDs_subjectR <- unique(pair_df[,c("R","pair.id")])[,"R"]
  # correlated pairs like (1,2) and (1,3), shared subject is on the same side
  shared.factor=1
  # correlated pairs like (1,2) and (2,4), shared subject is on either side
  switched.factor=1
  # self correlation, (1,2) and (1,2)
  self.factor=1
  # Compute column sums across rows of score matrix-like object for pairs with the same subject1/2
  Usum1.tmp <- rowsum(U,IDs_subjectL,reorder=FALSE)
  Usum2.tmp <- rowsum(U,IDs_subjectR,reorder=FALSE)
  Usum1  <- matrix(nrow = N, ncol = ncol(U),0)
  Usum2  <- matrix(nrow = N, ncol = ncol(U),0)
  Usum1[unique(IDs_subjectL),] <- Usum1.tmp
  Usum2[unique(IDs_subjectR),] <- Usum2.tmp
  UtUshared <- crossprod(Usum1) + crossprod(Usum2) 
  UtUswitched <- crossprod(Usum1,Usum2) + crossprod(Usum2,Usum1)
  UDiag<-crossprod(U) #Is counted twice as shared.factor, but needs to be counted as self.factor
  UtU<-shared.factor*UtUshared  + 
    (switched.factor)*UtUswitched +
    (self.factor-2*shared.factor)*(UDiag)
  
  # variance-covariance matrix of coefficients
  varcov <- U.diff.inv%*%UtU%*%U.diff.inv
  colnames(varcov) <- rownames(varcov) <- names(betas_hat)
  # replace the wrong var-cov matrix in ans
  ans$vbeta <- varcov
  
  ##---- Prepare output
  out <- glmFit # build on the skeleton of glm
  toDelete    <- c("R", "deviance", "aic", "null.deviance", "iter", 
                   "df.null", "converged", "boundary")
  out[match(toDelete,names(out))] <- NULL
  out$method <- "geese.fit"
  out$geese <- ans
  out$weights <- ans$weights
  out$coefficients <- ans$beta
  out$offset <- offset
  
  if (is.null(out$offset)){
    out$linear.predictors <- ans$model.matrix %*% ans$beta
  } else {
    out$linear.predictors <- out$offset + ans$model.matrix %*% ans$beta
  }
  
  out$fitted.values <- family(out)$linkinv(out$linear.predictors)
  out$modelInfo <- ans$model
  out$id.pair <- ans$id.pair
  out$call <- ans$call
  out$corstr    <- ans$model$corstr
  out$cor.link  <- ans$model$cor.link
  out$control   <- ans$control
  out$std.err   <- "san.se.modified"
  out$var <- varcov
  out$pair_table <- num_visits_per_pair
  out$model.matrix <- ans$model.matrix
  class(out)    <- c("lwo","geeglm", "gee", "glm", "lm")
  
  return(out)
}

#' @export
vcov.lwo <- function(object, ...){
  out <- object$geese$vbeta
  pn <- names(coef(object))
  dimnames(out) <- list(pn, pn)
  cat("\nVariance is re-estimated to account for correlation between pseudo pairs\n")
  out
}


#######################################
# summary function for object type lwo
#######################################
# object: model object of lwo class, from the lwo function
#' @export
summary.lwo <- function(object,...)
{
  z <- object
  if (is.null(z$terms)) 
    stop("invalid 'lwo' object:  no 'terms' component")
  ans <- z["call"]
  # constrcut coefficient matrix
  coefmat  <- data.frame(estimate = z$geese$beta, std.err=sqrt(diag(z$geese$vbeta)))
  coefmat$wald <- (coefmat$estimate / coefmat$std.err)^2
  coefmat$p <- pchisq(coefmat$wald, df=1,lower.tail = FALSE)
  names(coefmat) <- c("Estimate", "Std.err", "Wald", "Pr(>|W|)") ## Thanks, Achim
  ans$coeffcients <- coefmat
  ans$corstr <- z$corstr
  # correlation parameter estimate
  corrmat <- data.frame(Estimate = z$geese$alpha, Std.err=sqrt(diag(z$geese$valpha)))
  ans$corr <- corrmat
  ans$median_pair_table <- median(z$pair_table)
  class(ans) <- "summary.lwo"
  ans
}
# format printing
#' @export
print.summary.lwo <- function (x,
                               digits = max(3, getOption("digits") - 3),
                               quote = FALSE, prefix = "", ...) 
{
  if (is.null(digits)) 
    digits <- options()$digits
  else options(digits = digits)
  tmp=data.matrix(x$coeffcients)
  cat("\nCall:\n");   print(x$call)
  cat("\n Coefficients:\n");
  printCoefmat(tmp, digits = digits) 
  
  cat("\nTemporal Correlation Structure =", x$corstr, "\n")
  
  if (pmatch(x$corstr, "independence", 0) == 0) {
    cat("\nEstimated Correlation Parameters:\n")
    print(x$corr, digits = digits)
  }
  cat("Median Number of Visits Per Pseudo-pair  ", x$median_pair_table, "\n")
  invisible(x)
}

#############################################
# predict function for the lwo class object
#############################################
# object: model object of class lwo
# newdata: new data frame with all predictors in the model formula
# type: specify whether return estimated linear predictor (link) or probability (reponse)
# conf.int: whether would like associated confidence interval
# alpha: significance level for conf.int
#' @export
predict.lwo <- function (object, newdata = NULL, type = c("link", "response"), 
                         conf.int=FALSE, alpha=0.05,
                         na.action = na.pass, ...) 
{
  type <- match.arg(type)
  # form design matrix for the new data frame
  mod_mat <- model.matrix(update(formula(object), NULL ~ .), data = newdata)
  # get coefs
  coefs <- object$coefficients
  # get linear predictor
  lp <- crossprod(t(mod_mat),coefs)
  # get standard error for the linear predictor
  var <- object$geese$vbeta
  se.lp <- sqrt(diag(mod_mat%*%var%*%t(mod_mat)))
  z.alpha <- qnorm(p=1-alpha/2)
  lp.lower.CI <- lp-z.alpha*se.lp
  lp.upper.CI <- lp+z.alpha*se.lp
  if (conf.int)
  {
    if(type=="link")
    {
      ans <- cbind(lp.lower.CI,lp,lp.upper.CI)
      colnames(ans) <- c("lower.CI","estimate","upper.CI")
      return(ans)
    } else
    {
      ans <- plogis(cbind(lp.lower.CI,lp,lp.upper.CI))
      colnames(ans) <- c("lower.CI","estimate","upper.CI")
      return(ans)
    }
  } else
  {
    if(type=="link")
    {
      return(lp)
    } else 
    {
      return(plogis(lp))
    } 
  }
}

# Function to extract variable names from a formula
extract_variable_names <- function(formula) {
  # Get all terms from the formula
  terms_object <- terms(formula)
  
  # Extract the variables as they appear in the formula
  all_vars <- all.vars(formula)
  
  # Remove functions and transformations
  # Extract the variables directly involved in functions
  raw_vars <- attr(terms_object, "variables")[-1]  # Remove the response
  raw_vars <- sapply(raw_vars, function(x) {
    if (is.call(x)) all.vars(x) else as.character(x)
  })
  
  # Flatten and deduplicate the list of variables
  unique(unlist(raw_vars))
}
