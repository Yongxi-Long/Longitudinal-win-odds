# corstr a character string specifying the correlation
#'     structure. The following are permitted: '"independence"',
#'     '"exchangeable"', '"ar1"', '"unstructured"' and '"userdefined"'
# std.err para is disabled because this way we only take into account
#' temporal correlation not pairwise correlation
#' # scale.fix Omitting this argument and allowing gee to estimate the scale parameter generates a quasibinomial model,
#' #  which is an adjustment for overdispersion. Because binary data cannot be overdispersed, 
#' # a quasibinomial model is never appropriate for binary data. 
#' @import geepack
#' @export lwo
lwo <- function (formula, family = binomial(), data = parent.frame(), 
                  group, id_pair,id1,id2,corstr = "independence", 
          weights,subset,na.action, start = NULL, etastart, mustart, 
          offset, control = geese.control(...), method = "glm.fit", 
          contrasts = NULL, waves = NULL, 
          zcor = NULL, scale.fix = TRUE, 
          scale.value = 1,std.err="san.se", ...) 
{
  ##---- Prepare model matrix and response vector
  # returns a call in which all of the specified arguments are specified by their full names.
  # expand.dot: logical. Should arguments matching ... in the call be included or left as a ... argument?
  if (missing(data)) 
  data <- environment(formula)
  call <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)
  # find position of the following terms in mf
  m <- match(c("formula", "data","id_pair","id1","id2","group",
               "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0L)
  # reduce mf to only contain the above terms
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  # mf is the model data frame
  mf <- eval(mf, parent.frame())
  
  ##---- Construct X and Y
  mt <- attr(mf, "terms")
  Y  <- model.response(mf, "numeric")
  N <- NROW(Y)
  X  <- if (!is.empty.model(mt)) 
  {
    model.matrix(mt, mf, contrasts)
  }else 
  {matrix(, NROW(Y), 0)}
  # if user specify no intercept by add -1, then do nothing, otherwise need to change
  # individual level intercept to pair level intercept: 1 between group pair and 0 within group pair
  if (colnames(X)[1]=="(Intercept)")
  {
    group <- model.extract(mf,group)
    X[,1] <- group
    colnames(X)[1] <- "(Intercept.pair)"
  }

  ##---- Get family, for longitudinal win odds, only logit link is valid
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (family$family!="binomial")
  {
    warning("Longitudinal win odds is only valid with a logit link!")
  }
  
  ##---- Match method, can only return model data frame if user inputs "model.frame"
  if (identical(method, "model.frame")) 
    return(mf)
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
  
  ##---- Clustering variable
  id_pair <- model.extract(mf, id_pair)
  if (is.null(id_pair)) stop("pair id variable not found.")
  
  ##---- Waves
  waves <- model.extract(mf, waves)
  if (!is.null(waves)) waves <- as.factor(waves)
  
  ##---- Weights for the observations
  weights <- model.weights(mf)
  if (is.null(weights)) weights <- rep(1, N)
  
  ##---- Offset for the mean and the scale
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
  glmFit$terms <- mt
  glmFit$model <- mf
  modelmat <- model.matrix(glmFit)
  qqrr     <- qr(modelmat)
  if (qqrr$rank < ncol(modelmat)){
    print(head(modelmat))
    stop("Model matrix is rank deficient; geeglm can not proceed\n")
  }
  
  ans <- geese.fit(x=X, y=Y, id=id_pair, offset=offset, soffset=soffset,
                   weights=weights, waves = waves, 
                   zcor = zcor, corp = NULL, control = control, b = NULL, 
                   alpha = NULL, gm = NULL, family=family, mean.link = NULL, variance = NULL, 
                   cor.link = "identity", sca.link = "identity", link.same = TRUE, 
                   scale.fix = scale.fix, scale.value = 1, corstr=corstr)
  ans <- c(ans, list(call = call, formula = formula))
  class(ans)  <- "geese"
  ans$X       <- X
  ans$id_pair <- id_pair
  ans$weights <- weights
  
  #####################################################
  # update to geese, modify the standard error vbeta
  #####################################################
  # extract U evaluated at estimated regression coefs
  # N_pairs <- choose(N,2) # not necessarily if missed visits, some individuals never meet
  N_pairs <- length(unique(id_pair))
  betas_hat <- c(ans$beta)
  rho_hat <- ans$alpha
  x <- X
  y <- Y
  # fitted values
  mu <- plogis(c(x%*%betas_hat))
  num_visits_per_pair <- c(table(id_pair))
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
  pair_df <- data.frame(
    ID_subject1 = model.extract(mf,id1),
    ID_subject2 = model.extract(mf,id2),
    pair_ID = id_pair
  )
  IDs_subject1 <- unique(pair_df[,c("ID_subject1","pair_ID")])[,"ID_subject1"]
  IDs_subject2 <- unique(pair_df[,c("ID_subject2","pair_ID")])[,"ID_subject2"]
  # correlated pairs like (1,2) and (1,3), shared subject is on the same side
  shared.factor=1
  # correlated pairs like (1,2) and (2,4), shared subject is on either side
  switched.factor=1
  # self correlation, (1,2) and (1,2)
  self.factor=1
  # Compute column sums across rows of score matrix-like object for pairs with the same subject1/2
  Usum1.tmp <- rowsum(U,IDs_subject1,reorder=FALSE)
  Usum2.tmp <- rowsum(U,IDs_subject2,reorder=FALSE)
  Usum1  <- matrix(nrow = N, ncol = ncol(U),0)
  Usum2  <- matrix(nrow = N, ncol = ncol(U),0)
  Usum1[unique(IDs_subject1),] <- Usum1.tmp
  Usum2[unique(IDs_subject2),] <- Usum2.tmp
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
    out$linear.predictors <- ans$X %*% ans$beta
  } else {
    out$linear.predictors <- out$offset + ans$X %*% ans$beta
  }
  
  out$fitted.values <- family(out)$linkinv(out$linear.predictors)
  out$modelInfo <- ans$model
  out$id_pair <- ans$id_pair
  out$call <- ans$call
  out$corstr    <- ans$model$corstr
  out$cor.link  <- ans$model$cor.link
  out$control   <- ans$control
  out$std.err   <- "san.se.modified"
  out$var <- varcov
  out$pair_table <- num_visits_per_pair
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
