## LARS/LASSO stage 1 
## Without imputation (complete-case analysis)
## Written to rerun ARIES analyses 
## Yiwen Zhu
## Oct 23, 2018

## df: analytic sample with complete data
## outcome.vec: a vector of outcome names (e.g. all CpGs)
## adver: short name of adversity
## covars: a list of covariates
## hypos: additional hypotheses to be tested (aside from the individual time point)
## recency.vec: recency vector - weighting by age in months 
## inf.method: inference method, covTest, sI (selective inference)
## exposures: default to individual time points; otherwise specify a list of collapsed time points
## sample ID column must be named "ID"

## NOTE: Before running this, please check carefully that the variable type in the df 
## are set to the correct type(categorical/numeric)
## otherwise this won't run properly!! 

# ## personal package library path --- 
# x <- .libPaths()
# y <- .libPaths()
# x[1] <- y[2]
# x[2] <- y[1]
# .libPaths(x)




# loadlibs <- function(){
#   library(lars)
#   library(covTest)
#   library(dplyr)
#   library(selectiveInference)
# }
# suppressPackageStartupMessages(loadlibs())
# 
# rm(loadlibs)




select.LARS.complete <- function(df, outcome.vec, adver, covars, hypos, recency.vec = NULL, 
                                inf.method = "covTest", 
                                exposures = "default", 
                                FWL = FALSE) {
  
  DNAm <- df[, c("ID", outcome.vec)]
  #df <- df[, -which(colnames(df) %in% outcome.vec)]
  df <- df[, 1:86]
  
  ## data prep ---- 
  
  # define exposures 
  if (exposures == "default") { 
    exposures <- grep(paste0("^", adver), colnames(df), value = T)
  }
  
  # subset to those with complete data on this exposure and all covariates
  df <- df %>% filter(rowSums(is.na(.[,c(exposures,covars)])) == 0)
  DNAm <- DNAm[DNAm$ID %in% df$ID, ]
  
  
  # define sample size:
  n <- nrow(df) 
  
  # add ever/never exposed, accumulation, and recency 
  exp <- do.call("cbind", lapply(df[, exposures], function(x) as.numeric(as.character(x))))
  df <- df %>% 
    mutate(accumulation = rowSums(exp)) %>% 
    mutate(ever = ifelse(accumulation > 0, 1, 0))
  
  if ("recency" %in% hypos){
    df$recency <-  rowSums(exp %*% diag(recency.vec))
  } 
  
  
  ## select covariates -- 
  ## what are the categorical covariates
  cat.covars <- do.call("c", lapply(1:length(covars), function(x) {
    cov1 <- covars[x]
    if (is.factor(df[, cov1])) {
      return(cov1)
    } else {
      return(NULL)
    }
  }))
  
  ## number of hypotheses 
  n_hypo <- length(exposures) + length(hypos)
  
  ## dummy code the categorical covars
  cat.covariates <- model.matrix(as.formula(paste0("~", paste(cat.covars, collapse = "+"))), df)[,-1]
  covariates <- as.matrix(cbind(cat.covariates, as.matrix(df[,setdiff(covars, cat.covars)])))
  X_hypos <- do.call("cbind", lapply(df[, c(exposures, hypos)], function(x) as.numeric(as.character(x))))
  X_residual <- lm(X_hypos ~ covariates)$resid
  
  ## Normalize X_hypos for use in max-|t| calculations
  col_mean <- apply(X_residual, 2, mean)
  X_centered <- X_residual - rep(col_mean, rep(n, n_hypo)) #subtract mean
  col_sss <- apply(X_centered, 2, function(x) sqrt(sum(x^2)))
  X_normed <- X_centered / rep(col_sss, rep(n, n_hypo)) #divide by sqrt sum squares
  Xt <- t(X_normed)
  XtX <- Xt %*% X_normed
  
  ## create LASSO function to return results
  
  ## covariance test ---- 
  ## y: single locus DNAm 
  run.covTest <- function(i){
    y <- DNAm[, outcome.vec[i]]
    y <- as.numeric(as.character(y))
    
    ## Frish-Waugh-Lovell 
    if (FWL) {
      y <- lm(y~covariates)$resid
    }
    
    if (i %% 1000 == 0){
      print(i)
    }
    
    lasso <- lars(X_residual, y)
    # last.action <- length(lasso$actions)
    
    ## Covariance test results
    # d <- covTest3(lasso, X_residual, y, maxp=last.action)$results
    
    d <- covTest3(lasso, X_residual, y, maxp=1)$results
    d <- as.data.frame(d)
    row.names(d) <- NULL
    d$Variable_Name <- colnames(X_residual)[abs(d$Predictor_Number)]
    d$Improvement_R2 <- lasso$R2[2]
    ## Remove steps where hypo dropped:
    d <- d[!is.na(d$Drop_in_covariance), ]
    
    # ## Make vector of results, up to max # hypos added:
    # k <- nrow(d)
    # d2 <- as.data.frame(t(d[, c(4, 5, 3)]))
    # r <- as.character(unlist(unname(d2)))
    # r2 <- c(r, rep(NA, c((15*3) - length(r))))
    
  }
  
  ## selective inference ---- 
  # create a similar function using selectiveInference: fixedLassoInference
  run.sI <- function(i){
    
    y <- DNAm[, outcome.vec[i]]
    y <- as.numeric(as.character(y))
    
    ## Frish-Waugh-Lovell 
    if (FWL) {
      y <- lm(y~covariates)$resid
    }
    
    if (i %% 1000 == 0){
      print(i)
    }
    
    # prevent p-value calculations from rounding
    options(digits=10)
    
    # Apply LARS
    lasso <- lars(X_residual, y)
    
    # scale X_residual
    sumsq <- lasso$normx
    X_normed <- scale(X_residual, scale = sumsq)
    
    #Covariance test results
    # first significant one
    fli <- fixedLassoInf(X_normed, y, lasso$beta[2, ], lasso$lambda[2], tol.beta = 1e-22) # lowered the tolerance for determining if a coefficient is zero
    
    # Construct a similar data frame
    fli.df <- data.frame(matrix(NA, nrow = 1, ncol = 0))
    fli.df$Predictor_Number <- fli$vars
    fli.df$Beta <- fli$coef0
    fli.df$SE <- fli$sd
    fli.df$`P-value` <- fli$pv
    fli.df$Variable_Name <- colnames(X_residual)[abs(fli.df$Predictor_Number)]
    fli.df$R2 <- lasso$R2[[2]]
    fli.df <- cbind(fli.df, fli$ci)
    
    return(fli.df)
    
  }
  
  ### here add the other two/three methods in ---- 
  
  run.max.t <- function(i){
    y <- DNAm[, outcome.vec[i]]
    y <- as.numeric(as.character(y))
    
    ## Frish-Waugh-Lovell 
    if (FWL) {
      y <- lm(y~covariates)$resid
    }
    
    if (i %% 1000 == 0){
      print(i)
    }
    
    # prevent p-value calculations from rounding
    options(digits=10)
    
    # My calculations
    y_centered <- y-mean(y)
    Xty <- Xt %*% y_centered
    selection <- which.max(abs(Xty))
    absbeta <- abs(Xty[selection])
    s <- summary(lm(y_centered ~ X_normed))$sigma
    set.seed(y*n)
    p.mycalc <- 1 - pmvt(lower=-rep(absbeta,n_hypo),
                         upper= rep(absbeta,n_hypo),
                         delta= rep(0,n_hypo),
                         df= n-n_hypo, sigma= s^2 * XtX)
    
    maxt.df <- data.frame(matrix(NA, nrow = 1, ncol = 0))
    maxt.df$Predictor_Number <- selection
    maxt.df$pval <- p.mycalc
    maxt.df$Variable_Name <- colnames(X_residual)[abs(maxt.df$Predictor_Number)]
    
    return(maxt.df)
  }
  

  run.naive <- function(i){
    y <- DNAm[, outcome.vec[i]]
    y <- as.numeric(as.character(y))
    
    ## Frish-Waugh-Lovell 
    if (FWL) {
      y <- lm(y~covariates)$resid
    }
    
    if (i %% 1000 == 0){
      print(i)
    }
    
    # prevent p-value calculations from rounding
    options(digits=10)
    
    # My calculations
    y_centered <- y-mean(y)
    Xty <- Xt %*% y_centered
    selection <- which.max(abs(Xty))
    
    # Naive calculation
    p.naive <- summary(lm(y ~ X_normed[,selection]))$coef[2,4]
    p.bonf <- p.naive*n_hypo
    p.bonf[p.bonf > 1] <- 1
    
    naive.df <- data.frame(matrix(NA, nrow = 2, ncol = 0))
    naive.df$Predictor_Number <- selection
    naive.df$pval <- c(p.naive, p.bonf)
    naive.df$Variable_Name <- colnames(X_residual)[abs(naive.df$Predictor_Number)]
    naive.df$ptype <- c("naive", "naive.bonf")
    
    
    return(naive.df)
  }
  
  if (inf.method == "covTest") { 
    results <- as.data.frame(t(sapply(1:length(outcome.vec), run.covTest)))
    }
  
  if (inf.method == "sI") { 
    results <- as.data.frame(t(sapply(1:length(outcome.vec), run.sI)))
  }
  
  if (inf.method == "maxt") { 
    results <- as.data.frame(t(sapply(1:length(outcome.vec), run.max.t)))
  }
  
  if (inf.method == "naive") { 
    results <- lapply(1:length(outcome.vec), run.naive) %>% bind_rows()
  }
  
  return(results)
}






###-----------------------------------------------
## load modified covTest, without rounding p-values:
covTest3 <- function (fitobj, x, y, sigma.est = "full", status = NULL, maxp = min(nrow(x), 
                                                                                  ncol(x))) 
{
  s4 = substring(fitobj$call, 1, 4)[1]
  s7 = substring(fitobj$call, 1, 7)[1]
  s8 = substring(fitobj$call, 1, 8)[1]
  if (s4 == "lars" & s7 != "lars.en" & s8 != "lars.glm") {
    calltype = "lars"
    family = "gaussian"
  }
  if (s7 == "lars.en") {
    calltype = "lars.en"
    family = "gaussian"
  }
  if (s8 == "lars.glm") {
    calltype = "lars.glm"
    family = fitobj$family
  }
  if (family == "cox") 
    stop("Cox model not yet implemented")
  if (calltype == "lars") {
    if (fitobj$type != "LASSO") {
      stop("Call to Lars must use type='LASSO'")
    }
  }
  if (calltype == "lars") {
    type = "lar"
  }
  if (calltype == "lars.en") {
    type = "lars.en"
  }
  if (calltype == "lars.glm") {
    type = "lars.glm"
  }
  if (calltype == "lars.glm" & sigma.est == "full") {
    sigma.est = 1
    cat("glm model; sigma set to 1", fill = TRUE)
  }
  n = nrow(x)
  p = ncol(x)
  my = mean(y)
  lambda.min.ratio = ifelse(nrow(x) < ncol(x), 0.1, 1e-04)
  jlist = unlist(fitobj$act)
  if (type == "lar") 
    lamlist = c(fitobj$lambda, 0)
  if (type == "lars.en") 
    lamlist = c(fitobj$lambda, 0)
  if (type == "lars.glm") 
    lamlist = c(fitobj$lambda, 0)
  maxp.call = maxp
  maxp = length(jlist)
  maxp = min(maxp, which(lamlist == 0))
  maxp = min(maxp, maxp.call)
  jlist = jlist[1:maxp]
  cov0 = cov = sig = rep(NA, maxp)
  yy = y - my
  if (family == "binomial") {
    glmobj = glmnet(x, y, family = "binomial", standardize = fitobj$standardize, 
                    lambda.min.ratio = lambda.min.ratio)
  }
  if (family == "cox") {
    glmobj = glmnet(x, Surv(y, status), family = "cox", standardize = fitobj$standardize, 
                    lambda.min.ratio = lambda.min.ratio)
  }
  if (family == "cox") {
    junk = calcz02(x, y, status)
    sc = junk$sceta
    inf = junk$infeta
    miinf = misqrt((t(inf) + inf)/2)
    yy = t(sc) %*% miinf
  }
  for (j in 1:maxp) {
    if (jlist[j] > 0) {
      lambda = lamlist[j + 1]
      if (type == "lar") 
        yhat = predict(fitobj, x, s = lambda, type = "fit", 
                       mode = "lam")$fit
      if (type == "lars.en") 
        yhat = (1 + fitobj$lambda2) * predict.lars.en(fitobj, 
                                                      x, lambda)
      if (type == "lars.glm" & family == "binomial") {
        yhat = as.vector(predict(glmobj, x, type = "link", 
                                 s = lambda/n))
      }
      if (type == "lars.glm" & family == "cox") {
        yhat = as.vector(predict(glmobj, x, type = "link", 
                                 s = lambda/n))
      }
      cov[j] = sum(yy * yhat)
      if (j == 1) {
        cov0[j] = 0
      }
      if (j > 1) {
        tt0 = which(fitobj$beta[j, ] != 0)
        if (type == "lar") {
          aa = update(fitobj, x = x[, tt0, drop = F])
          yhat0 = predict(aa, x[, tt0], type = "fit", 
                          s = lambda, mode = "lam")$fit
        }
        if (type == "lars.en") {
          aa = update(fitobj, x = x[, tt0, drop = F])
          yhat0 = (1 + fitobj$lambda2) * predict.lars.en(aa, 
                                                         x[, tt0], lambda)
        }
        if (type == "lars.glm") {
          if (family == "binomial") {
            if (length(tt0) == 1) {
              tt0 = c(tt0, tt0)
            }
            glmobj0 = glmnet(x[, tt0, drop = F], y, family = "binomial", 
                             standardize = fitobj$standardize, lambda.min.ratio = lambda.min.ratio)
            yhat0 = as.vector(predict(glmobj0, x[, tt0, 
                                                 drop = F], type = "link", s = lambda/n))
          }
          if (family == "cox") {
            if (length(tt0) == 1) {
              tt0 = c(tt0, tt0)
            }
            glmobj0 = glmnet(x[, tt0, drop = F], Surv(y, 
                                                      status), family = "cox", standardize = fitobj$standardize, 
                             lambda.min.ratio = lambda.min.ratio)
            yhat0 = as.vector(predict(glmobj0, x[, tt0, 
                                                 drop = F], type = "link", s = lambda/n))
          }
        }
        cov0[j] = sum(yy * yhat0)
      }
    }
  }
  if (is.numeric((sigma.est))) {
    sigma = sigma.est
    sigma.type = "known"
    null.dist = "Exp(1)"
    if (sigma.est <= 0) {
      stop("sigma.est must be positive")
    }
  }
  if (sigma.est == "full") {
    if (nrow(x) < ncol(x) + 1) 
      stop("Number of observations must exceed number of variables,\nwhen sigma.est is `full; you need to specify a numeric value for sigma.est")
    sigma.type = "unknown"
    aaa = lsfit(x, y)
    sigma = sqrt(sum(aaa$res^2)/(n - p))
    np = n - p
    null.dist = paste("F(2,", as.character(np), ")", sep = "")
  }
  tt = ((cov - cov0)/sigma^2)
  if (sigma.type == "known") {
    out = cbind(jlist, tt, 1 - pexp(tt, 1))
    dimnames(out)[[2]] = c("Predictor_Number", "Drop_in_covariance", 
                           "P-value")
  }
  if (sigma.type == "unknown") {
    out = cbind(jlist, tt, 1 - pf(tt, 2, n - p))
    dimnames(out)[[2]] = c("Predictor_Number", "Drop_in_covariance", 
                           "P-value")
  }
  dimnames(out)[[1]] = rep("", nrow(out))
  return(list(results = out, 
              sigma = round(sigma,  4), 
              null.dist = null.dist))
}
