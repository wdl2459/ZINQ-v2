
#' Marginal tests for the logistic and quantile regression components
#'
#' @import quantreg
#' @import MASS
#'
#' @param formula.logistic The full model of logistic regression, e.g., Y ~ X + Y + Z, where Y is zero-inflated.
#' @param formula.quantile The full model of quantile regression, can be different from formula.logistic.
#' @param C The name(s) of clinical variable(s) of interest, e.g., "Condition" or c("Condition", "Batch").
#' @param y_CorD An indicator: use "D" if Y is count, a perturbation from U(0, 1) will be added to the response; use "C" if Y is continuous; default is "C".
#' @param data A data.frame: better cleaned and processed, use numeric for Y and binary covariates, use factor for multi-class discrete covariates.
#' @param taus A grid of quantile levels, e.g., 0.5 for the median, 0.75 for the 3rd quartile; default is c(0.1, 0.25, 0.5, 0.75, 0.9).
#' @param seed A seed for perturbation when \code{y_CorD} is "D"; default is 2020.
#'
#' @details
#' \itemize{
#'   \item Compositional data is regarded as continuous, determined by its support.
#'   \item \code{taus} is a tuning parameter that does not have an efficient selection process yet, try from coarsed to fine grids (e.g., seq(0.1, 0.9, by=0.2) to seq(0.1, 0.9, by=0.1)),
#'           or try adding more extreme levels (e.g., c(0.25, 0.5, 0.75) to c(0.1, 0.25, 0.5, 0.75, 0.9)), with a goal to keep type I error controlled and boost the power;
#'         for common taxa, start from the default; for rare taxa, start from c(0.25, 0.5, 0.75).
#'   \item Quantile rank-score test corrected for zero-inflation is used for the quantile regression component.
#'   \item If \code{C} is a single continuous or binary covariate, Wald test is used for the logistic regression component, else Rao's score test is used
#' }
#'
#' @return A list
#' \itemize{
#'   \item pvalue.logistic - A single p-value from the logistic regression component.
#'   \item pvalue.quantile - A length(\code{taus}) by 1 vector, a squence of p-values from the quantile regression component.
#'   \item Sigma.hat - A df x length(\code{taus}) by df x length(\code{taus}) matrix, where df is the dimension of \code{C}, the covariance matrix of quantile rank-scores.
#'   \item zerorate - The proportion of zeroes in Y.
#'   \item taus - The grid of quantile levels used.
#' }
#'
#' @references
#' \itemize{
#'   \item Ling, W. et al. (2020+). Powerful and robust non-parametric association testing for microbiome data via a zero-inflated quantile approach (ZINQ)
#'   \item Machado, J.A.F., Silva, J.S. (2005). Quantiles for counts. Journal of the American Statistical Association 100(472), 1226–1237.
#' }
#'
#'
#' @examples
#' library(quantreg)
#' n = 300
#' p <- function(x0, gam0=0.75, gam1=-0.15){
#'   lc = gam0 + gam1*x0
#'   exp(lc) / (1 + exp(lc))
#' }
#' x = c(rep(0, n), rep(1, n))
#' w = 0.5 + 1.5*x + (1+0.15*x)*rchisq(2*n,df=1)
#' b = rbinom(2*n, 1, p(x))
#' y = w*b
#' dat = data.frame(y, x)
#'
#' ZINQ_tests(formula.logistic=y~x, formula.quantile=y~x, C="x", data=dat)
#'
#' @export

library(logistf)

### marginal tests ###

ZINQ_tests <- function(formula.logistic, formula.quantile, C, y_CorD="C", data, taus=c(0.1, 0.25, 0.5, 0.75, 0.9), seed=2020){

  ## formulas

  # arrange logistic model
  mf.logistic = model.frame(formula.logistic, data=data)
  y = model.response(mf.logistic, "numeric")
  b = 1*(y != 0) ## v2: +ve -> non-zero
  formula.logistic = update(formula.logistic, b ~ .)
  data.logistic = cbind(data, b)
  mf.logistic = model.frame(formula.logistic, data=data.logistic)

  # locate C in logistic model
  namesx = all.vars(formula.logistic)[-1]
  condition.loc = which(namesx %in% C)

  # arrange logistic null model if C is not a single covariate
  if (length(condition.loc) > 1){
    mul.logistic = T
    namesx.null = setdiff(namesx, C)
    if (length(namesx.null) == 0){
      formula.logistic.null = as.formula( "b ~ 1" )
    } else formula.logistic.null = as.formula( paste( "b ~", paste(namesx.null, collapse = "+") ) )
    mf.logistic.null = model.frame(formula.logistic.null, data=data.logistic)
  } else mul.logistic = F


  # quantile model

  # determine elements in quantile model, create the "positive subset"
  namey = all.vars(formula.quantile)[1]
  namesx = all.vars(formula.quantile)[-1]
  namesx.score = setdiff(namesx, C)
  data.quantile = data[b==1, ]
  if (y_CorD == "D"){ # perturbation if response is count
    set.seed(seed)
    data.quantile[, namey] = dither(data.quantile[, namey], type = "right", value = 1)
  }

  # extract C
  formula.quantile = as.formula( paste( namey, "~", paste(C, collapse = "+") ) )
  mf.quantile = model.frame(formula.quantile, data=data.quantile)
  c = model.matrix(attr(mf.quantile, "terms"), data=mf.quantile)[, -1]
  if (is.null(dim(c))){ # determine whether C is a single covariate, also continuous or binary
    single_CorB=T
  } else single_CorB=F

  # arrange quantile null model, and extract Z
  if (length(namesx.score) == 0){
    formula.quantile = as.formula( paste( namey, "~ 1" ) )
  } else formula.quantile = as.formula( paste( namey, "~", paste(namesx.score, collapse = "+") ) )
  mf.quantile = model.frame(formula.quantile, data=data.quantile)
  z = model.matrix(attr(mf.quantile, "terms"), data=mf.quantile)


  ## set up parameters
  m = length(y) # total sample size
  width = length(taus) # size of tau
  zerorate = length(which(b == 0)) / m # rate of 0's


  ## compute p-values from the marginal tests

  if (single_CorB == T){ # when C is a single covariate, either continuous or binary

    # # logistic, wald test
    # mod.logistic = glm(mf.logistic, family=binomial(link = 'logit'))
    # pvalue.logistic = summary(mod.logistic)$coef[condition.loc+1, 4]

    # firth logistic, profile penalized log-likelihood test
    mod.logistic = logistf(mf.logistic)
    pvalue.logistic = mod.logistic$prob[condition.loc+1]
    
    # estimate quantiles of y|y>0 | H0
    rq0 = rq(mf.quantile, tau=taus)
    qpred0 = predict(rq0)

    # project C on the space of intercept and Z
    C.star = c - z %*% solve( (t(z) %*% z) ) %*% t(z) %*% c

    # compute the rank-score test stats, and its covariance matrix
    RS = unlist( lapply(1:width, function(kk){ sum( (taus[kk] - (data.quantile[, namey] < as.matrix(qpred0, ncol=width)[, kk]))*C.star ) / sqrt(m) }) ) ## v2: as.matrix -> matrix

    if (width == 1){
      cov.RS = taus*(1 - taus)
    } else {
      cov.RS = matrix(0, ncol=width, nrow=width)
      for (kk in 1:(width-1)){
        for (ll in (kk+1):width){
          cov.RS[kk, ll] = min(taus[kk], taus[ll]) - taus[kk]*taus[ll]
        }
      }
      cov.RS = cov.RS + t(cov.RS) + diag(taus*(1 - taus))
    }

    Sigma.hat = cov.RS * sum( C.star^2 ) / m
    if (width == 1){
      sigma.hat = sqrt( Sigma.hat )
    } else {
      sigma.hat = sqrt( diag(Sigma.hat) )
    }

    # marginal p-value in quantile regression
    pvalue.quantile = 2*( 1 - pnorm( abs( RS / sigma.hat ) ) )

  } else { # when C is a set of covariates, or a single covariate with multiple categories

    # # logistic, score test
    # if (mul.logistic != T){
    #   mod.logistic = glm(mf.logistic, family=binomial(link = 'logit'))
    #   pvalue.logistic = anova(mod.logistic, test="Rao")$`Pr(>Chi)`[condition.loc+1]
    # } else {
    #   mod.logistic = glm(mf.logistic, family=binomial(link = 'logit'))
    #   mod.logistic.null = glm(mf.logistic.null, family=binomial(link = 'logit'))
    #   pvalue.logistic = anova(mod.logistic.null, mod.logistic, test="Rao")$`Pr(>Chi)`[2]
    # }
    
    # firth logistic, profile penalized log-likelihood test
    if (mul.logistic != T){
      mod.logistic = logistf(mf.logistic)
      pvalue.logistic = anova(mod.logistic, test="Rao")$`Pr(>Chi)`[condition.loc+1]
    } else {
      mod.logistic = logistf(mf.logistic)
      mod.logistic.null = logistf(mf.logistic.null)
      pvalue.logistic = anova(mod.logistic, mod.logistic.null)$pval
    }

    # estimate quantiles of y|y>0 | H0
    rq0 = rq(mf.quantile, tau=taus)
    qpred0 = predict(rq0)

    # project C on the space of intercept and Z
    C.star = c - z %*% solve( (t(z) %*% z) ) %*% t(z) %*% c

    # compute the rank-score test stats, and its covariance matrix
    RS = lapply(1:width, function(kk){ apply( (taus[kk] - (data.quantile[, namey] < as.matrix(qpred0, ncol=width)[, kk]))*C.star, 2, sum ) / sqrt(m) })

    df = ncol(C.star)
    tmp = t(C.star) %*% C.star / m

    var.RS = NULL
    for (kk in 1:width){
      var.RS[[kk]] = taus[kk] * (1 - taus[kk]) * tmp
    }

    Sigma.hat  = NULL
    for (kk in 1:width){
      temp = NULL
      for (ll in 1:width){
        temp = cbind( temp, ( min(taus[kk], taus[ll]) - taus[kk]*taus[ll] )*tmp )
      }
      Sigma.hat = rbind(Sigma.hat, temp)
    }

    # marginal p-value in quantile regression
    pvalue.quantile = NULL
    for (kk in 1:width){
      stat = t( RS[[kk]] ) %*% solve(var.RS[[kk]]) %*% RS[[kk]]
      pvalue.quantile[kk] = 1 - pchisq(stat, df=df)
    }

  }

  return(list(pvalue.logistic=pvalue.logistic, pvalue.quantile=pvalue.quantile, Sigma.hat=Sigma.hat, zerorate=zerorate, taus=taus))

}

