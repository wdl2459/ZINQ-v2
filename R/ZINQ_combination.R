
#' Combine the marginal p-values
#'
#' @param input An output from \code{\link{ZINQ_tests}}.
#' @param method Combination method, "MinP" for MinP test, "Cauchy" for Cauchy combination test; default is "MinP".
#' @param taus A grid of quantile levels, must be a subset or equal to that from \code{input}; default is c(0.1, 0.25, 0.5, 0.75, 0.9).
#' @param M The number of MC draws from the joint distribution of quantile rank-scores when \code{method} is "MinP"; default is 10000.
#'
#' @details
#' \itemize{
#'   \item Please choose 'MinP' or 'Cauchy' for \code{method}, no other options.
#'   \item \code{taus} must be a subset or equal to the grid used to produce \code{input}.
#' }
#'
#' @return A pvalue, the final p-value of ZINQ
#'
#' @references
#' \itemize{
#'   \item Ling, W. et al. (2020+). Powerful and robust non-parametric association testing for microbiome data via a zero-inflated quantile approach (ZINQ).
#'   \item He, Z. et al. (2017). Unified sequence-based association tests allowing for multiple functional annotations and meta-analysis of noncoding variation in metabochip data. The American Journal of HumanGenetics 101(3), 340–352.
#'   \item Lee, S. et al. (2012). Optimal tests for rare variant effects in sequencing association studies. Biostatistics 13(4), 762–775.
#'   \item Liu, Y., Xie, J. (2019). Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures. Journal of the American Statistical Association, 1–18
#' }
#'
#' @examples
#' library(quantreg)
#' library(MASS)
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
#' result = ZINQ_tests(formula.logistic=y~x, formula.quantile=y~x, C="x", data=dat)
#' ZINQ_combination(result, method="Cauchy")
#'
#' @export

### tests combination ###

ZINQ_combination <- function(input, method="MinP", taus=c(0.1, 0.25, 0.5, 0.75, 0.9), M=10000){

  ## check

  # choose either MinP or Cauchy
  if (method != "MinP" & method != "Cauchy"){
    stop("Please choose 'MinP' or 'Cauchy', no other options.")
  }

  # taus and ind should match
  if (!all(taus %in% input$taus)){
    stop("taus should be a subset of that taus used to produce input.")
  }


  ## whether from single_CorB=T or not
  if (!is.null( ncol(input$Sigma.hat) )){
    if (length(input$pvalue.quantile) != ncol(input$Sigma.hat)){
      single_CorB = F
      df = ncol(input$Sigma.hat) / length(input$pvalue.quantile)
    } else single_CorB = T
  } else single_CorB = T


  ## compute the aggregate p-value
  width = length(taus)

  ind = match(taus, input$taus)
  pvalue.quantile = input$pvalue.quantile[ind]


  if (method == "MinP"){

    # t.obs in the minp test
    t.obs = min(input$pvalue.logistic, pvalue.quantile)

    if (single_CorB != F){

      if (is.null( ncol(input$Sigma.hat) )){
        Sigma.hat = input$Sigma.hat[ind]
      } else Sigma.hat = input$Sigma.hat[ind, ind]

      # the (1 - t.obs/2)th percentile of the statistics in quantile regression, normal
      if (is.null( ncol(input$Sigma.hat) )){
        sigma.hat = sqrt( Sigma.hat )
      } else sigma.hat = sqrt( diag(Sigma.hat) )
      qmin.quantile = qnorm((1-t.obs/2), mean=0, sd=sigma.hat)

      # MC, to estimate the probability that the absolute of joint statistics in quantile regression < each threshold
      beta.sim = mvrnorm(n=M, mu=rep(0, width), Sigma=Sigma.hat)
      prob.quantile = mean( apply(beta.sim, 1, function(z){ all( abs(z) < qmin.quantile ) }) )

    } else {

      index = unlist( lapply(ind, function(kk){ ((kk-1)*df+1):(kk*df) }) )
      Sigma.hat = input$Sigma.hat[index, index]

      # the (1 - t.obs)th percentile of the statistics in quantile regression, chisq
      qmin.quantile = rep(qchisq(1-t.obs, df=df), width)

      # MC, to estimate the probability that the joint statistics in quantile regression < each threshold
      beta.sim = mvrnorm(n=M, mu=rep(0, width*df), Sigma=Sigma.hat)

      prob.quantile = mean( apply(beta.sim, 1, function(z){
        obs.sim = NULL
        for (kk in 1:width){
          obs.sim[kk] = t( z[((kk-1)*df+1):(kk*df)] ) %*% solve(Sigma.hat[((kk-1)*df+1):(kk*df), ((kk-1)*df+1):(kk*df)]) %*% z[((kk-1)*df+1):(kk*df)]
        }
        all( obs.sim < qmin.quantile )
      }) )

    }

    # final p-value
    pvalue = 1 - (1 - t.obs) * prob.quantile

  } else {

    # weights for quantile tests
    w = taus*(taus <= 0.5) + (1-taus)*(taus > 0.5)
    w = w / sum(w) * (1 - input$zerorate)

    stats.cauchy = input$zerorate * tan( (0.5-input$pvalue.logistic)*pi ) + sum ( w * tan( (0.5-pvalue.quantile)*pi ) )

    # final p-value
    pvalue = 1 - pcauchy(stats.cauchy)

  }

  return(pvalue)

}

