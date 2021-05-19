
#' Sanity check before applying ZINQ
#'
#' @param tax_tab The taxa read count table (un-normalized), sample (row) by taxa (col).
#' @param metadata The metadata, sample (row) by variable (col).
#' @param C The name(s) of clinical variable(s) of interest, e.g., "Condition" or c("Condition", "Batch").
#'
#' @details
#' \itemize{
#'   \item It is recommended to do the sanity check before applying ZINQ. If it is necessary, warnings will be printed to guide the analysis using ZINQ.
#'   \item If library size is a confounder of the variable(s) of interest, ZINQ might not control type I error.
#'   \item If there are few non-zero read counts, use ZINQ with caution.
#'   \item ZINQ is not designed for perfect separation, e.g., there are all zeroes in one group (case or control).
#'   \item The sanity check is mainly about zero inflation. Most normalizations will keep the original zeroes, thus investigaing the un-normalized taxa read count table provides sufficient clues to use ZINQ.
#'    For normalizations not retaining the zeroes. e.g., CLR, results of the sanity check is not informative, one can apply ZINQ directly.
#' }
#'
#' @return Print warnings if necessary
#' \itemize{
#'   \item When library size is a confounder
#'   \item For each taxon, (1) when all read counts are zero, (2) when there are limited non-zero read counts (<30 or <15), (3) when there is a perfect separation w.r.t. the variable(s) of interest.
#' }
#'
#' @export

ZINQ_check <- function(tax_tab, metadata, C){

  # whether library size is a major confounder of components of C
  libsize = apply(tax_tab, 1, sum)
  C_expanded = model.frame(paste("~", C, collapse = "+"), data=metadata)

  p = NULL
  for (ii in 1:ncol(C_expanded)){
    if (any( (C_expanded[, ii] %in% c(0,1)) == F )){ # continuous variable
      p[ii] = cor.test(libsize, C_expanded[, ii])$p.value
    } else{ # expanded from a discrete variable
      p[ii] = t.test(libsize[C_expanded[, ii]==0], libsize[C_expanded[, ii]==1])$p.value
    }
  }

  if (any(p < 0.05)) warning("Library size is a confounder of the variable(s) of interest. ZINQ might not control type I error.
                             Use it with caution, or use other methods that treat zeroes and non-zeroes quantitatively, such as LDM.")


  ## check every taxon
  tax_names = colnames(tax_tab)
  if (is.null(tax_names) | length(unique(tax_names)) < length(tax_names)) tax_names = paste0("taxon", 1:ncol(tax_tab))

  for (ii in 1:ncol(tax_tab)){

    b = 1*(tax_tab[, ii] != 0)

    good_status = rep(T, 4)

    # whether there are non-zero outcome
    if (length(which(b == 1)) == 0) good_status[1] = F

    # whether num of non-zero outcome < 30
    if (length(which(b == 1)) < 30 & length(which(b == 1)) >= 15) good_status[2] = F

    # whether num of non-zero outcome < 15
    if (length(which(b == 1)) < 15) good_status[3] = F

    # whether the C of interest has non-zero variance in the positive subset
    test_dat = cbind(b, C_expanded)
    test_dat_nonzero = test_dat[test_dat$b==1, ]
    var_C_nonzero = apply(test_dat_nonzero[, -1, drop=F], 2, sd, na.rm=T)
    if (all(var_C_nonzero == 0)) good_status[4] = F

    if (any(good_status == F)){
      if (good_status[1] == F) warning(paste(tax_names[ii], ": There are no non-zero read counts. ZINQ will return singularity error."))
      if (good_status[2] == F) warning(paste(tax_names[ii], ": The number of non-zero read counts is less than 30. Central quantiles levels, such as taus=0.2,0.4,0.5,0.6,0.8, quartiles (or even median), and Cauchy combination test are recommended."))
      if (good_status[3] == F) warning(paste(tax_names[ii], ": The number of non-zero read counts is less than 15. ZINQ might not control type I error. Use ZINQ with caution (central quantiles levels, such as taus=quartiles (or even median), and Cauchy combination test are recommended), or use other methods, such as LDM."))
      if (good_status[4] == F) warning(paste(tax_names[ii], ": Only a single status of the variable(s) of interest shown in the non-zero subset. Since ZINQ is not designed for perfect separation, NA will be returned by ZINQ. Use other methods, such as LDM."))
    }

  }

}

