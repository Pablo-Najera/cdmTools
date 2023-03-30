#' @title Calculate corrected classification accuracy with multiple imputation
#'
#' @description This function calculates the test-, pattern-, and attribute-level classification accuracy indices based on integrated posterior probabilities from multiple imputed item parameters (Kreitchmann et al., 2022).
#' The classification accuracy indices are the ones developed by Iaconangelo (2017) and Wang et al. (2015).
#' It is only applicable to dichotomous attributes. The function is built upon the \code{CA} function from the \code{GDINA} package (Ma & de la Torre, 2020).
#'
#' @param fit A G-DINA model fit object from the \code{GDINA} package (Ma & de la Torre, 2020).
#' @param what What attribute estimates are used? The default is \code{"EAP"}.
#' @param R Number of bootstrap samples and imputations. The default is 500.
#' @param n.cores Number of processors to use to speed up multiple imputation. The default is 2.
#' @param verbose Show progress. The default is \code{TRUE}.
#' @param seed A seed for obtaining consistent results. If \code{NULL}, no seed is used. The default is \code{NULL}.
#'
#' @return \code{CA.MI} returns an object of class \code{CA}, with a list of elements:
#' \describe{
#' \item{\code{tau}}{Estimated test-level classification accuracy, see Iaconangelo (2017, Eq 2.2) (\code{vector}).}
#' \item{\code{tau_l}}{Estimated pattern-level classification accuracy, see Iaconangelo (2017, p. 13) (\code{vector}).}
#' \item{\code{tau_k}}{Estimated attribute-level classification accuracy, see Wang, et al (2015, p. 461 Eq 6) (\code{vector}).}
#' \item{\code{CCM}}{Conditional classification matrix, see Iaconangelo (2017, p. 13) (\code{matrix}).}
#' }
#'
#' @author {Rodrigo S. Kreitchmann, Universidad Autónoma de Madrid}
#'
#' @references
#' Iaconangelo, C.(2017). \emph{Uses of classification error probabilities in the three-step approach to estimating cognitive diagnosis models}. (Unpublished doctoral dissertation). New Brunswick, NJ: Rutgers University.
#'
#' Kreitchmann, R. S., de la Torre, J., Sorrel, M. A., Nájera, P., & Abad, F. J. (2022). Improving reliability estimation in cognitive diagnosis modeling. \emph{Behavior Research Methods}. https://doi.org/10.3758/s13428-022-01967-5
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R package for cognitive diagnosis modeling. \emph{Journal of Statistical Software}, \emph{93}(14). https://doi.org/10.18637/jss.v093.i14
#'
#' Wang, W., Song, L., Chen, P., Meng, Y., & Ding, S. (2015). Attribute-level and pattern-level classification consistency and accuracy indices for cognitive diagnostic assessment. \emph{Journal of Educational Measurement}, emph{52} , 457-476.
#'
#' @import foreach
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(GDINA)
#' dat <- sim10GDINA$simdat[1:100,]
#' Q <- sim10GDINA$simQ
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' ca.mi <- CA.MI(fit)
#' ca.mi
#' }

CA.MI <- function(fit, what = "EAP", R = 500, n.cores = 1, verbose = TRUE, seed = NULL){

  if(!inherits(fit, "GDINA")){stop("Error in CA.MI: fit must be of class 'GDINA'.")}
  if(!what %in% c("EAP", "MAP", "MLE")){stop("Error in CA.MI: what must be either 'EAP', 'MAP', or 'MLE'.")}
  if(any(apply(fit$options$dat, 2, stats::sd, na.rm = TRUE) == 0)){stop("Error in CA.MI: The data must not contain constant responses for an item.")}
  if((!is.numeric(R) & !is.double(R)) | length(R) > 1){stop("Error in CA.MI: R must be a unique numeric value.")}
  if(R < 1){stop("Error in CA.MI: R must be greater than 0.")}
  if((!is.numeric(n.cores) & !is.double(n.cores)) | length(n.cores) > 1){stop("Error in CA.MI: n.cores must be a unique numeric value.")}
  if(n.cores < 1){stop("Error in CA.MI: n.cores must be greater than 0.")}
  if(!is.logical(verbose)){stop("Error in CA.MI: verbose must be logical.")}
  if(!is.null(seed)){if((!is.numeric(seed) & !is.double(seed)) | length(seed) > 1){stop("Error in CA.MI: seed must be a unique numeric value.")}}
  if(GDINA::extract(fit, "ngroup") != 1) {stop("Error in CA.MI: only applicable for single group analysis.")}

  if(is.null(seed)){seed <- sample(1:100000, 1)}
  set.seed(seed)

  GDINA.options <- formals(GDINA::GDINA)
  GDINA.options <- GDINA.options[-c(1, 2, length(GDINA.options))]
  tmp <- as.list(fit$extra$call)[-c(1:3)]
  GDINA.options[names(GDINA.options) %in% names(tmp)] <- tmp
  GDINA.options$verbose <- 0
  GDINA.options$att.dist <- "fixed"
  GDINA.options$model <- fit$model
  GDINA.options$control <- list(maxitr = rep(0, nrow(fit$options$Q)))

  boot.sampling.dist <- bootSE.parallel(fit, bootsample = R, n.cores = n.cores, verbose = verbose, seed = seed)

  posterior.R <- array(NA, dim = c(dim(t(fit$technicals$logposterior.i)), R))
  if(verbose){cat("\n", "\n")}
  for(r in 1:R){
    if(verbose & r %% 100 == 0){cat("Imputing Item Parameters: Iteration", r, "of", R, "\r")}
    GDINA.options$catprob.parm <- boot.sampling.dist$boot.est$itemprob[[r]]
    GDINA.options$att.prior <- boot.sampling.dist$boot.est$lambda[[r]]
    fit.tmp <- do.call(GDINA::GDINA, c(list(dat = fit$options$dat, fit$options$Q), GDINA.options))
    posterior.R[,,r] <- as.matrix(exp(t(fit.tmp$technicals$logposterior.i)))
  }

  posterior <- apply(posterior.R, 1:2, mean)
  fit.MI <- fit
  fit.MI$technicals$logposterior.i <- t(log(posterior))

  p_c <- GDINA::extract(fit.MI, "posterior.prob")
  pp <- GDINA::personparm(fit, what = what)
  if (what == "MAP" || what == "MLE") {
    if (any(pp[, ncol(pp)]))
      warning(paste0(what, " estimates for some individuals have multiple modes.",
                     collapse = ""), call. = FALSE)
    pp <- as.matrix(pp[, -ncol(pp)])
  }
  mp <- GDINA::personparm(fit.MI, what = "mp")
  patt <- GDINA::extract(fit.MI, "attributepattern")
  gr <- cdmTools.matchMatrix(patt, pp)
  pseudo.gr <- setdiff(seq(nrow(patt)), unique(gr))
  gr <- c(gr, pseudo.gr)
  lab <- apply(patt, 1, paste0, collapse = "")
  post <- cbind(exp(t(GDINA::indlogPost(fit.MI))), matrix(0, nrow(patt),
                                                           length(pseudo.gr)))
  CCM <- cdmTools.aggregateCol(post, gr)/c(GDINA::extract(fit.MI, "nobs") *
                                            p_c)
  tau_c <- diag(CCM)
  tau <- sum(tau_c * c(p_c))
  tau_k <- colMeans(pp * mp + (1 - pp) * (1 - mp))
  names(tau_c) <- rownames(CCM) <- colnames(CCM) <- lab
  res <- list(tau = tau, tau_l = tau_c, tau_k = tau_k, CCM = CCM)
  class(res) <- "CA"

  return(res)
}
