#' @title Calculate standardized log-likelihood statistic (lZ) for person fit evaluation
#'
#' @description This function calculates the standardized log-likelihood statistic (lZ; Cui & Li, 2015); Drasgow et al. 1985) and the proposals for correcting its distribution discussed in Santos et al., 2019).
#'
#' @param fit A G-DINA model fit object from the \code{GDINA} package (Ma & de la Torre, 2020).
#' @param att.est What attribute estimates are used? The default is \code{"MLE"}.
#' @param sig.level Scalar numeric. Alpha level for decision. Default is 0.05.
#' @param p.adjust.method Scalar character. Correction method for p-values. Possible values include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none". See p.adjust function from the stats R package for additional details. Default is BH.
#'
#' @return \code{personFit} returns an object of class \code{personFit}, with a list of elements:
#' \describe{
#' \item{\code{stat}}{Person fit statistics (\code{data.frame}).}
#' \item{\code{p}}{p-values (two-sided test) for the person fit statistics (\code{data.frame}).}
#' \item{\code{sigp}}{Scalar vectors denoting the examinees for which the person fit statitic is significant (p-value) (\code{list}).}
#' \item{\code{sigadjp}}{Scalar vectors denoting the examinees for which the person fit statitic is significant (adjusted p-value) (\code{list}).}
#' }
#'
#' @author {Pablo Nájera, Universidad Autónoma de Madrid, \email{pablo.najera@uam.es}, \cr Miguel A. Sorrel, Universidad Autónoma de Madrid, \cr Kevin K. Santos, University of the Philippines}
#'
#' @references
#'
#' Cui, Y., & Li, J. (2015). Evaluating person fit for cognitive diagnostic assessment. \emp{Applied Psychological Measurement}, \emp{39}, 223–238. https://doi.org/10.1177/0146621614557272
#'
#' Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness measurement with polychotomous item response models and standardized indices. \emph{British Journal of Mathematical and Statistical Psychology}, \emph{38}, 67–86. https://psycnet.apa.org/doi/10.1111/j.2044-8317.1985.tb00817.x
#'
#' Santos, K. C. P., de la Torre, J., & von Davier, M. (2020). Adjusting person fit index for skewness in cognitive diagnosis modeling. \emp{Journal of Classification}, \emp{37}, 399-420. https://doi.org/10.1007/s00357-019-09325-5
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(GDINA)
#' dat <- sim10GDINA$simdat[1:20, ]
#' Q <- sim10GDINA$simQ
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' res.personFit <- personFit(fit)
#' res.personFit
#' }

personFit <- function(fit, att.est = "MLE", sig.level = 0.05, p.adjust.method = "BH"){
  dat <- fit$options$dat
  Q <- fit$options$Q
  N <- nrow(dat)
  J <- nrow(Q)
  K <- ncol(Q)
  LCprob <- fit$LC.prob
  att <- switch(att.est,
                "EAP" = as.matrix(GDINA::personparm(fit, what = c("EAP"))),
                "MAP" = as.matrix(GDINA::personparm(fit, what = c("MAP"))[,1:K]),
                "MLE" = as.matrix(GDINA::personparm(fit, what = c("MLE"))[,1:K]))
  prob <- t(LCprob[, match(apply(att, 1, paste, collapse = ""), apply(GDINA::attributepattern(K), 1, paste, collapse = ""))])
  
  l0 <- as.matrix(rowSums(dat * log(prob) + (1 - dat) * log(1 - prob)))
  mean.l0 <- rowSums(prob * log(prob) + (1 - prob) * log(1 - prob))
  var.l0 <- rowSums(prob * (1 - prob) * (log(prob/(1 - prob)))^2)
  lz <- (l0 - mean.l0)/sqrt(var.l0)
  skew.lz <- rowSums((1/sqrt(var.l0)^3) * prob * (1 - prob) * (1 - 2 * prob) * (log(prob) - log(1 - prob))^3)
  lz.p <- pnorm(lz)*2
  lz.p[lz.p > 1] <- 1
  
  v <- 8 / (skew.lz^2)
  b <- sqrt(var.l0/(2 * v))
  a <- -b * v - mean.l0
  lz.chi <- abs(l0 + a) / b
  lz.chi.p <- pchisq(lz.chi, df = v, lower.tail = FALSE)
  lz.chi.p[lz.chi.p > 1] <- 1
  
  lz.cf <- lz - skew.lz * (lz^2 - 1)/12
  lz.cf.p <- pnorm(lz.cf)*2
  lz.cf.p[lz.cf.p > 1] <- 1
  
  lz.edge.p <- lz.p - dnorm(lz, mean = 0, sd = 1) * (lz^2 - 1) * (skew.lz/6)
  lz.edge.p[lz.edge.p > 1] <- 1
  
  res.stat <- data.frame("lz" = lz, "lz.chi" = lz.chi, "lz.edge" = NA, "lz.cf" = lz.cf)
  res.p <- data.frame("lz" = lz.p, "lz.chi" = lz.chi.p, "lz.edge" = lz.edge.p, "lz.cf" = lz.cf.p)
  rownames(res.stat) <- rownames(res.p) <- 1:N
  res.adjp <- apply(res.p, 2, function(x) p.adjust(x, method = p.adjust.method))
  res.sigp <- res.p < sig.level
  res.sigadjp <- res.adjp < sig.level
  list.sigp <- apply(res.sigp, 2, function(x) as.numeric(which(x)))
  list.sigadjp <- apply(res.sigadjp, 2, function(x) as.numeric(which(x)))
  res <- list(stat = res.stat, p = res.p, adjp = res.adjp, sigp = list.sigp, sigadjp = list.sigadjp)
  class(res) <- "personFit"
  return(res)
}
