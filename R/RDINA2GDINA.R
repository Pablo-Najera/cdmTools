#' Translate RDINA object into GDINA object
#'
#' @description This function translates an object of class \code{RDINA} to an object of class \code{GDINA}, so that the estimated R-DINA object is compatible with most of the functions in the \code{GDINA} package (Ma & de la Torre, 2020), including model fit, item fit, and Q-matrix validation.
#'
#' @param fit An object of class \code{RDINA}.
#'
#' @return \code{RDINA2GDINA} returns an object of class \code{GDINA}. See the \code{GDINA} package for more information.
#'
#' @author {Pablo Nájera, Universidad Pontificia Comillas}
#'
#' @references
#' Ma, W., & de la Torre, J. (2020). GDINA: An R package for cognitive diagnosis modeling. \emph{Journal of Statistical Software}, \emph{93}(14). https://doi.org/10.18637/jss.v093.i14
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(GDINA)
#' dat <- sim30DINA$simdat
#' Q <- sim30DINA$simQ
#' fit1 <- RDINA(dat, Q)
#' fit2 <- RDINA2GDINA(fit1)
#' modelfit(fit2) # Model fit evaluation
#' itemfit(fit2) # Item fit evaluation
#' }
RDINA2GDINA <- function(fit){
  if(!inherits(fit, "RDINA")){stop("Error in RDINA2GDINA: fit must be of class 'RDINA'.")}
  dat <- fit$specifications$dat
  Q <- fit$specifications$Q
  N <- nrow(dat)
  J <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K
  phi <- fit$phi
  catprob.parm <- list()
  for(j in 1:J){
    lj <- 2^sum(Q[j,])
    if(fit$specifications$gate == "AND"){
      catprob.parm[[j]] <- c(rep(phi, lj - 1), 1 - phi)
    } else {
      catprob.parm[[j]] <- c(phi, rep(1 - phi, lj - 1))
    }
  }
  model <- ifelse(fit$specifications$gate == "AND", "DINA", "DINO")
  res <- GDINA::GDINA(dat, Q, model = model, catprob.parm = catprob.parm, control = list(maxitr = 0), verbose = 0)
  res$struc.parm <- as.numeric(fit$post.probs$lp)
  res$testfit[3] <- 1
  res$testfit[-3] <- fit$test.fit
  res$technicals$free.item.npar <- 1
  res$technicals$total.npar <- L
  res$extra$call <- "RDINA"
  return(res)
}
