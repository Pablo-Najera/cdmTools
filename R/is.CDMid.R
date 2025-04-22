#' Check whether a CDM is identified
#'
#' @description Uses a post-hoc simulation approach to check whether a cognitive diagnosis model is identified (i.e., all latent classes are distinguishable; de la Torre et al., 2023).
#'
#' @param fit An object of class RDINA or GDINA (Ma & de la Torre, 2020).
#' @param N A \emph{numeric} value that indicates the number of respondents to simulate. Default is 10000.
#' @param timesJ A \emph{numeric} value that indicates the number of times the test length is multiplied. Default is 20.
#' @param Wald A \code{logical} value that indicates whether the Wald method should be used to find the best model for each item (only applicable if fit is of class \code{GDINA}). Default is \code{FALSE}.
#' @param verbose A \code{logical} value that indicates whether information about the process should be printed or not. Default is \code{TRUE}.
#'
#' @return \code{is.CDMid} returns an object of class \code{is.CDMid}.
#' \describe{
#' \item{\code{total}}{Overall classification accuracy (CCA) and number of posterior multiple modes (PMM). A CCA = 1 indicates that all latent classes are identified (\code{vector}).}
#' \item{\code{class}}{Classification accuracy (CCA) and number of posterior multiple modes (PMM) for each latent class. A CCA = 1 indicates that the latent class is identified (\code{data.frame}).}
#' }
#'
#' @author {Pablo Nájera, Universidad Pontificia Comillas}
#'
#' @references
#' de la Torre, J., Sorrel, M. A., & Nájera, P. (2023, July). \emph{Cognitive diagnosis modeling}. Workshop at the VII International Psychometric Summer School "Applied Psychometrics in Psychology and Education". Yerevan, Armenia.
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(GDINA)
#' dat <- sim30GDINA$simdat
#' Q <- sim30GDINA$simQ
#' fit <- GDINA(dat, Q)
#' id <- is.CDMid(fit)
#' }
is.CDMid <- function(fit, N = 10000, timesJ = 20, Wald = FALSE, verbose = TRUE){

  # Argument control
  if(!inherits(fit, "GDINA")){
    if(Wald){
      Wald <- FALSE
      warning("Warning in CA.MI: Wald requires all items to be fitted by the G-DINA model. Wald turned to FALSE.")
    }
    if(inherits(fit, "RDINA")){
      fit <- RDINA2GDINA(fit)
    } else {
      stop("Error in CA.MI: fit must be of class 'RDINA' or 'GDINA'.")
    }
  } else {
    if(!all(fit$model == "GDINA") & Wald){
      Wald <- FALSE
      warning("Warning in CA.MI: Wald requires all items to be fitted by the G-DINA model. Wald turned to FALSE.")
    }
  }
  if(!is.numeric(N)){stop("Error in CA.MI: N must be numeric.")}
  if(!is.numeric(timesJ)){stop("Error in CA.MI: timesJ must be numeric.")}
  if(!is.logical(Wald)){stop("Error in CA.MI: Wald must be logical.")}
  if(!is.logical(verbose)){stop("Error in CA.MI: verbose must be logical.")}

  # Gather information from fitted model
  dat <- fit$options$dat
  Q <- fit$options$Q
  J <- nrow(Q)
  K <- ncol(Q)
  model.j <- fit$options$model

  # Select most appropriate model for each item with Wald method, and estimate item parameters with those models
  if(Wald){
    if(verbose){cat("\r", "1/6 | Applying Wald method to select the best model for each item...")}
    model.j <- GDINA::modelcomp(fit)$selected.model$models
  } else {
    if(verbose){cat("\r", "1/6 | Model comparison via Wald method has been disabled...")}
  }
  if(Wald){
    if(verbose){cat("\n", "2/6 | Estimating item and person parameters using the most appropriate model for each item...")}
  } else {
    if(verbose){cat("\n", "2/6 | Estimating item and person parameters using the fitted model...")}
  }
  fit <- GDINA::GDINA(dat, Q, model = model.j, verbose = 0)

  # Extract information from new fitted object and generate attribute profiles
  post <- as.numeric(fit$posterior.prob)
  parm <- fit$catprob.parm
  mi.att <- GDINA::attributepattern(K)[sample(1:2^K, size = N, replace = TRUE, prob = post),]
  mi.Q <- do.call("rbind", replicate(timesJ, Q, simplify = FALSE))
  mi.parm <- rep(fit$catprob.parm, timesJ)
  names(mi.parm) <- paste("Item", 1:length(mi.parm))

  # Simulate data using estimated item parameters and generated attribute profiles
  if(verbose){cat(paste("\n", "3/6 | Simulating N =", N, "responses to J =", nrow(mi.Q), "items...                                                     "))}
  mi.sim <- GDINA::simGDINA(N = N, Q = mi.Q, catprob.parm = mi.parm, attribute = mi.att)

  # Extract person parameters using the generating item parameters
  if(verbose){cat("\n", "4/6 | Computing attribute profiles based on the generating item parameters...                ")}
  mi.fit <- GDINA::GDINA(mi.sim$dat, mi.Q, catprob.parm = mi.sim$catprob.parm, control = list(maxitr = 0))

  # Compute conditional classification accuracy
  if(verbose){cat("\n", "5/6 | Computing conditional classification accuracy...                                       ")}
  totalRES <- c(CCA = mean(apply(GDINA::personparm(mi.fit, what = "MLE")[,1:K] == mi.att, 1, all)),
                PMM = mean(GDINA::personparm(mi.fit, what = "MLE")[,K+1]))
  classCCA <- stats::aggregate(apply(GDINA::personparm(mi.fit, what = "MLE")[,1:K] == mi.att, 1, all), by = list(apply(mi.att, 1, paste, collapse = "")), mean)
  classCCA <- classCCA[match(apply(GDINA::attributepattern(K), 1, paste, collapse = ""), classCCA[,1]),]
  classCCA[,1][is.na(classCCA[,1])] <- apply(GDINA::attributepattern(K), 1, paste, collapse = "")[is.na(classCCA[,1])]
  rownames(classCCA) <- 1:nrow(classCCA)
  names(classCCA) <- c("class", "CCA")
  classPMM <- stats::aggregate(GDINA::personparm(mi.fit, what = "MLE")[,K+1], by = list(apply(mi.att, 1, paste, collapse = "")), mean)
  classPMM <- classPMM[match(apply(GDINA::attributepattern(K), 1, paste, collapse = ""), classPMM[,1]),]
  classRES <- cbind(classCCA, PMM = classPMM[,2])

  # Return results
  if(verbose){
    if(all(totalRES["CCA"] == 1)){
      cat("\n", "6/6 | Done, and the model is identified!                                                                 ", "\n")
    } else {
      cat("\n", "6/6 | Done, but there are some identifiability issues...                                                 ", "\n")
    }
  }
  res <- list(total = totalRES, class = classRES)
  class(res) <- "is.CDMid"
  return(res)
}
