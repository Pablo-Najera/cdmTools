#' CDM fit comparison - dimensionality assessment method
#'
#' @description A procedure for determining the number of attributes underlying CDM using model fit comparison.
#' For each number of attributes under exploration, a Q-matrix is estimated from the data using the \emph{discrete factor loading} method (Wang, Song, & Ding, 2018), which can be further validated using the \emph{Hull} method (Nájera, Sorrel, de la Torre, & Abad, 2020).
#' Then, a CDM is fitted to the data using the resulting Q-matrix, and several fit indices are computed.
#' After the desired range of number of attributes has been explored, the fit indices are compared.
#' A suggested number of attributes is given for each fit index.
#' The AIC index should be preferred among the other fit indices.
#' For further details, see Nájera, Abad, & Sorrel (2020).
#' This function can be also used by directly providing different Q-matrices (instead of estimating them from the data) in order to compare their fit and select the most appropriate Q-matrix.
#' Note that, if Q-matrices are provided, this function will no longer serve as a dimensionality assessment method, but just as an automated model comparison procedure.
#'
#' @param dat A \emph{N} individuals x \emph{J} items (\code{matrix} or \code{data.frame}). Missing values need to be coded as \code{NA}.
#' @param exploreK Number of attributes to explore. The default is from 1 to 7 attributes.
#' @param Qs A list of Q-matrices to compare in terms of fit. If \code{Qs} is used, \code{exploreK} is ignored.
#' @param stop A fit index to use for stopping the procedure if a model leads to worse fit than a simpler one. This can be useful for saving time without exploring the whole exploreK when it is probable that the correct dimensionality has been already visited. It includes \code{"AIC"}, \code{"BIC"}, \code{"CAIC"}, \code{"SABIC"}, \code{"M2"}, \code{"SRMSR"}, \code{"RMSEA2"}, or \code{"sig.item.pairs"}. The latter represents the number of items that show bad fit with at least another item based on the transformed correlations (see \code{itemfit} function in the \code{GDINA} package; Ma & de la Torre, 2020). It can be also \code{"none"}, which means that the whole \code{exploreK} will be examined. The default is \code{"none"}.
#' @param val.Q Validate the estimated Q-matrices using the \emph{Hull} method? Note that validating the Q-matrix is expected to increase its quality, but the computation time will increase. The default is \code{TRUE}.
#' @param estQ.args A list of arguments for the \emph{discrete factor loading} empirical Q-matrix estimation method (see the \code{estQ} function):
#' \describe{
#' \item{\code{criterion}}{Dichotomization criterion to transform the factor loading matrix into the Q-matrix. The possible options include \code{"row"} (for row means), \code{"col"} (for column means), \code{"loaddiff"} (for the procedure based on loading differences), or a value between 0 and 1 (for a specific threshold). The default is \code{"row"}.}
#' \item{\code{cor}}{Type of correlations to use. It includes \code{"cor"} (for Pearson correlations) and \code{"tet"} (for tetrachoric/polychoric correlations), among others. See \code{fa} function from the \code{psych} R package for additional details. The default is \code{"tet"}.}
#' \item{\code{rotation}}{Rotation procedure to use. It includes \code{"oblimin"}, \code{"varimax"}, and \code{"promax"}, among others. An oblique rotation procedure is usually recommended. See \code{fa} function from the \code{psych} R package for additional details. The default is \code{"oblimin"}.}
#' \item{\code{fm}}{Factoring method to use. It includes \code{"uls"} (for unweighted least squares), \code{"ml"} (for maximum likelihood), and \code{"wls"} (for weighted least squares), among others. See \code{fa} function from the \code{psych} R package for additional details. The default is \code{"uls"}.}
#' }
#' @param valQ.args A list of arguments for the \emph{Hull} empirical Q-matrix validation method. Only applicable if \code{valQ = TRUE} (see the \code{valQ} function):
#' \describe{
#' \item{\code{index}}{What index to use. It includes \code{"PVAF"} or \code{"R2"}. The default is \code{"PVAF"}.}
#' \item{\code{iterative}}{(Iterative) implementation procedure. It includes \code{"none"} (for non-iterative), \code{"test"} (for test-level iterations), \code{"test.att"} (for test-level iterations modifying the least possible amount of q-entries in each iteration), and \code{"item"} (for item-level iterations). The default is \code{"test.att"}.}
#' \item{\code{maxitr}}{Maximum number of iterations if an iterative procedure has been selected. The default is 5.}
#' \item{\code{CDMconv}}{Convergence criteria for the CDM estimations between iterations (only if an iterative procedure has been selected). The default is 0.01.}
#' }
#' @param verbose Show progress? The default is \code{TRUE}.
#'
#' @return \code{modelcompK} returns an object of class \code{modelcompK}.
#' \describe{
#' \item{\code{sug.K}}{The suggested number of attributes for each fit index (\code{vector}). Only if \code{Qs = NULL}.}
#' \item{\code{sel.Q}}{The suggested Q-matrix for each fit index (\code{vector}).}
#' \item{\code{fit}}{The fit indices for each fitted model (\code{matrix}).}
#' \item{\code{exp.exploreK}}{Explored dimensionality (\code{vector}). It can be different from \code{exploreK} if \code{stop} has been used.}
#' \item{\code{usedQ}}{Q-matrices used to fit each model (\code{list}). They will be the estimated (and validated) Q-matrices if \code{Qs = NULL}. Otherwise, they will be \code{Qs}.}
#' \item{\code{specifications}}{Function call specifications (\code{list}).}
#' }
#'
#' @references
#' Ma, W., & de la Torre, J. (2020). GDINA: An R package for cognitive diagnosis modeling. \emph{Journal of Statistical Software}, \emph{93}(14). https://doi.org/10.18637/jss.v093.i14
#'
#' Nájera, P., Abad, F. J., & Sorrel, M. A. (2020). Determining the number of attributes in cognitive diagnosis modeling. \emph{Frontiers in Psychology}, \emph{12}:614470. https://doi.org/10.3389/fpsyg.2021.614470
#'
#' Nájera, P., Sorrel, M. A., de la Torre, J., & Abad, F. J. (2020). Balancing fit and parsimony to improve Q-matrix validation. \emph{British Journal of Mathematical and Statistical Psychology}. https://doi.org/10.1111/bmsp.12228
#'
#' Wang, W., Song, L., & Ding, S. (2018). An exploratory discrete factor loading method for Q-matrix specification in cognitive diagnosis models. In: M. Wilberg, S. Culpepper, R. Janssen, J. González, & D. Molenaar (Eds.), \emph{Quantitative Psychology. IMPS 2017. Springer Proceedings in Mathematics & Statistics} (Vol. 233, pp. 351-362). Springer.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(GDINA)
#' dat <- sim30GDINA$simdat
#' Q <- sim30GDINA$simQ
#'
#' #-------------------------------------
#' # Assess dimensionality from CDM data
#' #-------------------------------------
#' mcK <- modelcompK(dat = dat, exploreK = 1:7, stop = "AIC", val.Q = TRUE, verbose = TRUE)
#' mcK$sug.K # Check suggested number of attributes by each fit index
#' mcK$fit # Check fit indices for each K explored
#' sug.Q <- mcK$usedQ[[mcK$sug.K["AIC"]]] # Extract the suggested Q-matrix by the suggested dimentionality according to AIC
#' sug.Q <- orderQ(sug.Q, Q)$order.Q # Reorder Q-matrix attributes
#' mean(sug.Q == Q) # Check similarity with the generating Q-matrix
#'
#' #--------------------------------------------------
#' # Automatic fit comparison of competing Q-matrices
#' #--------------------------------------------------
#' trueQ <- Q
#' missQ1 <- missQ(Q, .10, seed = 123)$miss.Q
#' missQ2 <- missQ(Q, .20, seed = 456)$miss.Q
#' missQ3 <- missQ(Q, .30, seed = 789)$miss.Q
#' Qs <- list(trueQ, missQ1, missQ2, missQ3)
#' mc <- modelcompK(dat = dat, Qs = Qs, verbose = TRUE)
#' mc$sel.Q # Best-fitting Q-matrix for each fit index
#' mc$fit # Check fit indices for each Q explored
#' }
modelcompK <- function(dat, exploreK = 1:7, Qs = NULL, stop = "none", val.Q = TRUE, estQ.args = list(criterion = "row", cor = "tet", rotation = "oblimin", fm = "uls"), valQ.args = list(index = "PVAF", iterative = "test.att", maxitr = 5, CDMconv = 0.01), verbose = TRUE){
  if(!is.matrix(dat) & !is.data.frame(dat)){stop("Error in modelcompK: dat must be a matrix or data.frame.")}
  if((!is.numeric(exploreK) & !is.double(exploreK))){stop("Error in modelcompK: exploreK must be a numeric vector.")}
  if(any(exploreK < 1)){stop("Error in modelcompK: exploreK must be greater than 0.")}
  if(!is.null(Qs)){if(!is.list(Qs)){stop("Error in modelcompK: Qs must be a list.")}}
  if(!(stop %in% c("none", "AIC", "BIC", "CAIC", "SABIC", "M2", "SRMSR", "RMSEA2", "sig.item.pairs"))){stop("Error in modelcompK: stop must be 'none', 'AIC', 'BIC', 'CAIC', 'SABIC', 'M2', 'SRMSR', 'RMSEA2', 'sig.item.pairs'.")}
  if(!is.logical(val.Q)){stop("Error in modelcompK: val.Q must be logical.")}
  if(is.null(estQ.args$criterion)){estQ.args$criterion <- "row"}
  if(is.null(estQ.args$cor)){estQ.args$cor <- "tet"}
  if(is.null(estQ.args$rotation)){estQ.args$rotation <- "oblimin"}
  if(is.null(estQ.args$fm)){estQ.args$fm <- "uls"}
  if(!is.numeric(estQ.args$criterion) & !is.double(estQ.args$criterion) | length(estQ.args$criterion) > 1){
    if((estQ.args$criterion > 1 | estQ.args$criterion < 0) & !(estQ.args$criterion %in% c("row", "col", "loaddiff"))){stop("Error in modelcompK: estQ.args$criterion must be 'row', 'column', 'loaddiff', or a value between 0 and 1.")}
  }
  if(!(estQ.args$cor %in% c("tet", "cor"))){stop("Error in modelcompK: estQ.args$cor must be 'cor' or 'tet'.")}
  if(!(estQ.args$rotation %in% c("none", "varimax", "quartimax", "bentlerT", "equamax", "varimin", "geominT", "bifactor", "Promax", "promax", "oblimin", "simplimax", "bentlerQ", "geominQ", "biquartimin", "cluster"))){stop("Error in modelcompK: estQ.args$rotation must be 'none', 'varimax', 'quartimax', 'bentlerT', 'equamax', 'varimin', 'geominT', 'bifactor', 'Promax', 'promax', 'oblimin', 'simplimax', 'bentlerQ', 'geominQ', 'biquartimin', or 'cluster'.")}
  if(!(estQ.args$fm %in% c("minres", "uls", "ols", "wls", "gls", "pa", "ml", "minchi", "minrank", "old.min", "alpha"))){stop("Error in modelcompK: estQ.args$fm must be 'minres', 'uls', 'ols', 'wls', 'gls', 'pa', 'ml', 'minchi', 'minrank', 'old.min', or 'alpha'.")}
  if(is.null(valQ.args$index)){valQ.args$index <- "PVAF"}
  if(is.null(valQ.args$iterative)){valQ.args$iterative <- "test.att"}
  if(is.null(valQ.args$maxitr)){valQ.args$maxitr <- 5}
  if(is.null(valQ.args$CDMconv)){valQ.args$CDMconv <- 0.01}
  if(!(valQ.args$index %in% c("PVAF", "R2"))){stop("Error in modelcompK: valQ.args$index must be 'PVAF' or 'R2'.")}
  if(!(valQ.args$iterative %in% c("none", "test", "test.att", "item"))){stop("Error in modelcompK: valQ.args$iterative must be 'none', 'test', 'test.att', or 'item'.")}
  if((!is.numeric(valQ.args$maxitr) & !is.double(valQ.args$maxitr)) | length(valQ.args$maxitr) > 1){stop("Error in modelcompK: valQ.args$maxitr must be a unique numeric value.")}
  if((!is.numeric(valQ.args$CDMconv) & !is.double(valQ.args$CDMconv)) | length(valQ.args$CDMconv) > 1){stop("Error in modelcompK: valQ.args$CDMconv must be a unique numeric value.")}
  if(!is.logical(verbose)){stop("Error in modelcompK: verbose must be logical.")}

  N <- nrow(dat)
  J <- ncol(dat)

  Qs.pre <- Qs
  if(is.null(Qs)){
    fit.res <- matrix(NA, nrow = length(exploreK), ncol = 13)
    colnames(fit.res) <- c("logLik", "np", "AIC", "BIC", "CAIC", "SABIC", "M2", "M2.p", "SRMSR", "RMSEA2", "RMSEA2.low", "RMSEA2.high", "sig.item.pairs")
    rownames(fit.res) <- paste0("K", exploreK)
  } else {
    fit.res <- matrix(NA, nrow = length(Qs), ncol = 13)
    colnames(fit.res) <- c("logLik", "np", "AIC", "BIC", "CAIC", "SABIC", "M2", "M2.p", "SRMSR", "RMSEA2", "RMSEA2.low", "RMSEA2.high", "sig.item.pairs")
    rownames(fit.res) <- paste0("Q", 1:length(Qs))
  }
  fit.res <- as.data.frame(fit.res)
  valQ.conv <- sug.K <- c()

  if(is.null(Qs)){
    if(verbose){
      if(val.Q){cat("\n", "Estimating and validating Q-matrix with K =", exploreK, "\n")}
      if(!val.Q){cat("\n", "Estimating Q-matrix with K =", exploreK, "\n")}
    }
    Qs <- list()
    for(k in exploreK){
      est.Q.info <- estQ(r = dat, K = k, criterion = estQ.args$criterion, efa.args = list(cor = estQ.args$cor, rotation = estQ.args$rotation, fm = estQ.args$fm))
      est.Q <- est.Q.info$est.Q
      fit <- GDINA::GDINA(dat, est.Q, verbose = 0)
      mfit <- GDINA::modelfit(fit)
      if(is.null(mfit$M2)){mfit$M2 <- NA; mfit$M2.pvalue <- NA; mfit$M2.df <- NA; mfit$RMSEA2 <- NA; mfit$RMSEA2.CI <- c(NA, NA)}
      ifit <- GDINA::itemfit(fit)
      if(val.Q){
        if(k > 1){
          sug.Q.info <- valQ(fit, index = valQ.args$index, iterative = valQ.args$iterative, emptyatt = FALSE, maxitr = valQ.args$maxitr, CDMconv = valQ.args$CDMconv, verbose = FALSE)
          sug.Q <- sug.Q.info$sug.Q
          valQ.conv <- c(valQ.conv, sug.Q.info$convergence)
          Qs[[which(exploreK == k)]] <- sug.Q
          if(length(sug.Q.info$sugQ.fit) == 13){
            fit.res[which(exploreK == k),] <- sug.Q.info$sugQ.fit
          } else {
            sug.Q.info$sugQ.fit <- c(sug.Q.info$sugQ.fit[1:6], M2 = NA, M2.p = NA, sug.Q.info$sugQ.fit[7], RMSEA2 = NA, RMSEA2.low = NA, RMSEA2.high = NA, sug.Q.info$sugQ.fit[8])
            fit.res[which(exploreK == k),] <- sug.Q.info$sugQ.fit
          }
        } else {
          Qs[[which(exploreK == k)]] <- est.Q
          fit.res[which(exploreK == k),] <- c(GDINA::extract(fit, what = "logLik"), fit$testfit$npar, AIC(fit), BIC(fit), mfit$CAIC, mfit$SABIC, mfit$M2, mfit$M2.pvalue, mfit$SRMSR, mfit$RMSEA2, mfit$RMSEA2.CI[1], mfit$RMSEA2.CI[2], length(which(ifit$max.itemlevel.fit[,5] < 0.05)))
        }
      } else {
        Qs[[which(exploreK == k)]] <- est.Q
        fit.res[which(exploreK == k),] <- c(GDINA::extract(fit, what = "logLik"), fit$testfit$npar, AIC(fit), BIC(fit), mfit$CAIC, mfit$SABIC, mfit$M2, mfit$M2.pvalue, mfit$SRMSR, mfit$RMSEA2, mfit$RMSEA2.CI[1], mfit$RMSEA2.CI[2], length(which(ifit$max.itemlevel.fit[,5] < 0.05)))
      }
      if(verbose){cat("   k =", k, "explored | AIC =", round(fit.res$AIC[which(exploreK == k)]), "| BIC =", round(fit.res$BIC[which(exploreK == k)]), "\n")}
      if(k > 1 & stop != "none" & which(exploreK == k) > 1){
        if(fit.res[which(exploreK == k), stop] > fit.res[which(exploreK == k) - 1, stop]){break}
      }
    }
  } else {
    if(verbose){cat("\n", "Fitting CDM and extracting fit indices from", length(Qs), "Q-matrices", "\n")}
    for(q in 1:length(Qs)){
      fit <- GDINA::GDINA(dat, Qs[[q]], verbose = 0)
      mfit <- GDINA::modelfit(fit)
      if(is.null(mfit$M2)){mfit$M2 <- NA; mfit$M2.pvalue <- NA; mfit$M2.df <- NA; mfit$RMSEA2 <- NA; mfit$RMSEA2.CI <- c(NA, NA)}
      ifit <- GDINA::itemfit(fit)
      fit.res$logLik[q] <- GDINA::extract(fit, what = "logLik")
      fit.res$np[q] <- fit$testfit$npar
      fit.res$AIC[q] <- AIC(fit)
      fit.res$BIC[q] <- BIC(fit)
      fit.res$CAIC[q] <- mfit$CAIC
      fit.res$SABIC[q] <- mfit$SABIC
      fit.res$M2[q] <- mfit$M2
      fit.res$M2.p[q] <- mfit$M2.pvalue
      fit.res$SRMSR[q] <- mfit$SRMSR
      fit.res$RMSEA2[q] <- mfit$RMSEA2
      fit.res$RMSEA2.low[q] <- mfit$RMSEA2.CI[1]
      fit.res$RMSEA2.high[q] <- mfit$RMSEA2.CI[2]
      fit.res$sig.item.pairs[q] <- length(which(ifit$max.itemlevel.fit[,5] < 0.05))
      if(verbose){cat(paste0("   Q", q), "explored | AIC =", round(fit.res$AIC[q]), "| BIC =", round(fit.res$BIC[q]), "\n")}
    }
  }
  rm.row <- which(apply(fit.res, 1, function(x) all(is.na(x))))
  if(length(rm.row) > 0){fit.res <- fit.res[-rm.row,]}
  names(Qs) <- rownames(fit.res)
  sel.Q.pre <- fit.res[,-c(1:2)]
  if(is.null(Qs.pre)){
    sel.Q <- paste0("K", exploreK[apply(sel.Q.pre, 2, which.min)])
    sel.Q[6] <- paste0("K", exploreK[which.max(sel.Q.pre$M2.p)])
  } else{
    sel.Q <- rownames(sel.Q.pre)[apply(sel.Q.pre, 2, which.min)]
    sel.Q[6] <- rownames(sel.Q.pre)[which.max(sel.Q.pre$M2.p)]
  }
  names(sel.Q) <- colnames(sel.Q.pre)
  if(val.Q){
    if(is.null(Qs.pre)){
      sug.K <- exploreK[apply(sel.Q.pre, 2, which.min)]
      sug.K[6] <- exploreK[which.max(sel.Q.pre$M2.p)]
      names(sug.K) <- names(sel.Q)
    }
  }

  spec <- list(dat = dat, exploreK = exploreK, Qs = Qs.pre, stop = stop, val.Q = val.Q, estQ.args = estQ.args, valQ.args = valQ.args, verbose = verbose)
  if(is.null(Qs.pre)){
    res <- list(sug.K = sug.K, sel.Q = sel.Q, fit = fit.res, usedQ = Qs, specifications = spec)
  } else {
    res <- list(sel.Q = sel.Q, fit = fit.res, usedQ = Qs, specifications = spec)
  }
  class(res) <- "modelcompK"
  return(res)
}
