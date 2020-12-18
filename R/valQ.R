#' Empirical Q-matrix validation
#'
#' @description Empirical Q-matrix validation using the \emph{Hull} method (Najera, Sorrel, de la Torre, & Abad, 2020a).
#' The procedure can be used either with the PVAF (de la Torre & Chiu, 2016) or McFadden's pseudo R-squared (McFadden, 1974).
#' The PVAF is recommended (Najera, Sorrel, de la Torre, & Abad, 2020a).
#' Note that the pseudo R-squared might not be computationally feasible for highly dimensional Q-matrices, say more than 10 attributes.
#' Different iterative implementations are available, such as the test-level implementation (see Terzi & de la Torre, 2018), attribute-test-level implementation (Najera, Sorrel, de la Torre, & Abad, 2020), and item-level implementation (Najera, Sorrel, de la Torre, & Abad, 2020b).
#' If an iterative implementation is used, the \code{GDINA} R package (Ma & de la Torre, 2020) is used for the calibration of the CDMs.
#'
#' @param fit A G-DINA model fit object from the \code{GDINA} package (Ma & de la Torre, 2020).
#' @param index What index to use. It includes \code{"PVAF"} or \code{"R2"}. The default is \code{"PVAF"}.
#' @param iterative (Iterative) implementation procedure. It includes \code{"none"} (for non-iterative), \code{"test"} (for test-level iterations), \code{"test.att"} (for attribute-test-level), and \code{"item"} (for item-level iterations). The default is \code{"test.att"}.
#' @param emptyatt Is it possible for the suggested Q-matrix to have an empty attribute (i.e., an attribute not measured by any item)? Although rarely, it is possible for iterative procedures to provide a suggested Q-matrix in which one or more attributes are empty. This might indicate that the original Q-matrix had more attributes than necessary. If \code{FALSE}, then at least one item (i.e., the one that is most likely) will measure each attribute in the suggested Q-matrix. The default is \code{TRUE}.
#' @param maxitr Maximum number of iterations if an iterative procedure has been selected. The default is 100.
#' @param CDMconv Convergence criteria for the CDM estimations between iterations (only if an iterative procedure has been selected). The default is 0.0001.
#' @param verbose Print information after each iteration if an iterative procedure is used. The default is \code{TRUE}.
#'
#' @return \code{valQ} returns an object of class \code{valQ}.
#' \describe{
#' \item{\code{sug.Q}}{Suggested Q-matrix (\code{matrix}).}
#' \item{\code{Q}}{Original Q-matrix (\code{matrix}).}
#' \item{\code{sugQ.fit}}{Several fit indices from the model obtained with the suggested Q-matrix (\code{vector}).}
#' \item{\code{index}}{PVAF or pseudo R-squared (depending on which one was used) for each item (\code{matrix}).}
#' \item{\code{iter.Q}}{Q-matrices used in each iteration (\code{list}). Provided only if an iterative procedure has been used.}
#' \item{\code{iter.index}}{PVAF or pseudo R-squared (depending on which one was used) for each item in each iteration (\code{list}). Provided only if an iterative procedure has been used.}
#' \item{\code{n.iter}}{Number of iterations used (\code{double}). Provided only if an iterative procedure has been used.}
#' \item{\code{convergence}}{Convergence information (\code{double}). It can be 1 (convergence), 2 (lack of convergence: maximum number of iterations achieved), 3 (lack of convergence: empty attribute obtained), and 4 (lack of convergence: loop Q-matrices). Provided only if an iterative procedure has been used.}
#' \item{\code{time}}{Initial and finish time (\code{vector}).}
#' \item{\code{time.used}}{Total computation time (\code{difftime}).}
#' \item{\code{specifications}}{Function call specifications (\code{list}).}
#' }
#'
#' @references
#' de la Torre, J., & Chiu, C.-Y. (2016). A general method of empirical Q-matrix validation. \emph{Psychometrika}, \emph{81}, 253-273. DOI: 10.1007/s11336-015-9467-8.
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R package for cognitive diagnosis modeling. \emph{Journal of Statistical Software}, \emph{93}(14). DOI: 10.18637/jss.v093.i14.
#'
#' McFadden, D. (1974). Conditional logit analysis of qualitative choice behavior. In P. Zarembka (Ed.), \emph{Frontiers in Economics} (pp. 105-142). Academic Press.
#'
#' Najera, P., Sorrel, M. A., de la Torre, J., & Abad, F. J. (2020a). Balancing fit and parsimony to improve Q-matrix validation. \emph{British Journal of Mathematical and Statistical Psychology}.
#'
#' Najera, P., Sorrel, M. A., de la Torre, J., & Abad, F. J. (2020b). Improving robustness in Q-matrix validation using an iterative and dynamic procedure. \emph{Applied Psychological Measurement}. DOI: 10.1177/0146621620909904.
#'
#' Terzi, R., & de la Torre, J. (2018). An iterative method for empirically-based Q-matrix validation. \emph{International Journal of Assessment Tools in Education}, \emph{5}, 248-262. DOI: 10.21449/ijate.407193.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(GDINA)
#' dat <- sim30GDINA$simdat
#' Q <- sim30GDINA$simQ # Generating Q-matrix
#' miss.Q <- missQ(Q = Q, qjk = .30, retainJ = 5, seed = 123)$miss.Q # Misspecified Q-matrix
#' fit <- GDINA(dat, miss.Q) # GDINA object
#' sug.Q <- valQ(fit = fit, verbose = TRUE) # Hull method for Q-matrix validation
#' mean(sug.Q$sug.Q == Q) # Check similarity with the generating Q-matrix
#' }
valQ <- function(fit, index = "PVAF", iterative = "test.att", emptyatt = TRUE, maxitr = 100, CDMconv = 0.0001, verbose = TRUE){
  if(!(class(fit) == "GDINA")){stop("Error in valQ: fit must be of class 'GDINA'.")}
  if(!(index %in% c("PVAF", "R2"))){stop("Error in valQ: index must be 'PVAF' or 'R2'.")}
  if(!(iterative %in% c("none", "test", "test.att", "item"))){stop("Error in valQ: iterative must be 'none', 'test', 'test.att', or 'item'.")}
  if(!is.logical(emptyatt)){stop("Error in valQ: emptyatt must be logical.")}
  if((!is.numeric(maxitr) & !is.double(maxitr)) | length(maxitr) > 1){stop("Error in valQ: maxitr must be a unique numeric value.")}
  if((!is.numeric(CDMconv) & !is.double(CDMconv)) | length(CDMconv) > 1){stop("Error in valQ: CDMconv must be a unique numeric value.")}
  if(!is.logical(verbose)){stop("Error in valQ: verbose must be logical.")}

  dat <- fit$options$dat
  Q <- as.matrix(fit$options$Q)
  J <- nrow(Q)
  K <- ncol(Q)
  q.pattern <- GDINA::attributepattern(K)[-1,]
  iter.Q <- iter.index <- list()
  iter.Q[["0"]] <- test.Q <- prov.Q <- Q

  i <- 1
  init.time <- Sys.time()

  if(iterative == "none"){
    M <- extract.Hull(fit, index, R2.aux = NULL)[[1]]
    tmp <- select.q(M, NULL)
    prov.Q <- tmp$sug.Q
    convergence <- 1 # YES convergence
    fit <- GDINA::GDINA(dat, prov.Q, verbose = 0, control = list(conv.crit = CDMconv))
  } else{ # iterative
    while(i <= maxitr){
      M <- extract.Hull(fit, index, R2.aux = NULL)[[1]]
      tmp <- select.q(M, best.pos = NULL)
      iter.index[[as.character(i - 1)]] <- tmp$bKj$best.M
      hit <- hit.cand <- which(sapply(1:J, function(j) any(prov.Q[j,] != tmp$sug.Q[j,])))
      if(length(hit.cand) == 0){i <- i - 1; convergence <- 1; break} # YES convergence
      if(length(hit.cand) > 1){
        Kj.prov <- rowSums(prov.Q[hit.cand,, drop = F])
        Kj.sug <- rowSums(tmp$sug.Q[hit.cand,, drop = F])
        if(iterative == "test.att"){
          att.change <- ifelse(Kj.prov < Kj.sug, 1, ifelse(Kj.prov == Kj.sug, 0, -1))
          new.att <- Kj.prov + att.change
          for(r in 1:length(hit)){
            att.r <- new.att[r]
            while(is.na(tmp$H[att.r, hit[r]])){att.r <- att.r + att.change[r]}
            q.att <- q.pattern[tmp$bKj$best.pos[att.r, hit[r]],]
            tmp$sug.Q[hit[r],] <- q.att
          }
        }
        if(iterative == "item"){
          prov.q <- match(apply(prov.Q[hit.cand,, drop = F], 1, paste, collapse = ""), apply(q.pattern, 1, paste, collapse = ""))
          prov.q.in.best.q <- sapply(1:length(hit.cand), function(j) prov.q[j] %in% tmp$bKj$best.pos[,hit.cand, drop = F][,j])
          M.diff <- c()
          for(j in 1:length(hit.cand)){
            if(prov.q.in.best.q[j]){
              M.diff <- c(M.diff, tmp$H[Kj.sug[j], hit.cand[j]] - tmp$H[Kj.prov[j], hit.cand[j]])
            } else {
              tmp.best <- tmp$bKj$best.pos
              tmp.best[Kj.prov[j], hit.cand[j]] <- prov.q[j]
              tmpj <- select.q(M, tmp.best)
              M.diff <- c(M.diff, tmp$H[Kj.sug[j], hit.cand[j]] - tmpj$H[Kj.prov[j], hit.cand[j]])
            }
          }
          hit <- hit.cand[which.max(M.diff)]
        }
      }
      test.Q[hit,] <- tmp$sug.Q[hit,]
      if(any(sapply(1:i, function(x) mean(test.Q == iter.Q[[x]])) == 1)){convergence <- 2; break} # NO convergence (loop)
      if(any(colSums(test.Q) == 0)){convergence <- 3; break} # NO convergence (null attribute)
      prov.Q <- iter.Q[[as.character(i)]] <- test.Q
      if(verbose){cat("Iteration =", sprintf("%03d", i), "| Item/s modified =", paste(hit, collapse = ", "), "\n")}
      fit <- GDINA::GDINA(dat, prov.Q, verbose = 0, control = list(conv.crit = CDMconv))
      i <- i + 1
    }
  }

  n.iter <- i
  if(n.iter > maxitr){convergence <- 4} # NO convergence (max.iter reached)
  if(emptyatt & convergence == 3){
    sug.Q <- test.Q
  } else {
    sug.Q <- prov.Q
  }

  mfit <- GDINA::modelfit(fit)
  ifit <- GDINA::itemfit(fit)
  sugQ.fit <- c(logLik = GDINA::extract(fit, what = "logLik"), np = fit$testfit$npar, AIC = AIC(fit), BIC = BIC(fit), CAIC = mfit$CAIC, SABIC = mfit$SABIC, M2 = mfit$M2, M2.pvalue = mfit$M2.pvalue, SRMSR = mfit$SRMSR, RMSEA2 = mfit$RMSEA2, RMSEA2.CI = mfit$RMSEA2.CI[1], RMSEA2.CI = mfit$RMSEA2.CI[2], sig.item.pairs = length(which(ifit$max.itemlevel.fit[,5] < 0.05)))
  finish.time <- Sys.time()
  time.used <- difftime(finish.time, init.time, units = "secs")
  spec <- list(fit = fit, index = index, iterative = iterative, emptyatt = emptyatt, maxitr = maxitr, CDMconv = CDMconv, verbose = verbose)
  res <- list(sug.Q = sug.Q, Q = Q, sugQ.fit = sugQ.fit, index = tmp$bKj$best.M, iter.Q = iter.Q, iter.index = iter.index, n.iter = n.iter, convergence = convergence, time = c(initial = init.time, finish = finish.time), time.used = time.used, specifications = spec)
  class(res) <- "valQ"
  return(res)
}
