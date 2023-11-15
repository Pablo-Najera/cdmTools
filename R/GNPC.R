#' General nonparametric classification method
#'
#' @description Attribute profile estimation using the \emph{general nonparametric classification method} (GNPC; Chiu, Sun, & Bian, 2018).
#' The GNPC can be considered as a robust alternative to the parametric G-DINA model with low sample sizes.
#' The \code{AlphaNP} function from the \code{NPCD} package (Zheng & Chiu, 2019; Chiu, Sun, & Bian, 2018) using weighted Hamming distances is used to initiate the procedure.
#'
#' @param dat A \emph{N} individuals x \emph{J} items (\code{matrix} or \code{data.frame}). Missing values need to be coded as \code{NA}. Caution is advised if missing data are present.
#' @param Q A \emph{J} items x \emph{K} attributes Q-matrix (\code{matrix} or \code{data.frame}).
#' @param initiate Should the conjunctive (\code{"AND"}) or disjunctive (\code{"OR"}) NPC be used to initiate the procedure? Default is \code{"AND"}.
#' @param min.change Minimum proportion of modified attribute profiles to use as a stopping criterion. Default is .001.
#' @param maxitr Maximum number of iterations. Default is 1000.
#' @param verbose Print information after each iteration. Default is \code{TRUE}.
#'
#' @return \code{GNPC} returns an object of class \code{GNPC}.
#' \describe{
#' \item{\code{alpha.est}}{Estimated attribute profiles (\code{matrix}).}
#' \item{\code{loss.matrix}}{The distances between the weighted ideal responses from each latent class (rows) and examinees' observed responses (columns) (\code{matrix}).}
#' \item{\code{eta.w}}{The weighted ideal responses for each latent class (rows) on each item (columns) (\code{matrix}).}
#' \item{\code{w}}{The estimated weights, used to compute the weighted ideal responses (\code{matrix}).}
#' \item{\code{n.ite}}{Number of iterations required to achieve convergence (\code{double}).}
#' \item{\code{hist.change}}{Proportion of modified attribute profiles in each iteration (\code{vector}).}
#' \item{\code{specifications}}{Function call specifications (\code{list}).}
#' }
#'
#' @author {Pablo Nájera, Universidad Pontificia Comillas}
#'
#' @references
#' Chiu, C.-Y., & Douglas, J. (2013). A nonparametric approach to cognitive diagnosis by proximity to ideal response patterns. \emph{Journal of Classification}, \emph{30}, 225-250. DOI: 10.1007/s00357-013-9132-9
#'
#' Chiu, C.-Y., Sun, Y., & Bian, Y. (2018). Cognitive diagnosis for small education programs: The general nonparametric classification method. \emph{Psychometrika}, \emph{83}, 355-375. DOI: 10.1007/s11336-017-9595-4
#'
#' Zheng, Y., & Chiu, C.-Y. (2019). \emph{NPCD: Nonparametric methods for cognitive diagnosis}. R package version 1.0-11. https://cran.r-project.org/web/packages/NPCD/.
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(GDINA)
#' Q <- sim30GDINA$simQ # Q-matrix
#' K <- ncol(Q)
#' J <- nrow(Q)
#' set.seed(123)
#' GS <- data.frame(guessing = rep(0.1, J), slip = rep(0.1, J))
#' sim <- simGDINA(200, Q, GS)
#' simdat <- sim$dat # Simulated data
#' simatt <- sim$attribute # Generating attributes
#' fit.GNPC <- GNPC(simdat, Q) # Apply the GNPC method
#' ClassRate(fit.GNPC$alpha.est, simatt) # Check classification accuracy
#' }
GNPC <- function(dat, Q, initiate = "AND", min.change = 0.001, maxitr = 1000, verbose = TRUE){
  K <- ncol(Q)
  J <- nrow(Q)
  N <- nrow(dat)
  Pj <- colMeans(dat)
  lc <- GDINA::attributepattern(K); rownames(lc) <- 1:nrow(lc)
  lc.ch <- apply(lc, 1, function(x) paste(x, collapse = ""))
  L <- length(lc.ch)

  if(N > 1){
    att.clas.wHD <- cdmTools.AlphaNP(dat, Q, initiate, method = "Weighted")$alpha.est
  } else {
    att.clas.wHD <- cdmTools.AlphaNP(dat, Q, initiate)$alpha.est
  }
  change.h <- c()
  change <- 1
  ite <- 0

  while(change >= min.change){
    ite <- ite + 1
    prev.att.clas.wHD <- att.clas.wHD
    w <- eta.w <- matrix(NA, nrow = L, ncol = J, dimnames = list(lc.ch, paste("J", 1:J, sep = "")))
    for(j in 1:J){
      kj <- which(Q[j,] == 1)
      Kj <- length(kj)
      groups <- GDINA::attributepattern(Kj)
      for(lg in 1:nrow(groups)){
        eta_k.c <- eta_k.d <- c()
        for(k in 1:Kj){
          eta_k.c <- c(eta_k.c, groups[lg, k]^Q[j, kj[k]])
          eta_k.d <- c(eta_k.d, (1 - groups[lg, k])^Q[j, kj[k]])
        }
        eta.c <- prod(eta_k.c)
        eta.d <- 1 - prod(eta_k.d)
        Ninlg <- sapply(1:N, function(i) all(att.clas.wHD[i, kj, drop = F] == groups[lg,]))
        dat.lg <- dat[Ninlg,,drop = FALSE]
        Cl <- nrow(dat.lg)
        if(is.na(sum(dat.lg[,j] - eta.d) / (Cl * (eta.c - eta.d)))){
          wlj <- switch(initiate, "AND" = 1, "OR" = 0)
        } else if(sum(dat.lg[,j] - eta.d) / (Cl * (eta.c - eta.d)) == Inf){
          wlj <- 1
        } else if(sum(dat.lg[,j] - eta.d) / (Cl * (eta.c - eta.d)) == -Inf){
          wlj <- 0
        } else {
          wlj <- sum(dat.lg[,j] - eta.d) / (Cl * (eta.c - eta.d))
        }
        lg2lc <- sapply(1:L, function(l) all(lc[l, kj, drop = F] == groups[lg,]))
        eta.w[lg2lc, j] <- wlj * eta.c + (1 - wlj) * eta.d
        if(!is.na(sum(dat.lg[,j] - eta.d) / (Cl * (eta.c - eta.d)))){
          w[lg2lc, j] <- wlj
        } else {
          w[lg2lc, j] <- NA
        }
      }
    }

    att.clas.wHD <- matrix(NA, N, K, dimnames = list(1:N, paste("A", 1:K, sep = "")))
    for(i in 1:N){
      wHD <- sweep(eta.w[,!is.na(dat[i,])], 2, dat[i,][!is.na(dat[i,])])^2
      wHD <- rowSums(wHD)
      att.clas.wHD[i,] <- as.double(substring(names(wHD)[which.min(wHD)], seq(1, K, 1), seq(1, K, 1)))
    }

    i.change <- as.vector(which(rowSums(abs(att.clas.wHD - prev.att.clas.wHD)) > 0))
    change <- length(which(rowSums(abs(att.clas.wHD - prev.att.clas.wHD)) > 0)) / N
    change.h <- c(change.h, change)

    if(ite > 1){
      if(identical(i.change, prev.i.change)){
        if(identical(att.clas.wHD[i.change,], two.prev.att.clas.wHD[prev.i.change,])){break}
      }
    }

    two.prev.att.clas.wHD <- prev.att.clas.wHD
    prev.i.change <- i.change
    if(verbose){cat(paste0("\r", "Iter = ", ite, " | Change rate = ", round(change, 3), "          "))}

    if(length(change.h) > 5){if(length(unique(change.h[(length(change.h) - 4):length(change.h)])) < 2){break}}
    if(length(change.h) > 10){if(length(unique(change.h[(length(change.h) - 9):length(change.h)])) < 3){break}}
    if(length(change.h) > 15){if(length(unique(change.h[(length(change.h) - 14):length(change.h)])) < 4){break}}
    if(length(change.h) > 20){if(length(unique(change.h[(length(change.h) - 19):length(change.h)])) < 5){break}}
    if(length(change.h) > 25){if(length(unique(change.h[(length(change.h) - 24):length(change.h)])) < 6){break}}
    if(ite == maxitr){break}
  }

  w <- t(w)
  eta.w <- t(eta.w)
  d <- apply(dat, 1, function(i) apply(i - eta.w, 2, function(l) sum(l^2)))
  specs <- list(dat = dat, Q = Q, initiate = initiate, min.change = min.change, maxitr = maxitr, verbose = verbose)
  res <- list(alpha.est = att.clas.wHD, loss.matrix = d, eta.w = eta.w, w = w, n.ite = ite, hist.change = change.h, specifications = specs)
  class(res) <- "GNPC"
  return(res)
}
