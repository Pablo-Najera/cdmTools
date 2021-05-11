#' Reorder Q-matrix columns
#'
#' @description Reorders Q-matrix columns according to a target matrix (e.g., another Q-matrix).
#' Specifically, it provides a reordered Q-matrix which columns show the lowest possible average Tucker index congruent coefficient with the target columns.
#' Reordering a Q-matrix is alike relabeling the attributes and it does not change the model.
#' Useful for simulation studies (e.g., comparing a validated Q-matrix with the generating Q-matrix).
#'
#' @param Q A \emph{J} items × \emph{K} attributes Q-matrix (\code{matrix} or \code{data.frame}). This is the Q-matrix that will be reordered.
#' @param target A \emph{J} items × \emph{K} attributes Q-matrix (\code{matrix} or \code{data.frame}). This could be the "true", generating Q-matrix.
#'
#' @return \code{orderQ} returns an object of class \code{orderQ}.
#' \describe{
#' \item{\code{order.Q}}{The reordered Q-matrix (\code{matrix}).}
#' \item{\code{configs}}{Comparison information between the different column configurations of the Q-matrix and the target Q-matrix, including the average absolute difference and the average Tucker index of factor congruence (\code{matrix}). The function will not look for all possible specifications if a perfect match is found.}
#' \item{\code{specifications}}{Function call specifications (\code{list}).}
#' }
#'
#' @export
#'
#' @examples
#' library(GDINA)
#' dat <- sim30GDINA$simdat
#' Q <- sim30GDINA$simQ
#' sugQ1 <- estQ(r = dat, K = 5) # Estimate Q-matrix
#' sugQ1$est.Q <- orderQ(sugQ1$est.Q, Q)$order.Q # Reorder Q-matrix attributes
#' mean(sugQ1$est.Q == Q) # Check similarity with the generating Q-matrix
orderQ <- function(Q, target){
  if(!is.matrix(Q) & !is.data.frame(Q)){stop("Error in orderQ: Q must be a matrix or data.frame.")}
  if(!is.matrix(target) & !is.data.frame(target)){stop("Error in orderQ: target must be a matrix or data.frame.")}

  mQ <- Q
  mT <- target

  if(!is.matrix(mQ)){mQ <- as.matrix(mQ)}
  if(!is.matrix(mT)){mT <- as.matrix(mT)}

  K <- ncol(mQ)
  order.k <- combinat::permn(1:K)
  order.k.col <- as.numeric(sapply(order.k, paste, collapse = ""))
  MAD <- CC <- c()
  for(i in 1:length(order.k)){
    tmp.MAD <- round(mean(apply(abs(mT - mQ[,order.k[[i]]]), 2, mean)), 3)
    tmp.CC <- round(mean(diag((t(mT) %*% mQ[,order.k[[i]]]) / (sqrt(colSums(mT^2) %o% colSums(mQ[,order.k[[i]]]^2))))), 3)
    if(is.na(tmp.MAD)){tmp.MAD <- 999}
    if(is.na(tmp.CC)){tmp.CC <- 999}
    MAD <- c(MAD, tmp.MAD)
    CC <- c(CC, tmp.CC)
    if(MAD[i] == 0 | CC[i] == 1){break}
  }
  order.Q <- mQ[,order.k[[which.max(CC)]]]
  rownames(order.Q) <- 1:nrow(order.Q)
  configs <- cbind(order = order.k.col[1:length(MAD)], MAD, CC)
  configs <- configs[order(configs[,"order"]),,drop = FALSE]
  rownames(configs) <- 1:nrow(configs)

  spec <- list(Q = Q, target = target)
  res <- list(order.Q = order.Q, configs = configs, specifications = spec)
  class(res) <- "orderQ"
  return(res)
}
