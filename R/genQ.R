#' Generate Q-matrix
#'
#' @description Generates a Q-matrix.
#' The criteria from Chen, Liu, Xu, & Ying (2015) and Xu & Shang (2018) can be used to generate identifiable Q-matrices.
#' Only binary Q-matrix are supported so far.
#' Useful for simulation studies.
#'
#' @param J Number of items.
#' @param K Number of attributes.
#' @param Kj A vector specifying the number (or proportion, if summing up to 1) of items measuring 1, 2, 3, ..., attributes. The first element of the vector determines the number (or proportion) of items measuring 1 attribute, and so on. See \code{Examples}.
#' @param I Number of identity matrices to include in the Q-matrix (up to column permutation). The default is 2.
#' @param min.JK Minimum number of items measuring each attribute. It can be overwritten by \code{I}, if \code{I} is higher than \code{min.JK}. The default is 3.
#' @param max.Kcor Maximum allowed tetrachoric correlation among the columns to avoid overlapping (Nájera, Sorrel, de la Torre, & Abad, 2020). The default is 1.
#' @param Qid Assure that the generated Q-matrix is identifiable. It includes \code{"none"} (for no identifiability assurance), \code{"DINA"}, \code{"DINO"}, or \code{"others"} (for other CDMs identifiability). The default is \code{"none"}.
#' @param seed A seed for obtaining consistent results. If \code{NULL}, no seed is used. The default is \code{NULL}.
#'
#' @return \code{genQ} returns an object of class \code{genQ}.
#' \describe{
#' \item{\code{gen.Q}}{The generated Q-matrix (\code{matrix}).}
#' \item{\code{JK}}{Number of items measuring each attribute (\code{vector}).}
#' \item{\code{Kcor}}{Tetrachoric correlations among the columns (\code{matrix}).}
#' \item{\code{is.Qid}}{Is the generated Q-matrix identifiable under the DINA/DINO models or others CDMs? (\code{vector}).}
#' \item{\code{specifications}}{Function call specifications (\code{list}).}
#' }
#'
#' @references
#' Chen, Y., Liu, J., Xu, G., & Ying, Z. (2015). Statistical analysis of Q-matrix based diagnostic classification models. \emph{Journal of the American Statistical Association}, \emph{110}, 850-866. https://doi.org/10.1080/01621459.2014.934827
#'
#' Nájera, P., Sorrel, M. A., de la Torre, J., & Abad, F. J. (2020). Balancing fit and parsimony to improve Q-matrix validation. \emph{British Journal of Mathematical and Statistical Psychology}. https://doi.org/10.1111/bmsp.12228
#'
#' Xu, G., & Shang, Z. (2018). Identifying latent structures in restricted latent class models. \emph{Journal of the American Statistical Association}, \emph{113}, 1284-1295. https://doi.org/10.1080/01621459.2017.1340889
#'
#' @export
#'
#' @examples
#' Kj <- c(15, 10, 0, 5) # 15 one-att, 10 2-atts, 0 3-atts, and 5 four-atts items
#' Q <- genQ(J = 30, K = 4, Kj = Kj, Qid = "others", seed = 123)
genQ <- function(J, K, Kj, I = 2, min.JK = 3, max.Kcor = 1, Qid = "none", seed = NULL){
  if((!is.numeric(J) & !is.double(J)) | length(J) > 1){stop("Error in genQ: J must be a unique numeric value.")}
  if((!is.numeric(K) & !is.double(K)) | length(K) > 1){stop("Error in genQ: K must be a unique numeric value.")}
  entries <- J * K
  if(length(Kj) > K){stop("Error in genQ: the length of Kj cannot be greater than K.")}
  if(any(Kj < 1 & Kj > 0)){
    Kj.type <- "prop"
    if(sum(Kj) != 1){stop("Error in genQ: Kj must sum up to 1 when proportions are provided.")}
  } else {
    Kj.type <- "num"
    if(sum(Kj) != J){stop("Error in genQ: Kj must sum up to J when integers are provided.")}
  }
  if(I * K > J){stop("Error in genQ: a higher J is required to be able to provide the desired number of identity matrices.")}
  if(I * K + (min.JK - I) * K > J){stop("Error in genQ: a higher J is required to be able to provide the desired number of identity matrices and number of items measuring each attribute.")}
  if((!is.numeric(max.Kcor) & !is.double(max.Kcor)) | length(max.Kcor) > 1){stop("Error in genQ: max.Kcor must be a unique numeric value.")}
  if(max.Kcor > 1 | max.Kcor < 0){stop("Error in genQ: max.Kcor must be a value between 0 and 1.")}
  if(!(Qid %in% c("none", "DINA", "DINO", "others"))){stop("Error in genQ: Qid must be 'none', 'DINA', 'DINO', or 'others'.")}
  if(!is.null(seed)){if((!is.numeric(seed) & !is.double(seed)) | length(seed) > 1){stop("Error in genQ: seed must be a unique numeric value.")}}
  idQ.DINA <- idQ.others <- FALSE

  if(!is.null(seed)){set.seed(seed)}
  Q <- matrix(rep(diag(1, K), I), ncol = K, byrow = T)
  pat <- GDINA::attributepattern(K)
  for(k in 1:length(Kj)){assign(paste0("pat", k), rbind(pat[which(rowSums(pat) == k),]))}
  Jk <- c()
  if(Kj.type == "prop"){
    for(k in 1:length(Kj)){Jk <- c(Jk, floor(J * Kj[k]))}
    while(sum(Jk) < J){
      if(nrow(Q) - Jk[1] == 1){
        Jk[1] <- Jk[1] + 1
      } else if(nrow(Q) - Jk[1] > 1){
        stop("Error in genQ: More items or less 1-attribute q-vectors required to generate the Q-matrix.")
      } else {
        tmp <- sample(1:length(Jk), 1)
        Jk[tmp] <- Jk[tmp] + 1
      }
    }
  } else {
    Jk <- Kj
  }

  for(k in 1:length(Kj)){assign(paste0("replace", k), ifelse(Jk[k] > nrow(get(paste0("pat", k))), T, F))}
  if(Qid %in% c("DINA", "DINO") & Jk[1] < (2 * K)){warning("Warning in genQ: The Q-matrix cannot be identified for the provided CDMs unless two identity matrices are included.")}
  if(Qid %in% c("others") & Jk[1] < (3 * K)){warning("Warning in genQ: The Q-matrix cannot be identified for the provided CDMs unless three identity matrices are included.")}
  Jk[1] <- Jk[1] - nrow(Q)
  for(k in 1:length(Kj)){Q <- rbind(Q, rbind(get(paste0("pat", k))[sample(nrow(get(paste0("pat", k))), size = Jk[k], replace = get(paste0("replace", k))),]))}
  corQ <- sirt::tetrachoric2(Q)$rho
  if(Qid == "none"){
    while(any(colSums(Q[(K*I + 1):J,]) < min.JK) | any(corQ[lower.tri(corQ)] > max.Kcor)){
      Q <- matrix(rep(diag(1, K), I), ncol = K, byrow = T)
      for(k in 1:length(Kj)){Q <- rbind(Q, rbind(get(paste0("pat", k))[sample(nrow(get(paste0("pat", k))), size = Jk[k], replace = get(paste0("replace", k))),]))}
      idQ.DINA <- is.Qid(Q, model = "DINA", verbose = FALSE)$id.Q
      idQ.others <- is.Qid(Q, model = "others", verbose = FALSE)$id.Q
      corQ <- sirt::tetrachoric2(Q)$rho
    }
  } else if(Qid == "DINA" | Qid == "DINO"){
    while(any(colSums(Q[(K*I + 1):J,]) < min.JK) | any(corQ[lower.tri(corQ)] > max.Kcor) | idQ.DINA == FALSE){
      Q <- matrix(rep(diag(1, K), I), ncol = K, byrow = T)
      for(k in 1:length(Kj)){Q <- rbind(Q, rbind(get(paste0("pat", k))[sample(nrow(get(paste0("pat", k))), size = Jk[k], replace = get(paste0("replace", k))),]))}
      idQ.DINA <- is.Qid(Q, model = "DINA", verbose = FALSE)$id.Q
      idQ.others <- is.Qid(Q, model = "others", verbose = FALSE)$id.Q
      corQ <- sirt::tetrachoric2(Q)$rho
    }
  } else if(Qid == "others"){
    while(any(colSums(Q[(K*I + 1):J,]) < min.JK) | any(corQ[lower.tri(corQ)] > max.Kcor) | idQ.others == FALSE){
      Q <- matrix(rep(diag(1, K), I), ncol = K, byrow = T)
      for(k in 1:length(Kj)){Q <- rbind(Q, rbind(get(paste0("pat", k))[sample(nrow(get(paste0("pat", k))), size = Jk[k], replace = get(paste0("replace", k))),]))}
      idQ.DINA <- is.Qid(Q, model = "DINA", verbose = FALSE)$id.Q
      idQ.others <- is.Qid(Q, model = "others", verbose = FALSE)$id.Q
      corQ <- sirt::tetrachoric2(Q)$rho
    }
  }
  idQ.DINA <- is.Qid(Q, model = "DINA", verbose = FALSE)$id.Q
  idQ.others <- is.Qid(Q, model = "others", verbose = FALSE)$id.Q
  corQ <- sirt::tetrachoric2(Q)$rho
  Jk[1] <- Jk[1] + (I * K)

  spec <- list(J = J, K = K, Kj = Kj, I = I, min.JK = min.JK, max.Kcor = max.Kcor, Qid = Qid, seed = seed)
  res <- list(gen.Q = Q, JK = Jk, Kcor = round(corQ, 3), is.Qid = c("DINA/O" = idQ.DINA, "others" = idQ.others), specifications = spec)
  class(res) <- "genQ"
  return(res)
}
