#' Introduce random misspecifications in Q-matrix
#'
#' @description Introduces random misspecifications in a Q-matrix.
#' Only binary Q-matrix are supported so far.
#' Useful for simulation studies.
#'
#' @param Q A \emph{J} items x \emph{K} attributes Q-matrix (\code{matrix} or \code{data.frame}).
#' @param qjk Number (or proportion, if lower than 1) of q-entries to modify in the Q-matrix.
#' @param retainJ Number of items to retain (i.e., not modify) in the Q-matrix. It will retain the first \code{retainJ} items. It is useful for assuring the completeness of the misspecified Q-matrix if the first items conform one or more identity matrices. The default is 0.
#' @param Qid Assure that the generated Q-matrix is identifiable. It includes \code{"none"} (for no identifiability assurance), \code{"DINA"}, \code{"DINO"}, or \code{"others"} (for other CDMs identifiability). The default is \code{"none"}.
#' @param seed A seed for obtaining consistent results. If \code{NULL}, no seed is used. The default is \code{NULL}.
#'
#' @return \code{missQ} returns an object of class \code{missQ}.
#' \describe{
#' \item{\code{miss.Q}}{The misspecified Q-matrix (\code{matrix}).}
#' \item{\code{Q}}{The input (true) Q-matrix (\code{matrix}).}
#' \item{\code{JK}}{Number of items measuring each attribute (\code{vector}).}
#' \item{\code{Kcor}}{Tetrachoric correlations among the columns (\code{matrix}).}
#' \item{\code{is.Qid}}{Is the generated Q-matrix identifiable under the DINA/DINO models or others CDMs? (\code{vector}).}
#' \item{\code{specifications}}{Function call specifications (\code{list}).}
#' }
#'
#' @references
#' Xu, G., & Shang, Z. (2018). Identifying latent structures in restricted latent class models. \emph{Journal of the American Statistical Association}, \emph{113}, 1284-1295. DOI: 10.1080/01621459.2017.1340889.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Kj <- c(15, 10, 0, 5) # 15 one-attribute, 10 two-attribute, 0 three attribute, and 5 four-attribute items
#' Q <- genQ(J = 30, K = 4, Kj = Kj, Qid = "others", seed = 123)
#' miss.Q <- missQ(Q = Q$gen.Q, qjk = .20, retainJ = 4, seed = 123)
#' }
missQ <- function(Q, qjk, retainJ = 0, Qid = "none", seed = NULL){
  if(!is.matrix(Q) & !is.data.frame(Q)){stop("Error in missQ: Q must be a matrix or data.frame.")}
  J <- nrow(Q)
  K <- ncol(Q)
  if((!is.numeric(qjk) & !is.double(qjk)) | length(qjk) > 1){stop("Error in missQ: qjk must be a unique numeric value.")}
  if(qjk < 1){
    qjk.type <- "prop"
    if(qjk < 0){stop("Error in missQ: qjk must be between 0 and 1 when proportions are provided.")}
  } else {
    qjk.type <- "num"
    if(qjk > (K * J) | qjk < 0){stop("Error in missQ: qjk must be between 0 and total number of q-entries when integers are provided.")}
  }
  if((!is.numeric(retainJ) & !is.double(retainJ)) | length(retainJ) > 1){stop("Error in missQ: retainJ must be a unique numeric value.")}
  if(retainJ > J | retainJ < 0){stop("Error in missQ: retainJ must be between 0 and J.")}
  if(!(Qid %in% c("none", "DINA", "DINO", "others"))){stop("Error in genQ: Qid must be 'none', 'DINA', 'DINO', or 'others'.")}
  if(!is.null(seed)){if((!is.numeric(seed) & !is.double(seed)) | length(seed) > 1){stop("Error in genQ: seed must be a unique numeric value.")}}
  idQ.DINA <- idQ.others <- FALSE

  if(!is.null(seed)){set.seed(seed)}
  if(qjk.type == "prop"){
    qjk.mod <- J * K * qjk
  } else {
    qjk.mod <- qjk
  }
  if(((J - retainJ) * K) < qjk.mod){stop("Error in missQ: lower retainJ or lower qjk required.")}
  if(is.null(retainJ)){
    qjk.miss <- sample(x = 1:(J*K), size = qjk.mod)
  } else {
    qjk.miss <- sample(x = setdiff(1:(J*K), 1:(retainJ*K)), size = qjk.mod)
  }

  cQ <- c(as.matrix(t(Q)))
  cmiss.Q <- cQ
  for(qq in qjk.miss) {
    if(cQ[qq] == 0) {cmiss.Q[qq] <- 1}
    if(cQ[qq] == 1) {cmiss.Q[qq] <- 0}
  }
  miss.Q <- matrix(data = cmiss.Q, nrow = J, ncol = K, byrow = TRUE)
  idQ.DINA <- is.Qid(miss.Q, model = "DINA", verbose = FALSE)$id.Q
  idQ.others <- is.Qid(miss.Q, model = "others", verbose = FALSE)$id.Q

  if(any(rowSums(miss.Q) == 0) & Qid == "none"){
    cmiss.Q <- c(as.matrix(t(miss.Q)))
    q0 <- which(rowSums(miss.Q) == 0)
    for(q in q0){
      qq0_1 <- (q - 1)*K + 1
      qq0_K <- q*K
      qq <- qq0_1:qq0_K
      qjk_mod <- qjk_mods <- qjk.miss[qjk.miss >= qq0_1 & qjk.miss <= qq0_K]
      if(length(qjk_mod) > 1){qjk_mod <- sample(qjk_mod, size = 1)}
      cmiss.Q[qjk_mod] <- 1
      if(length(qq) > length(qjk_mods)){
        qjk_new <- sample(x = setdiff(qq, qjk_mods), size = 1)
        cmiss.Q[qjk_new] <- 1
      }
    }
    miss.Q <- matrix(data = cmiss.Q, nrow = J, ncol = K, byrow = TRUE)
    idQ.DINA <- is.Qid(miss.Q, model = "DINA", verbose = FALSE)$id.Q
    idQ.others <- is.Qid(miss.Q, model = "others", verbose = FALSE)$id.Q
  } else if(any(rowSums(miss.Q) == 0) | (Qid == "DINA" | Qid == "DINO")){
    while(!idQ.DINA | any(rowSums(miss.Q) == 0)){
      cmiss.Q <- c(as.matrix(t(miss.Q)))
      q0 <- which(rowSums(miss.Q) == 0)
      for(q in q0){
        qq0_1 <- (q - 1)*K + 1
        qq0_K <- q*K
        qq <- qq0_1:qq0_K
        qjk_mod <- qjk_mods <- qjk.miss[qjk.miss >= qq0_1 & qjk.miss <= qq0_K]
        if(length(qjk_mod) > 1){qjk_mod <- sample(qjk_mod, size = 1)}
        cmiss.Q[qjk_mod] <- 1
        if(length(qq) > length(qjk_mods)){
          qjk_new <- sample(x = setdiff(qq, qjk_mods), size = 1)
          cmiss.Q[qjk_new] <- 1
        }
      }
      miss.Q <- matrix(data = cmiss.Q, nrow = J, ncol = K, byrow = TRUE)
      idQ.DINA <- is.Qid(miss.Q, model = "DINA", verbose = FALSE)$id.Q
      idQ.others <- is.Qid(miss.Q, model = "others", verbose = FALSE)$id.Q
    }
  } else if(any(rowSums(miss.Q) == 0) | Qid == "others"){
    while(!idQ.others | any(rowSums(miss.Q) == 0)){
      cmiss.Q <- c(as.matrix(t(miss.Q)))
      q0 <- which(rowSums(miss.Q) == 0)
      for(q in q0){
        qq0_1 <- (q - 1)*K + 1
        qq0_K <- q*K
        qq <- qq0_1:qq0_K
        qjk_mod <- qjk_mods <- qjk.miss[qjk.miss >= qq0_1 & qjk.miss <= qq0_K]
        if(length(qjk_mod) > 1){qjk_mod <- sample(qjk_mod, size = 1)}
        cmiss.Q[qjk_mod] <- 1
        if(length(qq) > length(qjk_mods)){
          qjk_new <- sample(x = setdiff(qq, qjk_mods), size = 1)
          cmiss.Q[qjk_new] <- 1
        }
      }
      miss.Q <- matrix(data = cmiss.Q, nrow = J, ncol = K, byrow = TRUE)
      idQ.DINA <- is.Qid(miss.Q, model = "DINA", verbose = FALSE)$id.Q
      idQ.others <- is.Qid(miss.Q, model = "others", verbose = FALSE)$id.Q
    }
  }
  JK <- colSums(miss.Q)
  Kcor <- round(sirt::tetrachoric2(miss.Q)$rho, 3)

  spec <- list(Q = Q, qjk = qjk, retainJ = retainJ, Qid = Qid, seed = seed)
  res <- list(miss.Q = miss.Q, Q = Q, JK = JK, Kcor = Kcor, is.Qid = c("DINA/O" = idQ.DINA, "others" = idQ.others), specifications = spec)
  class(res) <- "missQ"
  return(res)
}
