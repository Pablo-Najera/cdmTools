#' @title Check whether a Q-matrix is identifiable
#'
#' @description Checks whether a Q-matrix fulfills the conditions for strict and generic identifiability according to Gu & Xu (2021).
#'
#' @param Q A \emph{J} items x \emph{K} attributes Q-matrix (\code{matrix} or \code{data.frame}).
#' @param model CDM to be considered. It includes \code{"DINA"}, \code{"DINO"}, or \code{"others"} (for other CDMs; e.g., G-DINA, A-CDM).
#'
#' @return \code{is.Qid} returns an object of class \code{is.Qid}.
#' \describe{
#' \item{\code{strict}}{Is the Q-matrix strictly identifiable? (\code{logical}).}
#' \item{\code{generic}}{Is the Q-matrix generically identifiable?  (\code{logical}).}
#' \item{\code{conditions}}{Identifiability criteria and whether they are fulfilled or not (\code{vector}).}
#' \item{\code{specifications}}{Function call specifications (\code{list}).}
#' }
#'
#' @author {Pablo Nájera, Universidad Autónoma de Madrid \cr Miguel A. Sorrel, Universidad Autónoma de Madrid}
#'
#' @references
#' Gu, Y., & Xu, G. (2021). Sufficient and necessary conditions for the identifiability of the Q-matrix. \emph{Statistica Sinica}, \emph{31}, 449-472. https://www.jstor.org/stable/26969691
#'
#' @export
#'
#' @examples
#' Kj <- c(15, 10, 0, 5)
#' Q <- genQ(J = 30, K = 4, Kj = Kj, Qid = "others", seed = 123)$gen.Q
#' idQ <- is.Qid(Q, model = "DINA")
is.Qid <- function(Q, model){
  if(!is.matrix(Q) & !is.data.frame(Q)){stop("Error in is.Qid: Q must be a matrix or data.frame.")}
  if(!(model %in% c("others", "DINA", "DINO"))){stop("Error in is.Qid: model must be 'DINA', 'DINO', or 'others'.")}
  spec <- list(Q = Q, model = model)
  Q <- as.matrix(Q)
  K <- ncol(Q)
  J <- nrow(Q)
  rownames(Q) <- 1:J
  qI <- c()
  for(k in 1:ncol(Q)){
    tmp <- Q[Q[,k] == 1 & rowSums(Q) == 1, , drop = FALSE]
    if(nrow(tmp) == 0){next}
    qI <- c(qI, rownames(tmp[1, , drop = FALSE]))
  }
  Q1 <- Q[which(rowSums(Q) == 1), , drop = FALSE]
  QP <- Q[!rownames(Q) %in% qI,]
  I1 <- all(colSums(Q1) >= 1)
  I2 <- all(colSums(Q1) >= 2)
  I3 <- all(colSums(Q1) >= 3)
  Kj3 <- all(colSums(Q) >= 3)
  A <- I1 # Condition A: Completeness
  B <- nrow(unique(t(QP))) == K # Condition B: Distinctiveness
  C <- Kj3 # Condition C: Repetition
  if(A & B & C){
    res <- list(strict = TRUE, generic = TRUE, conditions = data.frame(completeness = A, distinctiveness = B, repetition = C), specifications = spec)
    class(res) <- "is.Qid"
    return(res)
  }
  if(model %in% c("DINA", "DINO") & !A){
    res <- list(strict = FALSE, generic = FALSE, conditions = data.frame(completeness = A, distinctiveness = B, repetition = C), specifications = spec)
    class(res) <- "is.Qid"
    return(res)
  }
  if(model %in% c("DINA", "DINO") & !B){
    res <- list(strict = FALSE, generic = FALSE, conditions = data.frame(completeness = A, distinctiveness = B, repetition = C), specifications = spec)
    class(res) <- "is.Qid"
    return(res)
  }
  if(model %in% c("DINA", "DINO") & !C){
    if(any(colSums(Q) == 1)){
      res <- list(strict = FALSE, generic = FALSE, conditions = data.frame(completeness = A, distinctiveness = B, repetition = C), specifications = spec)
      class(res) <- "is.Qid"
      return(res)
    }
    Ca <- any(rowSums(Q) == K)
    if(Ca){
      res <- list(strict = FALSE, generic = FALSE, conditions = data.frame(completeness = A, distinctiveness = B, repetition = C), specifications = spec)
      class(res) <- "is.Qid"
      return(res)
    }
    dupl.Q1 <- as.numeric(which(colSums(Q1[which(duplicated(Q1)),]) > 0))
    dupl.Q1 <- as.numeric(names(which(table(c(dupl.Q1, which(colSums(Q) == 2))) > 1)))
    Cb <- FALSE
    if(length(dupl.Q1) > 0){
      for(k in dupl.Q1){
        Qk <- Q[-as.numeric(rownames(Q1[Q1[,k] == 1,])), -k]
        qI <- c()
        for(k2 in 1:ncol(Qk)){qI <- c(qI, rownames(Qk[Qk[,k2] == 1 & rowSums(Qk) == 1, , drop = FALSE][1, , drop = FALSE]))}
        Qk1 <- Qk[which(rowSums(Qk) == 1), , drop = FALSE]
        QkP <- Qk[!rownames(Qk) %in% qI,]
        Cb1.A <- all(colSums(Qk1) >= 1)
        Cb1.B <- nrow(unique(t(QkP))) == (K - 1)
        Cb1.C <- all(colSums(Qk) >= 3)
        if(Cb1.A & Cb1.B & Cb1.C){
          res <- list(strict = FALSE, generic = TRUE, conditions = data.frame(completeness = A, distinctiveness = B, repetition = C), specifications = spec)
          class(res) <- "is.Qid"
          return(res)
        }
        Cb2 <- all(colSums(Qk1) >= 2)
        if(Cb2){
          res <- list(strict = FALSE, generic = TRUE, conditions = data.frame(completeness = A, distinctiveness = B, repetition = C), specifications = spec)
          class(res) <- "is.Qid"
          return(res)
        }
      }
    }
    uniq.Q1 <- as.numeric(which(colSums(unique(Q1)) == 1))
    Cc <- FALSE
    for(k in uniq.Q1){
      j1 <- as.numeric(rownames(Q1[Q1[,k] == 1, , drop = FALSE][1, , drop = FALSE]))
      jj2 <- as.numeric(rownames(Q[Q[,k] == 1, , drop = FALSE][-j1, , drop = FALSE]))
      for(j2 in jj2){
        Qjk <- Q[-c(j1, j2), -k]
        qI <- c()
        for(k2 in 1:ncol(Qjk)){qI <- c(qI, rownames(Qjk[Qjk[,k2] == 1 & rowSums(Qjk) == 1, , drop = FALSE][1, , drop = FALSE]))}
        Qjk1 <- Qjk[which(rowSums(Qjk) == 1), , drop = FALSE]
        QjkP <- Qjk[!rownames(Qjk) %in% qI,]
        Cc.A <- all(colSums(Qjk1) >= 1)
        Cc.B <- nrow(unique(t(QjkP))) == (K - 1)
        Cc.C <- all(colSums(Qjk) >= 3)
        if(Cc.A & Cc.B & Cc.C){
          res <- list(strict = FALSE, generic = TRUE, conditions = data.frame(completeness = A, distinctiveness = B, repetition = C), specifications = spec)
          class(res) <- "is.Qid"
          return(res)
        }
      }
    }
  } else {
    if(J <= K * 2){stop("The Q-matrix does not have enough items to achieve identifiability.")}
    combs1 <- utils::combn(J, K)
    for(cc1 in 1:ncol(combs1)){
      c1 <- combs1[,cc1]
      combs2 <- utils::combn(setdiff(1:J, c1), K)
      for(cc2 in 1:ncol(combs2)){
        c2 <- combs2[,cc2]
        if((min(c2) - max(c1)) <= 0){next}
        QI1 <- Q[c(c1),]
        QI2 <- fungible::faAlign(QI1, Q[c(c2),])$F2
        QI <- Q[c(c1, c2),]
        D <- all(colSums(QI) >= 2) & !identical(QI1, QI2) # Condition D: Generic completeness
        Qstar <- Q[-c(c1, c2),, drop = FALSE]
        E <- all(colSums(Qstar) >= 1) # Condition E: Generic repetition
        if(D & E){
          res <- list(strict = FALSE, generic = TRUE, conditions = data.frame(completeness = A, distinctiveness = B, repetition = C, generic.completeness = D, generic.repetition = E), specifications = spec)
          class(res) <- "is.Qid"
          return(res)
        }
      }
    }
  }
  if(model %in% c("DINA", "DINO")){
    res <- list(strict = FALSE, generic = FALSE, conditions = data.frame(completeness = A, distinctiveness = B, repetition = C), specifications = spec)
    class(res) <- "is.Qid"
    return(res)
  } else {
    res <- list(strict = FALSE, generic = FALSE, conditions = data.frame(completeness = A, distinctiveness = B, repetition = C, generic.completeness = D, generic.repetition = E), specifications = spec)
    class(res) <- "is.Qid"
    return(res)
  }
}
