#' @title Check whether a Q-matrix is identifiable
#'
#' @description Checks whether a Q-matrix is complete (Köhn & Chiu, 2017, 2018) and identifiable according to the criteria from Chen, Liu, Xu, & Ying (2015) and Xu & Shang (2018).
#'
#' @param Q A \emph{J} items x \emph{K} attributes Q-matrix (\code{matrix} or \code{data.frame}).
#' @param model CDM to be considered. It includes \code{"DINA"}, \code{"DINO"}, or \code{"others"} (for other CDMs; e.g., G-DINA, A-CDM). The default is \code{"others"}.
#' @param verbose Should a message about the identifiability of the Q-matrix be printed? The default is \code{TRUE}.
#'
#' @return \code{is.Qid} returns an object of class \code{is.Qid}.
#' \describe{
#' \item{\code{id.Q}}{Is the Q-matrix identifiable? (\code{logical}).}
#' \item{\code{comp.Q}}{Is the Q-matrix complete?  (\code{logical}).}
#' \item{\code{criteria.Qid}}{Identifiability criteria and whether they are fulfilled or not (\code{vector}).}
#' \item{\code{message}}{A message about the identifiability of the Q-matrix and references (\code{string}).}
#' \item{\code{specifications}}{Function call specifications (\code{list}).}
#' }
#'
#' @author {Pablo Nájera, Universidad Autónoma de Madrid \cr Miguel A. Sorrel, Universidad Autónoma de Madrid}
#'
#' @references
#' Chen, Y., Liu, J., Xu, G., & Ying, Z. (2015). Statistical analysis of Q-matrix based diagnostic classification models. \emph{Journal of the American Statistical Association}, \emph{110}, 850-866. https://doi.org/10.1080/01621459.2014.934827
#'
#' Köhn, H.-F., & Chiu, C.-Y. (2017). A procedure for assessing the completeness of the Q-matrices of cognitively diagnostic tests. \emph{Psychometrika}, \emph{82}, 112-132. https://doi.org/10.1007/s11336-016-9536-7
#'
#' Köhn, H.-F., & Chiu, C.-Y. (2018). How to build a complete Q-matrix for a cognitively diagnostic test. \emph{Journal of Classification}, \emph{35}, 273-299. https://doi.org/10.1007/s00357-018-9255-0
#'
#' Xu, G., & Shang, Z. (2018). Identifying latent structures in restricted latent class models. \emph{Journal of the American Statistical Association}, \emph{113}, 1284-1295. https://doi.org/10.1080/01621459.2017.1340889
#'
#' @export
#'
#' @examples
#' Kj <- c(15, 10, 0, 5)
#' Q <- genQ(J = 30, K = 4, Kj = Kj, Qid = "others", seed = 123)$gen.Q
#' idQ <- is.Qid(Q)
is.Qid <- function(Q, model = "others", verbose = TRUE){
  if(!is.matrix(Q) & !is.data.frame(Q)){stop("Error in is.Qid: Q must be a matrix or data.frame.")}
  if(!(model %in% c("others", "DINA", "DINO"))){stop("Error in is.Qid: model must be 'DINA', 'DINO', or 'others'.")}
  if(!is.logical(verbose)){stop("Error in is.Qid: verbose must be logical.")}

  K <- ncol(Q)
  J <- nrow(Q)
  Q1 <- Q[which(rowSums(Q) == 1),,drop = FALSE]
  I1 <- all(colSums(Q1) >= 1)
  comp.Q <- ifelse(I1, TRUE, FALSE)
  I2 <- all(colSums(Q1) >= 2)
  I3 <- all(colSums(Q1) >= 3)
  Kj2 <- all(colSums(Q) >= 2)
  Kj3 <- all(colSums(Q) >= 3)
  if(model == "DINA" | model == "DINO"){
    id.Q <- ifelse(I2 & Kj3, TRUE, FALSE)
    info <- ifelse(id.Q, "The Q-matrix will be identifiable if (a) all latent classes have a non-zero probability and (b) all items have discriminating power. See Chen, Liu, Xu, & Ying (2015) for further details.", "The Q-matrix is not identifiable.")
  } else if(model == "others"){
    id.Q <- ifelse(I3, TRUE, FALSE)
    info <- ifelse(id.Q, "The Q-matrix will be identifiable if (a) all latent classes have a non-zero probability and (b) the monotonicity constraint is fulfilled. See Xu & Shang (2018) for further details.", "The Q-matrix is not identifiable.")
  }
  criteria.Qid <- c(I1, I2, I3, Kj2, Kj3)
  names(criteria.Qid) <- c("1 identity matrix", "2 identity matrices", "3 identity matrices", "2 items per attribute", "3 items per attribute")

  if(verbose){cat("\r", "In is.Qid:", info, "\n")}
  spec <- list(Q = Q, model = model, verbose = verbose)
  res <- list(id.Q = id.Q, comp.Q = comp.Q, criteria.Qid = criteria.Qid, message = info, specifications = spec)
  class(res) <- "is.Qid"
  return(res)
}
