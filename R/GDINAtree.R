#' @title G-DINA tree
#'
#' @description Implementation of a model-based recursive partitioning algorithm (Zeileis et al., 2008) for the G-DINA model (Nájera et al., in press) to detect subpopulations based on
#' different item parameters, underlying model, or Q-matrix. The function is based on the \code{mob} function of the package \code{partykit} package (Hothorn & Zeileis, 2015), and
#' estimates CDMs via the \emph{GDINA} package (MA & de la Torre, 2020).
#'
#' @param dat A \emph{N} individuals x \emph{J} items (\code{matrix} or \code{data.frame}). Missing values need to be coded as \code{NA}.
#' @param covariates A \emph{N} individuals x \emph{M} covariates (\code{matrix} or \code{data.frame}). Missing values need to be coded as \code{NA}.
#' @param Q A \emph{J} items x \emph{K} attributes Q-matrix (\code{matrix} or \code{data.frame}).
#' @param model CDM to be estimated. The possible options include "GDINA","DINA","DINO","ACDM","LLM", "RRUM", "MSDINA", "BUGDINO", "SISM", and "UDF". See the \code{GDINA} package (Ma & de la Torre, 2020) for more information. Default is \code{"GDINA"}.
#' @param maxdepth An \code{integer} indicating the maximum depth of a tree.
#' @param minsize An \code{integer} indicating the minimum number of observations in a node.
#' @param alpha Nominal significance level.
#' @param bonferroni A \code{logical} value indicating whether the Bonferroni correction should be applied.
#'
#' @return \code{GDINAtree} returns an object of class \code{GDINAtree}.
#' \describe{
#' \item{\code{tree}}{A \code{modelparty} object containing the main results from the G-DINA tree.}
#' \item{\code{fit.nodes}}{A \code{list} containing the models fitted in each node (each of them being a \code{GDINA} object).}
#' \item{\code{specifications}}{Function call specifications (\code{list}).}
#' }
#'
#' @author {David Goretzko, Goethe University Frankfurt, \cr Philipp Sterner, LMU Munich, \cr Pablo Nájera, Universidad Pontificia Comillas}
#'
#' @references
#' Hothorn, T., & Zeileis, A. (2015). partykit: A Modular Toolkit for Recursive Partytioning in R. \emph{Journal of Machine Learning Research}, \emph{16}, 3905-3909. https://jmlr.org/papers/v16/hothorn15a.html
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R package for cognitive diagnosis modeling. \emph{Journal of Statistical Software}, \emph{93}(14). https://doi.org/10.18637/jss.v093.i14
#'
#' Zeileis, A., Hothorn, T., & Hornik, K. (2008). Model-Based Recursive Partitioning. \emph{Journal of Computational and Graphical Statistics}, \emph{17}(2), 492-514. https://doi.org/10.1198/106186008X319331
#'
#' @export
#'
#' @examples
#' library(GDINA)
#' dat <- as.data.frame(rbind(sim30GDINA$simdat, sim30DINA$simdat))
#' Q <- sim30GDINA$simQ
#' set.seed(1713)
#' dat$y <- rep(1:2, each = nrow(dat)/2)
#' dat$cov1 <- round(runif(nrow(dat)), 2)
#' dat$cov2 <- round(runif(nrow(dat)), 2)
#' dat$cov3 <- sample(1:2, nrow(dat), replace = TRUE)
#' dat$cov4 <- sample(1:4, nrow(dat), replace = TRUE)
#' fit <- GDINAtree(dat = dat[,1:30], covariates = dat[,31:35], Q = Q)
#' plot(fit)
GDINAtree <- function(dat, covariates, Q, model = "GDINA", maxdepth = 3, minsize = 100, alpha = 0.05, bonferroni = TRUE, ...){

  #--------------------
  # Auxiliary function
  #--------------------

  GDINAtree_fit <- function(q){
    function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ..., estfun = TRUE, object = FALSE){
      mod <- GDINA::GDINA(dat = y, Q = q, model = model, start = start)
      list(coefficients = unlist(mod$catprob.parm), objfun = -as.numeric(logLik(mod)),
           estfun = do.call(cbind, GDINA::score(mod, parm = "prob")[-length(GDINA::score(mod))]), object = NULL)
    }
  }

  #-------------
  # G-DINA tree
  #-------------

  tree <- partykit::mob(formula = as.formula(paste(paste(colnames(dat), collapse = "+"), paste("~"), paste(colnames(covariates), collapse = "+"))),
                        data = cbind(dat, covariates),
                        fit = GDINAtree_fit(Q),
                        control = partykit::mob_control(ytype = "data.frame", alpha = alpha, bonferroni = bonferroni, minsize = minsize, maxdepth = maxdepth, vcov = "opg", ...))
  n.nodes <- length(partykit::nodeids(tree))
  fit.nodes <- list()
  for(n in 1:n.nodes){
    dat.n <- dat[as.numeric(names(predict(tree[n], type = "node"))),]
    fit.nodes[[n]] <- GDINA::GDINA(dat = dat.n, Q = Q, model = model, verbose = 0)
  }

  #-----------------
  # Return outcomes
  #-----------------

  specs <- list(dat = dat, covariates = covariates, Q = Q, model = model,
                maxdepth = maxdepth, minsize = minsize, alpha = alpha, bonferroni = bonferroni)
  res <- list(tree = tree, fit.nodes = fit.nodes, specifications = specs)
  class(res) <- "GDINAtree"
  return(res)
}
