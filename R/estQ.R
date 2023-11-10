#' @title Empirical Q-matrix estimation
#'
#' @description Empirical Q-matrix estimation based on the \emph{discrete factor loading} method (Wang, Song, & Ding, 2018) as used in Nájera, Abad, and Sorrel (2021).
#' Apart from the conventional dichotomization criteria, the procedure based on loading differences described in Garcia-Garzon, Abad, and Garrido (2018) is also available.
#' Furthermore, the bagging bootstrap implementation (Xu & Shang, 2018) can be applied; it is recommended when working with small sample sizes.
#' The \code{psych} package (Revelle, 2020) is used for estimating the required exploratory factor analysis (EFA).
#'
#' @param r A correlation matrix or raw data (\code{matrix} or \code{data.frame}). If a correlation matrix is used, it must have dimensions \emph{J} items × \emph{J} items. Please note that tetrachoric or polychoric correlations should be used when working with dichotomous or polytomous items, respectively. If raw data is used, it must have dimensions \emph{N} individuals × \emph{J} items. Missing values need to be coded as \code{NA}.
#' @param K Number of attributes to use.
#' @param n.obs Number of individuals if \code{r} is a correlation matrix. If \code{n.obs} is provided, \code{r} will be treated as a correlation matrix. Use \code{NULL} if \code{r} is raw data. The default is \code{NULL}.
#' @param criterion Dichotomization criterion to transform the factor loading matrix into the Q-matrix. The possible options include \code{"row"} (for row means), \code{"col"} (for column means), \code{"loaddiff"} (for the procedure based on loading differences), or a value between 0 and 1 (for a specific threshold). The default is \code{"row"}.
#' @param boot Apply the bagging bootstrap implementation? Only available if \code{r} is raw data. If \code{FALSE}, the EFA will be applied once using the whole sample size. If \code{TRUE}, several EFAs will be applied with different subsamples; the estimated Q-matrix will be dichotomized from the bootstrapped Q-matrix, but the EFA fit indices, factor loadings, and communalities will be computed from the EFA with the whole sample size. The default is \code{FALSE}.
#' @param efa.args A list of arguments for the EFA estimation:
#' \describe{
#' \item{\code{cor}}{Type of correlations to use. It includes \code{"cor"} (for Pearson correlations) and \code{"tet"} (for tetrachoric/polychoric correlations), among others. See \code{fa} function from the \code{psych} R package for additional details. The default is \code{"tet"}.}
#' \item{\code{rotation}}{Rotation procedure to use. It includes \code{"oblimin"}, \code{"varimax"}, and \code{"promax"}, among others. An oblique rotation procedure is usually recommended. See \code{fa} function from the \code{psych} R package for additional details. The default is \code{"oblimin"}.}
#' \item{\code{fm}}{Factoring method to use. It includes \code{"uls"} (for unweighted least squares), \code{"ml"} (for maximum likelihood), and \code{"wls"} (for weighted least squares), among others. See \code{fa} function from the \code{psych} R package for additional details. The default is \code{"uls"}.}
#' }
#' @param boot.args A list of arguments for the bagging bootstrap implementation (ignored if \code{boot = FALSE}):
#' \describe{
#' \item{\code{N}}{Sample size (or proportion of the total sample size, if lower than 1) to use in each bootstrap replication. The default is .8.}
#' \item{\code{R}}{Number of bootstrap replications. The default is 100.}
#' \item{\code{verbose}}{Show progress? The default is \code{TRUE}.}
#' \item{\code{seed}}{A seed for obtaining consistent results. If \code{NULL}, no seed is used. The default is \code{NULL}.}
#' }
#'
#' @return \code{estQ} returns an object of class \code{estQ}.
#' \describe{
#' \item{\code{est.Q}}{Estimated Q-matrix (\code{matrix}).}
#' \item{\code{efa.loads}}{Factor loading matrix (\code{matrix}).}
#' \item{\code{efa.comm}}{EFA communalities (\code{vector}).}
#' \item{\code{efa.fit}}{EFA model fit indices (\code{vector}).}
#' \item{\code{boot.Q}}{Bagging bootstrap Q-matrix before dichotomization. Only if \code{boot = TRUE} (\code{matrix}).}
#' \item{\code{is.Qid}}{Q-matrix identifiability information (\code{list}).}
#' \item{\code{specifications}}{Function call specifications (\code{list}).}
#' }
#'
#' @author {Pablo Nájera, Universidad Autónoma de Madrid}
#'
#' @references
#' Garcia-Garzon, E., Abad, F. J., & Garrido, L. E. (2018). Improving bi-factor exploratory modelling: Empirical target rotation based on loading differences. \emph{Methodology}, \emph{15}, 45–55. https://doi.org/10.1027/1614-2241/a000163
#'
#' Nájera, P., Abad, F. J., & Sorrel, M. A. (2021). Determining the number of attributes in cognitive diagnosis modeling. \emph{Frontiers in Psychology}, \emph{12}:614470. https://doi.org/10.3389/fpsyg.2021.614470
#'
#' Revelle, W. (2019). \emph{psych: Procedures for Psychological, Psychometric, and Personality Research}. R package version 1.9.12. https://CRAN.R-project.org/package=psych.
#'
#' Wang, W., Song, L., & Ding, S. (2018). An exploratory discrete factor loading method for Q-matrix specification in cognitive diagnosis models. In: M. Wilberg, S. Culpepper, R. Janssen, J. Gonzalez, & D. Molenaar (Eds.), \emph{Quantitative Psychology. IMPS 2017. Springer Proceedings in Mathematics & Statistics} (Vol. 233, pp. 351–362). Springer.
#'
#' Xu, G., & Shang, Z. (2018). Identifying latent structures in restricted latent class models. \emph{Journal of the American Statistical Association}, \emph{113}, 1284–1295. https://doi.org/10.1080/01621459.2017.1340889
#'
#' @export
#'
#' @examples
#' library(GDINA)
#' dat <- sim30GDINA$simdat
#' Q <- sim30GDINA$simQ
#'
#' #------------------------------
#' # Using default specifications
#' #------------------------------
#' sugQ1 <- estQ(r = dat, K = 5) # Estimate Q-matrix
#' sugQ1$est.Q <- orderQ(sugQ1$est.Q, Q)$order.Q # Reorder Q-matrix attributes
#' mean(sugQ1$est.Q == Q) # Check similarity with the generating Q-matrix
#'
#' #------------------------------------
#' # Using the bagging bootstrap method
#' #------------------------------------
#' # In boot.args argument, R >= 100 is recommended (R = 20 is here used for illustration purposes)
#' sugQ2 <- estQ(r = dat, K = 5, boot = TRUE, boot.args = list(R = 20, seed = 123)) # Estimate Q-matrix
#' sugQ2$est.Q <- orderQ(sugQ2$est.Q, Q)$order.Q # Reorder Q-matrix attributes
#' sugQ2$boot.Q # Proportion of replicas a q-entry was specified in the estimated Q-matrix
#' mean(sugQ2$est.Q == Q) # Check similarity with the generating Q-matrix
estQ <- function(r, K, n.obs = NULL, criterion = "row", boot = FALSE, efa.args = list(cor = "tet", rotation = "oblimin", fm = "uls"), boot.args = list(N = .8, R = 100, verbose = TRUE, seed = NULL)){
  if(!is.matrix(r) & !is.data.frame(r)){stop("Error in estQ: r must be a matrix or data.frame.")}
  if((!is.numeric(K) & !is.double(K)) | length(K) > 1){stop("Error in estQ: K must be a unique numeric value.")}
  if(K < 1){stop("Error in estQ: K must be greater than 0.")}
  if(!is.null(n.obs) & !is.numeric(n.obs) & !is.double(n.obs) | length(n.obs) > 1){stop("Error in estQ: n.obs must be NULL or a unique numeric value.")}
  N <- ifelse(is.null(n.obs), nrow(r), n.obs)
  if(!is.numeric(criterion) & !is.double(criterion) | length(criterion) > 1){
    if((criterion > 1 | criterion < 0) & !(criterion %in% c("row", "col", "loaddiff"))){stop("Error in estQ: criterion must be 'row', 'col', 'loaddiff', or a value between 0 and 1.")}
  }
  if(!is.logical(boot) | length(boot) > 1){stop("Error in estQ: boot must be logical.")}
  doboot <- boot
  if(!is.null(n.obs) & boot){
    warning("Warning in estQ: The bagging bootstrap method can be applied only when r is raw data (n.obs is NULL).")
    doboot <- FALSE
  } else if(is.null(n.obs) & boot){
    doboot <- TRUE
  }
  if(is.null(efa.args$cor)){efa.args$cor <- "tet"}
  if(is.null(efa.args$rotation)){efa.args$cor <- "oblimin"}
  if(is.null(efa.args$fm)){efa.args$cor <- "uls"}
  if(!(efa.args$cor %in% c("tet", "cor"))){stop("Error in estQ: efa.args$cor must be 'cor' or 'tet'.")}
  if(!(efa.args$rotation %in% c("none", "varimax", "quartimax", "bentlerT", "equamax", "varimin", "geominT", "bifactor", "Promax", "promax", "oblimin", "simplimax", "bentlerQ", "geominQ", "biquartimin", "cluster"))){stop("Error in estQ: efa.args$rotation must be 'none', 'varimax', 'quartimax', 'bentlerT', 'equamax', 'varimin', 'geominT', 'bifactor', 'Promax', 'promax', 'oblimin', 'simplimax', 'bentlerQ', 'geominQ', 'biquartimin', or 'cluster'.")}
  if(!(efa.args$fm %in% c("minres", "uls", "ols", "wls", "gls", "pa", "ml", "minchi", "minrank", "old.min", "alpha"))){stop("Error in estQ: efa.args$fm must be 'minres', 'uls', 'ols', 'wls', 'gls', 'pa', 'ml', 'minchi', 'minrank', 'old.min', or 'alpha'.")}
  if(is.null(boot.args$N)){boot.args$N <- 0.8}
  if(is.null(boot.args$R)){boot.args$R <- 100}
  if(is.null(boot.args$verbose)){boot.args$verbose <- TRUE}
  if((!is.numeric(boot.args$N) & !is.double(boot.args$N)) | length(boot.args$N) > 1){stop("Error in estQ: boot.args$N must be a unique numeric value.")}
  if(boot.args$N < 1){
    boot.N <- "prop"
    if(boot.args$N > 1 | boot.args$N < 0){stop("Error in estQ: boot.args$N must be between 0 and 1 when a proportion is provided.")}
  } else {
    boot.N <- "num"
    if(boot.args$N > N){stop("Error in estQ: boot.args$N must be lower than the number of observations when a number is provided.")}
  }
  if((!is.numeric(boot.args$R) & !is.double(boot.args$R)) | length(boot.args$R) > 1){stop("Error in estQ: boot.args$R must be a unique numeric value.")}
  if(!is.logical(boot.args$verbose)){stop("Error in estQ: boot.args$verbose must be logical.")}
  if(!is.null(boot.args$seed)){if((!is.numeric(boot.args$seed) & !is.double(boot.args$seed)) | length(boot.args$seed) > 1){stop("Error in estQ: boot.args$seed must be a unique numeric value.")}}

  J <- ncol(r)
  boot.Q <- NULL
  if(!is.null(boot.args$seed)){set.seed(boot.args$seed)}

  if(is.null(n.obs)){
    if(efa.args$cor == "cor"){
      use.r <- stats::cor(r)
    } else if(efa.args$cor == "tet"){
      use.r <- sirt::tetrachoric2(r, progress = FALSE)$rho
    }
  } else {
    use.r <- r
  }

  efa <- cdmTools.fa(r = use.r, n.obs = N, nfactors = K, rotate = efa.args$rotation, cor = efa.args$cor, fm = efa.args$fm)
  efa.comm <- efa$communality
  efa.loads <- round(matrix(efa$loadings, nrow = J, ncol = K), 3)

  Vaccounted <- ifelse(K == 1, efa$Vaccounted[2,], efa$Vaccounted[3, K])
  efa.fit <- c(NumFactors = K, ChiSq = efa$STATISTIC, ChiSq.PVAL = efa$PVAL, TLI = efa$TLI, efa$RMSEA[1], BIC = efa$BIC, MeanComm = mean(efa$communality), VarAccounted = Vaccounted, dof = efa$dof)

  if(K == 1){
    est.Q <- matrix(1, nrow = J)
  } else {
    if(doboot == FALSE){
      if(criterion == "row"){
        est.Q <- t(apply(efa.loads, 1, function(x) as.numeric(x >= mean(x))))
      } else if(criterion == "col"){
        est.Q <- apply(efa.loads, 2, function(x) as.numeric(x >= mean(x)))
      } else if(criterion == "loaddiff"){
        A <- efa.loads
        hA <- efa.comm
        B <- (t(sapply(1:nrow(A), function(j) A[j,] * as.vector(solve(sqrt(hA[j]))))))^2
        B.j <- matrix(1:nrow(B), nrow = nrow(B), ncol = ncol(B))
        D.j <- apply(B, 2, function(f) order(f, decreasing = T))
        C <- apply(B, 2, function(f) sort(f, decreasing = T))
        D <- apply(C, 2, function(f) abs(diff(f)))
        jump <- sapply(1:ncol(D), function(f) utils::tail(D[D[,f] > mean(D[,f]), f], 1))
        cut <- abs(sapply(1:ncol(D), function(f) A[D.j[which(D[,f] == jump[f]), f], f]))
        est.Q <- sapply(1:ncol(A), function(f) ifelse(abs(A[,f]) >= cut[f], 1, 0))
      } else {
        est.Q <- matrix(as.numeric(efa.loads >= criterion), nrow = J, ncol = K)
      }
    } else if(doboot == TRUE){
      boot.Q <- matrix(0, nrow = J, ncol = K)
      for(i in 1:boot.args$R){
        if(boot.N == "prop"){
          N.boot <- sort(sample(1:N, round(N * boot.args$N), replace = F))
        } else if(boot.N == "num"){
          N.boot <- sort(sample(1:N, boot.args$N, replace = F))
        }
        if(efa.args$cor == "cor"){
          use.r <- stats::cor(r[N.boot,])
        } else if(efa.args$cor == "tet"){
          use.r <- sirt::tetrachoric2(r[N.boot,], progress = FALSE)$rho
        }
        efa.boot <- cdmTools.fa(r = use.r, n.obs = N, nfactors = K, rotate = efa.args$rotation, cor = efa.args$cor, fm = efa.args$fm)
        efa.loads.boot <- round(matrix(efa.boot$loadings, nrow = J, ncol = K), 3)

        if(criterion == "row"){
          tmp.Q <- t(apply(efa.loads.boot, 1, function(x) as.numeric(x >= mean(x))))
          tmp.Q <- orderQ(tmp.Q, boot.Q)$order.Q
          boot.Q <- boot.Q + tmp.Q
        } else if(criterion == "col"){
          tmp.Q <- apply(efa.loads.boot, 2, function(x) as.numeric(x >= mean(x)))
          tmp.Q <- orderQ(tmp.Q, boot.Q)$order.Q
          boot.Q <- boot.Q + tmp.Q
        } else if(criterion == "loaddiff"){
          A <- efa.loads.boot
          hA <- efa.comm
          B <- (t(sapply(1:nrow(A), function(j) A[j,] * as.vector(solve(sqrt(hA[j]))))))^2
          B.j <- matrix(1:nrow(B), nrow = nrow(B), ncol = ncol(B))
          D.j <- apply(B, 2, function(f) order(f, decreasing = T))
          C <- apply(B, 2, function(f) sort(f, decreasing = T))
          D <- apply(C, 2, function(f) abs(diff(f)))
          jump <- sapply(1:ncol(D), function(f) utils::tail(D[D[,f] > mean(D[,f]), f], 1))
          cut <- abs(sapply(1:ncol(D), function(f) A[D.j[which(D[,f] == jump[f]), f], f]))
          tmp.Q <- sapply(1:ncol(A), function(f) ifelse(abs(A[,f]) >= cut[f], 1, 0))
          tmp.Q <- orderQ(tmp.Q, boot.Q)$order.Q
          boot.Q <- boot.Q + tmp.Q
        } else {
          tmp.Q <- matrix(as.numeric(efa.loads.boot >= criterion), nrow = J, ncol = K)
          tmp.Q <- orderQ(tmp.Q, boot.Q)$order.Q
          boot.Q <- boot.Q + tmp.Q
        }
        if(boot.args$verbose){cat("\r", "In estQ: Iteration", i, "out of", boot.args$R)}
      }
      boot.Q <- boot.Q / boot.args$R
      est.Q <- ifelse(boot.Q >= 0.5, 1, 0)
    }
  }

  if(any(rowSums(est.Q) == 0)){warning("Warning in estQ: Item(s) {", which(rowSums(est.Q) == 0), "} are not measuring any attributes. Check factor loadings for the item(s).")}
  if(any(colSums(est.Q) == 0)){warning("Warning in estQ: Attribute(s) {", which(colSums(est.Q) == 0), "} are not measured by any items. Check factor loadings for the attribute(s) or consider reducing K.")}

  idQ.DINA <- is.Qid(est.Q, model = "DINA")
  idQ.others <- is.Qid(est.Q, model = "others")

  spec <- list(r = r, K = K, n.obs = n.obs, criterion = criterion, efa.args = list(cor = efa.args$cor, rotation = efa.args$rotation, fm = efa.args$fm))
  res <- list(est.Q = est.Q, efa.loads = efa.loads, efa.comm = efa.comm, efa.fit = efa.fit, boot.Q = boot.Q, is.Qid = list("DINA" = idQ.DINA, "others" = idQ.others), specifications = spec)
  class(res) <- "estQ"
  return(res)
}
