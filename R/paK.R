#' @title Parallel analysis - dimensionality assessment method
#'
#' @description Parallel analysis with column permutation (i.e., resampling) as used in Nájera, Abad, & Sorrel (2021).
#' It is recommended to use principal components, Pearson correlations, and mean criterion (Garrido, Abad, & Ponsoda, 2013; Nájera, Abad, & Sorrel, 2021).
#' The parallel analysis based on principal axis factor analysis is conducted using the \code{fa.parallel} function of the \code{psych} R package (Revelle, 2020).
#' The tetrachoric correlations are efficiently estimated using the \code{sirt} R package (Robitzsch, 2020).
#' The graph is made with the \code{ggplot2} package (Wickham et al., 2020).
#'
#' @param dat A \emph{N} individuals x \emph{J} items (\code{matrix} or \code{data.frame}). Missing values need to be coded as \code{NA}.
#' @param R Number of resampled datasets (i.e., replications) to generate. The default is 100.
#' @param fa Extraction method to use. It includes \code{"pc"} (for principal components analysis), \code{"fa"} (for principal axis factor analysis), and \code{"both"}. The default is \code{"pc"}.
#' @param cor What type of correlations to use. It includes \code{"cor"} (for Pearson correlations), \code{"tet"} (for tetrachoric/polychoric correlations), and \code{"both"}. The default is \code{"both"}.
#' @param cutoff What criterion to use as the cutoff. It can be \code{"mean"} (for the average generated eigenvalues) or a value between 0 and 100 (for a percentile). A vector with several criteria can be used. The default is \code{"mean"}.
#' @param fm Factoring method to use. It includes \code{"uls"} (for unweighted least squares), \code{"ml"} (for maximum likelihood), and \code{"wls"} (for weighted least squares), among others. The default is \code{"uls"}.
#' @param plot Print the parallel analysis plot? Note that the plot might be messy if many variants are requested. The default is \code{TRUE}.
#' @param verbose progress. The default is \code{TRUE}.
#' @param seed A seed for obtaining consistent results. If \code{NULL}, no seed is used. The default is \code{NULL}.
#'
#' @return \code{paK} returns an object of class \code{paK}.
#' \describe{
#' \item{\code{sug.K}}{The suggested number of attributes for each variant (\code{vector}).}
#' \item{\code{e.values}}{The sample and reference eigenvalues (\code{matrix}).}
#' \item{\code{plot}}{The parallel analysis plot. Only if \code{plot = TRUE} (\code{plot}).}
#' \item{\code{specifications}}{Function call specifications (\code{list}).}
#' }
#'
#' @author {Pablo Nájera, Universidad Autónoma de Madrid \cr Miguel A. Sorrel, Universidad Autónoma de Madrid \cr Francisco J. Abad, Universidad Autónoma de Madrid}
#'
#' @references
#' Garrido, L. E., Abad, F. J., & Ponsoda, V. (2013). A new look at Horn's parallel analysis with ordinal variables. \emph{Psychological Methods}, \emph{18}, 454-474. https://doi.org/10.1037/a0030005
#'
#' Nájera, P., Abad, F. J., & Sorrel, M. A. (2021). Determining the number of attributes in cognitive diagnosis modeling. \emph{Frontiers in Psychology}, \emph{12}:614470. https://doi.org/10.3389/fpsyg.2021.614470
#'
#' Revelle, W. (2019). \emph{psych: Procedures for Psychological, Psychometric, and Personality Research}. R package version 1.9.12. https://CRAN.R-project.org/package=psych.
#'
#' Robitzsch, A. (2020). \emph{sirt: Supplementary Item Response Theory Models}. R package version 3.9-4. https://CRAN.R-project.org/package=sirt.
#'
#' Wickham, H., et al. (2020). \emph{ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics}. R package version 3.3.2. https://CRAN.R-project.org/package=ggplot2.
#'
#' @export
#'
#' @examples
#' library(GDINA)
#' dat <- sim30GDINA$simdat
#' Q <- sim30GDINA$simQ
#' # In paK, R = 100 is recommended (R = 30 is here used for illustration purposes)
#' pa.K <- paK(dat = dat, R = 30, fa = "pc", cutoff = c("mean", 95), plot = TRUE, seed = 123)
#' pa.K$sug.K # Check suggested number of attributes by each parallel analysis variant
#' pa.K$e.values # Check eigenvalues
#' pa.K$plot # Show parallel analysis plot
paK <- function(dat, R = 100, fa = "pc", cor = "both", cutoff = "mean", fm = "uls", plot = TRUE, verbose = TRUE, seed = NULL){
  if(!is.matrix(dat) & !is.data.frame(dat)){stop("Error in paK: dat must be a matrix or data.frame.")}
  if((!is.numeric(R) & !is.double(R)) | length(R) > 1){stop("Error in paK: R must be a unique numeric value.")}
  if(R < 1){stop("Error in paK: R must be greater than 0.")}
  if(!(fa %in% c("pc", "fa", "both"))){stop("Error in paK: fa must be 'pc', 'fa', or 'both'.")}
  if(!(cor %in% c("tet", "cor", "both"))){stop("Error in paK: cor must be 'cor', 'tet', or 'both'.")}
  for(i in cutoff){
    if(suppressWarnings(!is.na(as.numeric(i)))){i <- as.numeric(i)}
    if((i > 100 | i < 0) & !(i %in% c("mean"))){stop("Error in paK: cutoff must be 'mean' or a value between 0 and 100.")}
    if(i < 1){warning("Warning in paK: Values lower than 1 found in cutoff. Note that the percentiles are requested by using integers (e.g., 95 = percentile 95).")}
  }
  if(!(fm %in% c("minres", "uls", "ols", "wls", "gls", "pa", "ml", "minchi", "minrank", "old.min", "alpha"))){stop("Error in paK: fm must be 'minres', 'uls', 'ols', 'wls', 'gls', 'pa', 'ml', 'minchi', 'minrank', 'old.min', or 'alpha'.")}
  if(!is.logical(verbose)){stop("Error in paK: verbose must be logical.")}
  if(!is.null(seed)){if((!is.numeric(seed) & !is.double(seed)) | length(seed) > 1){stop("Error in paK: seed must be a unique numeric value.")}}
  if(sum(is.na(cor(dat, use = "pair"))) > 0) stop("Error in paK: Parallel analysis cannot be computed when NAs are found in dat correlation matrix.")

  if(!is.null(seed)){set.seed(seed)}
  J <- ncol(dat)
  N <- nrow(dat)

  if(fa == "both"){fa <- c("pc", "fa")}
  if(cor == "both"){cor <- c("cor", "tet")}

  resample.pc.r <- resample.fa.r <- resample.pc.p <- resample.fa.p <- matrix(NA, R, J)
  for(i in 1:R){
    resample <- sapply(1:J, function(j) sample(dat[,j], length(dat[,j]), F))
    if("cor" %in% cor){
      tmp.r <- cdmTools.fa(resample, nfactors = 1, cor = "cor", fm = fm)
      if("fa" %in% fa){resample.fa.r[i,] <- tmp.r$values}
      if("pc" %in% fa){resample.pc.r[i,] <- tmp.r$e.values}
    }
    if("tet" %in% cor){
      tmp.p <- cdmTools.fa(resample, nfactors = 1, cor = "tet", fm = fm)
      if("fa" %in% fa){resample.fa.p[i,] <- tmp.p$values}
      if("pc" %in% fa){resample.pc.p[i,] <- tmp.p$e.values}
    }
    if(verbose){cat("\r", "In paK: Iteration", i, "out of", R)}
  }
  if(verbose){cat("\n")}

  conds <- expand.grid(fa = fa, cor = cor, cutoff = cutoff)
  N.conds <- nrow(conds)
  reference <- as.data.frame(matrix(NA, nrow = N.conds, ncol = J, dimnames = list(apply(conds, 1, paste, collapse = "."), 1:J)))
  if("cor" %in% cor){
    if("fa" %in% fa){
      for(i in cutoff){
        if(i == "mean"){
          reference["fa.cor.mean",] <- colMeans(resample.fa.r)
        } else {
          reference[paste(c("fa.cor", i), collapse = "."),] <- apply(resample.fa.r, 2, function(x) stats::quantile(x, as.numeric(i) / 100, na.rm = T))
        }
      }
    }
    if("pc" %in% fa){
      for(i in cutoff){
        if(i == "mean"){
          reference["pc.cor.mean",] <- colMeans(resample.pc.r)
        } else {
          reference[paste(c("pc.cor", i), collapse = "."),] <- apply(resample.pc.r, 2, function(x) stats::quantile(x, as.numeric(i) / 100, na.rm = T))
        }
      }
    }
  }
  if("tet" %in% cor){
    if("fa" %in% fa){
      for(i in cutoff){
        if(i == "mean"){
          reference["fa.tet.mean",] <- colMeans(resample.fa.p)
        } else {
          reference[paste(c("fa.tet", i), collapse = "."),] <- apply(resample.fa.p, 2, function(x) stats::quantile(x, as.numeric(i) / 100, na.rm = T))
        }
      }
    }
    if("pc" %in% fa){
      for(i in cutoff){
        if(i == "mean"){
          reference["pc.tet.mean",] <- colMeans(resample.pc.p)
        } else {
          reference[paste(c("pc.tet", i), collapse = "."),] <- apply(resample.pc.p, 2, function(x) stats::quantile(x, as.numeric(i) / 100, na.rm = T))
        }
      }
    }
  }

  if("cor" %in% cor){cor.matrix <- cor(dat, use = "pair")}
  if("tet" %in% cor){tet.matrix <- cdmTools.tetrachoric2(dat, method = "Bo")$rho}
  dat.eigen <- as.data.frame(matrix(NA, nrow = 4, ncol = J, dimnames = list(c("dat.fa.cor", "dat.fa.tet", "dat.pc.cor", "dat.pc.tet"), 1:J)))
  if("fa" %in% fa){
    if("cor" %in% cor){dat.eigen["dat.fa.cor",] <- cdmTools.fa(dat, nfactors = 1, cor = "cor", fm = fm)$values}
    if("tet" %in% cor){dat.eigen["dat.fa.tet",] <- cdmTools.fa(dat, nfactors = 1, cor = "tet", fm = fm)$values}
  }
  if("pc" %in% fa){
    if("cor" %in% cor){dat.eigen["dat.pc.cor",] <- eigen(cor.matrix, only.values = T)$values}
    if("tet" %in% cor){dat.eigen["dat.pc.tet",] <- eigen(tet.matrix, only.values = T)$values}
  }
  dat.eigen <- stats::na.omit(dat.eigen)

  nF <- as.data.frame(matrix(NA, nrow = 1, ncol = N.conds, dimnames = list(1, apply(conds, 1, paste, collapse = "."))))
  if("cor" %in% cor){
    if("fa" %in% fa){
      for(i in cutoff){
        if(i == "mean"){
          nF[,"fa.cor.mean"] <- which(dat.eigen["dat.fa.cor",] - as.numeric(reference["fa.cor.mean",]) < 0)[1] - 1
        } else {
          nF[,paste(c("fa.cor", i), collapse = ".")] <- which(dat.eigen["dat.fa.cor",] - as.numeric(reference[paste(c("fa.cor", i), collapse = "."),]) < 0)[1] - 1
        }
      }
    }
    if("pc" %in% fa){
      for(i in cutoff){
        if(i == "mean"){
          nF[,"pc.cor.mean"] <- which(dat.eigen["dat.pc.cor",] - as.numeric(reference["pc.cor.mean",]) < 0)[1] - 1
        } else {
          nF[,paste(c("pc.cor", i), collapse = ".")] <- which(dat.eigen["dat.pc.cor",] - as.numeric(reference[paste(c("pc.cor", i), collapse = "."),]) < 0)[1] - 1
        }
      }
    }
  }
  if("tet" %in% cor){
    if("fa" %in% fa){
      for(i in cutoff){
        if(i == "mean"){
          nF[,"fa.tet.mean"] <- which(dat.eigen["dat.fa.tet",] - as.numeric(reference["fa.tet.mean",]) < 0)[1] - 1
        } else {
          nF[,paste(c("fa.tet", i), collapse = ".")] <- which(dat.eigen["dat.fa.tet",] - as.numeric(reference[paste(c("fa.tet", i), collapse = "."),]) < 0)[1] - 1
        }
      }
    }
    if("pc" %in% fa){
      for(i in cutoff){
        if(i == "mean"){
          nF[,"pc.tet.mean"] <- which(dat.eigen["dat.pc.tet",] - as.numeric(reference["pc.tet.mean",]) < 0)[1] - 1
        } else {
          nF[,paste(c("pc.tet", i), collapse = ".")] <- which(dat.eigen["dat.pc.tet",] - as.numeric(reference[paste(c("pc.tet", i), collapse = "."),]) < 0)[1] - 1
        }
      }
    }
  }
  nF <- nF[,!is.na(nF)]
  names.sug.K <- names(nF)
  sug.K <- as.numeric(nF)
  names(sug.K) <- names.sug.K
  row.names(reference) <- paste0("resample.", row.names(reference))
  dat.eigen <- dat.eigen[order(row.names(dat.eigen)),]
  reference <- reference[order(row.names(reference)),]
  e.values <- rbind(dat.eigen, reference)

  P <- dataeigen <- NULL
  if(plot){
    df.e.values <- data.frame(dataeigen = rep(row.names(e.values), each = J), J = rep(1:J, nrow(e.values)), eigen = as.numeric(t(e.values)))
    point <- ifelse(grepl("dat", df.e.values$dataeigen), 16, 4)
    P <- ggplot2::ggplot(df.e.values, ggplot2::aes(x = J, y = eigen, color = dataeigen)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(shape = point) +
      ggplot2::scale_y_continuous("Eigenvalues") +
      ggplot2::scale_x_continuous("Component number") +
      ggplot2::theme_bw()
  }

  spec <- list(dat = dat, R = R, fa = fa, cor = cor, cutoff = cutoff, fm = fm, verbose = verbose, seed = seed)
  res <- list(sug.K = sug.K, e.values = e.values, plot = P, specifications = spec)
  class(res) <- "paK"
  return(res)
}
