cdmTools.fa <- function(r, nfactors = 1, n.obs = NA, n.iter = 1, rotate = "oblimin",
                        scores = "regression", residuals = FALSE, SMC = TRUE,
                        covar = FALSE, missing = FALSE, impute = "median",
                        min.err = 0.001, max.iter = 50, symmetric = TRUE, warnings = TRUE,
                        fm = "minres", alpha = 0.1, p = 0.05, oblique.scores = FALSE,
                        np.obs = NULL, use = "pairwise", cor = "cor",
                        correct = 0.5, weight = NULL, ...){
  cl <- match.call()
  if (psych::isCorrelation(r)) {
    if (is.na(n.obs) && (n.iter > 1))
      stop("You must specify the number of subjects if giving a correlation matrix and doing confidence intervals")
  }
  f <- cdmTools.fac(r = r, nfactors = nfactors, n.obs = n.obs, rotate = rotate,
                    scores = scores, residuals = residuals, SMC = SMC, covar = covar,
                    missing = missing, impute = impute, min.err = min.err,
                    max.iter = max.iter, symmetric = symmetric, warnings = warnings,
                    fm = fm, alpha = alpha, oblique.scores = oblique.scores,
                    np.obs = np.obs, use = use, cor = cor, correct = correct,
                    weight = weight, ... = ...)
  fl <- f$loadings
  nvar <- dim(fl)[1]
  if (n.iter > 1) {
    if (is.na(n.obs)) {
      n.obs <- f$n.obs
    }
    replicates <- list()
    rep.rots <- list()
    replicateslist <- parallel::mclapply(1:n.iter, function(x) {
      if (psych::isCorrelation(r)) {
        mu <- rep(0, nvar)
        eX <- eigen(r)
        X <- matrix(stats::rnorm(nvar * n.obs), n.obs)
        X <- t(eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*% t(X))
      }
      else {
        X <- r[sample(n.obs, n.obs, replace = TRUE),
        ]
      }
      fs <- cdmTools.fac(X, nfactors = nfactors, rotate = rotate,
                         scores = "none", SMC = SMC, missing = missing,
                         impute = impute, min.err = min.err, max.iter = max.iter,
                         symmetric = symmetric, warnings = warnings, fm = fm,
                         alpha = alpha, oblique.scores = oblique.scores,
                         np.obs = np.obs, use = use, cor = cor, correct = correct,
                         ... = ...)
      if (nfactors == 1) {
        replicates <- list(loadings = fs$loadings)
      }
      else {
        t.rot <- psych::target.rot(fs$loadings, fl)
        if (!is.null(fs$Phi)) {
          phis <- fs$Phi
          replicates <- list(loadings = t.rot$loadings,
                             phis = phis[lower.tri(t.rot$Phi)])
        }
        else {
          replicates <- list(loadings = t.rot$loadings)
        }
      }
    })
    replicates <- matrix(unlist(replicateslist), nrow = n.iter,
                         byrow = TRUE)
    means <- colMeans(replicates, na.rm = TRUE)
    sds <- apply(replicates, 2, stats::sd, na.rm = TRUE)
    if (length(means) > (nvar * nfactors)) {
      means.rot <- means[(nvar * nfactors + 1):length(means)]
      sds.rot <- sds[(nvar * nfactors + 1):length(means)]
      ci.rot.lower <- means.rot + stats::qnorm(p/2) * sds.rot
      ci.rot.upper <- means.rot + stats::qnorm(1 - p/2) * sds.rot
      ci.rot <- data.frame(lower = ci.rot.lower, upper = ci.rot.upper)
    }
    else {
      rep.rots <- NULL
      means.rot <- NULL
      sds.rot <- NULL
      z.rot <- NULL
      ci.rot <- NULL
    }
    means <- matrix(means[1:(nvar * nfactors)], ncol = nfactors)
    sds <- matrix(sds[1:(nvar * nfactors)], ncol = nfactors)
    tci <- abs(means)/sds
    ptci <- 1 - stats::pnorm(tci)
    if (!is.null(rep.rots)) {
      tcirot <- abs(means.rot)/sds.rot
      ptcirot <- 1 - stats::pnorm(tcirot)
    }
    else {
      tcirot <- NULL
      ptcirot <- NULL
    }
    ci.lower <- means + stats::qnorm(p/2) * sds
    ci.upper <- means + stats::qnorm(1 - p/2) * sds
    ci <- data.frame(lower = ci.lower, upper = ci.upper)
    class(means) <- "loadings"
    colnames(means) <- colnames(sds) <- colnames(fl)
    rownames(means) <- rownames(sds) <- rownames(fl)
    f$cis <- list(means = means, sds = sds, ci = ci, p = 2 *
                    ptci, means.rot = means.rot, sds.rot = sds.rot, ci.rot = ci.rot,
                  p.rot = ptcirot, Call = cl, replicates = replicates,
                  rep.rots = rep.rots)
    results <- f
    results$Call <- cl
    class(results) <- c("psych", "fa.ci")
  }
  else {
    results <- f
    results$Call <- cl
    class(results) <- c("psych", "fa")
  }
  return(results)
}
cdmTools.fac <- function(r, nfactors = 1, n.obs = NA, rotate = "oblimin",
                         scores = "tenBerge", residuals = FALSE, SMC = TRUE,
                         covar = FALSE, missing = FALSE, impute = "median",
                         min.err = 0.001, max.iter = 50, symmetric = TRUE, warnings = TRUE,
                         fm = "minres", alpha = 0.1, oblique.scores = FALSE,
                         np.obs = NULL, use = "pairwise", cor = "cor",
                         correct = 0.5, weight = NULL, ...){
  cl <- match.call()
  control <- NULL
  "fit.residuals" <- function(Psi, S, nf, S.inv = NULL,
                              fm) {
    diag(S) <- 1 - Psi
    if (!is.null(S.inv))
      sd.inv <- diag(1/diag(S.inv))
    eigens <- eigen(S)
    eigens$values[eigens$values < .Machine$double.eps] <- 100 *
      .Machine$double.eps
    if (nf > 1) {
      loadings <- eigens$vectors[, 1:nf] %*% diag(sqrt(eigens$values[1:nf]))
    }
    else {
      loadings <- eigens$vectors[, 1] * sqrt(eigens$values[1])
    }
    model <- loadings %*% t(loadings)
    switch(fm, wls = {
      residual <- sd.inv %*% (S - model)^2 %*% sd.inv
    }, gls = {
      residual <- (S.inv %*% (S - model))^2
    }, uls = {
      residual <- (S - model)^2
    }, ols = {
      residual <- (S - model)
      residual <- residual[lower.tri(residual)]
      residual <- residual^2
    }, minres = {
      residual <- (S - model)
      residual <- residual[lower.tri(residual)]
      residual <- residual^2
    }, old.min = {
      residual <- (S - model)
      residual <- residual[lower.tri(residual)]
      residual <- residual^2
    }, minchi = {
      residual <- (S - model)^2
      residual <- residual * np.obs
      diag(residual) <- 0
    })
    error <- sum(residual)
  }
  "fit" <- function(S, nf, fm, covar) {
    if (is.logical(SMC)) {
      S.smc <- psych::smc(S, covar)
    }
    else {
      S.smc <- SMC
    }
    upper <- max(S.smc, 1)
    if ((fm == "wls") | (fm == "gls")) {
      S.inv <- solve(S)
    }
    else {
      S.inv <- NULL
    }
    if (!covar && (sum(S.smc) == nf) && (nf > 1)) {
      start <- rep(0.5, nf)
    }
    else {
      start <- diag(S) - S.smc
    }
    if (fm == "ml" || fm == "mle") {
      res <- stats::optim(start, FAfn, FAgr, method = "L-BFGS-B",
                   lower = 0.005, upper = upper, control = c(list(fnscale = 1,
                                                                  parscale = rep(0.01, length(start))), control),
                   nf = nf, S = S)
    }
    else {
      if (fm == "ols") {
        if (is.logical(SMC)) {
          start <- diag(S) - psych::smc(S, covar)
        }
        else {
          start <- SMC
        }
        res <- stats::optim(start, FA.OLS, method = "L-BFGS-B",
                     lower = 0.005, upper = upper, control = c(list(fnscale = 1,
                                                                    parscale = rep(0.01, length(start)))), nf = nf,
                     S = S)
      }
      else {
        if ((fm == "minres") | (fm == "uls")) {
          start <- diag(S) - psych::smc(S, covar)
          res <- stats::optim(start, fit.residuals, gr = FAgr.minres,
                       method = "L-BFGS-B", lower = 0.005,
                       upper = upper, control = c(list(fnscale = 1,
                                                       parscale = rep(0.01, length(start)))),
                       nf = nf, S = S, fm = fm)
        }
        else {
          start <- psych::smc(S, covar)
          res <- stats::optim(start, fit.residuals, gr = FAgr.minres2,
                       method = "L-BFGS-B", lower = 0.005,
                       upper = upper, control = c(list(fnscale = 1,
                                                       parscale = rep(0.01, length(start)))),
                       nf = nf, S = S, S.inv = S.inv, fm = fm)
        }
      }
    }
    if ((fm == "wls") | (fm == "gls") | (fm ==
                                         "ols") | (fm == "uls") | (fm == "minres") |
        (fm == "old.min")) {
      Lambda <- FAout.wls(res$par, S, nf)
    }
    else {
      Lambda <- FAout(res$par, S, nf)
    }
    result <- list(loadings = Lambda, res = res, S = S)
  }
  FAfn <- function(Psi, S, nf) {
    sc <- diag(1/sqrt(Psi))
    Sstar <- sc %*% S %*% sc
    E <- eigen(Sstar, symmetric = TRUE, only.values = TRUE)
    e <- E$values[-(1:nf)]
    e <- sum(log(e) - e) - nf + nrow(S)
    -e
  }
  FAgr <- function(Psi, S, nf) {
    sc <- diag(1/sqrt(Psi))
    Sstar <- sc %*% S %*% sc
    E <- eigen(Sstar, symmetric = TRUE)
    L <- E$vectors[, 1:nf, drop = FALSE]
    load <- L %*% diag(sqrt(pmax(E$values[1:nf] - 1, 0)),
                       nf)
    load <- diag(sqrt(Psi)) %*% load
    g <- load %*% t(load) + diag(Psi) - S
    diag(g)/Psi^2
  }
  FAgr.minres2 <- function(Psi, S, nf, S.inv, fm) {
    sc <- diag(1/sqrt(Psi))
    Sstar <- sc %*% S %*% sc
    E <- eigen(Sstar, symmetric = TRUE)
    L <- E$vectors[, 1:nf, drop = FALSE]
    load <- L %*% diag(sqrt(pmax(E$values[1:nf] - 1, 0)),
                       nf)
    load <- diag(sqrt(Psi)) %*% load
    g <- load %*% t(load) + diag(Psi) - S
    if (fm == "minchi") {
      g <- g * np.obs
    }
    diag(g)/Psi^2
  }
  FAgr.minres <- function(Psi, S, nf, fm) {
    Sstar <- S - diag(Psi)
    E <- eigen(Sstar, symmetric = TRUE)
    L <- E$vectors[, 1:nf, drop = FALSE]
    load <- L %*% diag(sqrt(pmax(E$values[1:nf], 0)), nf)
    g <- load %*% t(load) + diag(Psi) - S
    diag(g)
  }
  FAout <- function(Psi, S, q) {
    sc <- diag(1/sqrt(Psi))
    Sstar <- sc %*% S %*% sc
    E <- eigen(Sstar, symmetric = TRUE)
    L <- E$vectors[, 1L:q, drop = FALSE]
    load <- L %*% diag(sqrt(pmax(E$values[1L:q] - 1, 0)),
                       q)
    diag(sqrt(Psi)) %*% load
  }
  FAout.wls <- function(Psi, S, q) {
    diag(S) <- diag(S) - Psi
    E <- eigen(S, symmetric = TRUE)
    L <- E$vectors[, 1L:q, drop = FALSE] %*% diag(sqrt(pmax(E$values[1L:q,
                                                                     drop = FALSE], 0)), q)
    return(L)
  }
  "MRFA" <- function(S, nf) {
    com.glb <- psych::glb.algebraic(S)
    L <- FAout.wls(1 - com.glb$solution, S, nf)
    h2 <- com.glb$solution
    result <- list(loadings = L, communality = h2)
  }
  FA.OLS <- function(Psi, S, nf) {
    E <- eigen(S - diag(Psi), symmetric = T)
    U <- E$vectors[, 1:nf, drop = FALSE]
    D <- E$values[1:nf, drop = FALSE]
    D[D < 0] <- 0
    if (length(D) < 2) {
      L <- U * sqrt(D)
    }
    else {
      L <- U %*% diag(sqrt(D))
    }
    model <- L %*% t(L)
    diag(model) <- diag(S)
    return(sum((S - model)^2)/2)
  }
  FAgr.OLS <- function(Psi, S, nf) {
    E <- eigen(S - diag(Psi), symmetric = TRUE)
    U <- E$vectors[, 1:nf, drop = FALSE]
    D <- E$values[1:nf]
    D[D < 0] <- 0
    L <- U %*% diag(sqrt(D))
    model <- L %*% t(L)
    g <- diag(Psi) - diag(S - model)
    diag(g)/Psi^2
  }
  if (fm == "mle" || fm == "MLE" || fm == "ML")
    fm <- "ml"
  if (!any(fm %in% (c("pa", "alpha", "minrank",
                      "wls", "gls", "minres", "minchi",
                      "uls", "ml", "mle", "ols", "old.min")))) {
    message("factor method not specified correctly, minimum residual (unweighted least squares  used")
    fm <- "minres"
  }
  x.matrix <- r
  n <- dim(r)[2]
  if (!psych::isCorrelation(r) & !psych::isCovariance(r)) {
    matrix.input <- FALSE
    n.obs <- dim(r)[1]
    if (missing) {
      x.matrix <- as.matrix(x.matrix)
      miss <- which(is.na(x.matrix), arr.ind = TRUE)
      if (impute == "mean") {
        item.means <- colMeans(x.matrix, na.rm = TRUE)
        x.matrix[miss] <- item.means[miss[, 2]]
      }
      else {
        item.med <- apply(x.matrix, 2, stats::median, na.rm = TRUE)
        x.matrix[miss] <- item.med[miss[, 2]]
      }
    }
    np.obs <- psych::pairwiseCount(r)
    if (covar) {
      cor <- "cov"
    }
    switch(cor, cor = {
      if (!is.null(weight)) {
        r <- psych::cor.wt(r, w = weight)$r
      } else {
        r <- stats::cor(r, use = use)
      }
    }, cov = {
      r <- stats::cov(r, use = use)
      covar <- TRUE
    }, wtd = {
      r <- psych::cor.wt(r, w = weight)$r
    }, spearman = {
      r <- cor(r, use = use, method = "spearman")
    }, kendall = {
      r <- cor(r, use = use, method = "kendall")
    }, tet = {
      r <- sirt::tetrachoric2(r)$rho
    }, poly = {
      r <- sirt::polychoric2(r)$rho
    }, tetrachoric = {
      r <- sirt::tetrachoric2(r)$rho
    }, polychoric = {
      r <- sirt::polychoric2(r)$rho
    }, mixed = {
      r <- psych::mixedCor(r, use = use, correct = correct)$rho
    }, Yuleb = {
      r <- psych::YuleCor(r, , bonett = TRUE)$rho
    }, YuleQ = {
      r <- psych::YuleCor(r, 1)$rho
    }, YuleY = {
      r <- psych::YuleCor(r, 0.5)$rho
    })
  }
  else {
    matrix.input <- TRUE
    if (fm == "minchi") {
      if (is.null(np.obs)) {
        fm <- "minres"
        message("factor method minchi does not make sense unless we know the sample size, minres used instead")
      }
    }
    if (is.na(n.obs) && !is.null(np.obs))
      n.obs <- max(as.vector(np.obs))
    if (!is.matrix(r)) {
      r <- as.matrix(r)
    }
    if (!covar) {
      r <- stats::cov2cor(r)
    }
  }
  if (!residuals) {
    result <- list(values = c(rep(0, n)), rotation = rotate,
                   n.obs = n.obs, np.obs = np.obs, communality = c(rep(0,
                                                                       n)), loadings = matrix(rep(0, n * n), ncol = n),
                   fit = 0)
  }
  else {
    result <- list(values = c(rep(0, n)), rotation = rotate,
                   n.obs = n.obs, np.obs = np.obs, communality = c(rep(0,
                                                                       n)), loadings = matrix(rep(0, n * n), ncol = n),
                   residual = matrix(rep(0, n * n), ncol = n), fit = 0,
                   r = r)
  }
  if (is.null(SMC))
    SMC = TRUE
  r.mat <- r
  Phi <- NULL
  colnames(r.mat) <- rownames(r.mat) <- colnames(r)
  if (any(is.na(r))) {
    bad <- TRUE
    tempr <- r
    wcl <- NULL
    while (bad) {
      wc <- table(which(is.na(tempr), arr.ind = TRUE))
      wcl <- c(wcl, as.numeric(names(which(wc == max(wc)))))
      tempr <- r[-wcl, -wcl]
      if (any(is.na(tempr))) {
        bad <- TRUE
      }
      else {
        bad <- FALSE
      }
    }
    cat("\nLikely variables with missing values are ",
        colnames(r)[wcl], " \n")
    stop("I am sorry: missing values (NAs) in the correlation matrix do not allow me to continue.\nPlease drop those variables and try again.")
  }
  if (is.logical(SMC)) {
    if (SMC) {
      if (nfactors <= n) {
        diag(r.mat) <- psych::smc(r, covar = covar)
      }
      else {
        if (warnings) {
          message("In fa, too many factors requested for this number of variables to use SMC for communality estimates, 1s are used instead")
        }
      }
    }
    else {
      diag(r.mat) <- 1
    }
  }
  else {
    diag(r.mat) <- SMC
  }
  orig <- diag(r)
  comm <- sum(diag(r.mat))
  err <- comm
  i <- 1
  comm.list <- list()
  if (fm == "alpha") {
    i <- 1
    e.values <- eigen(r, symmetric = symmetric)$values
    H2 <- diag(r.mat)
    while (err > min.err) {
      r.mat <- stats::cov2cor(r.mat)
      eigens <- eigen(r.mat, symmetric = symmetric)
      loadings <- eigens$vectors[, 1:nfactors, drop = FALSE] %*%
        diag(sqrt(eigens$values[1:nfactors, drop = FALSE]))
      model <- loadings %*% t(loadings)
      newH2 <- H2 * diag(model)
      err <- sum(abs(H2 - newH2))
      r.mat <- r
      diag(r.mat) <- newH2
      H2 <- newH2
      i <- i + 1
      if (i > max.iter) {
        if (warnings) {
          message("maximum iteration exceeded")
        }
        err <- 0
      }
    }
    loadings <- sqrt(H2) * loadings
    eigens <- sqrt(H2) * eigens$vaues
    comm1 <- sum(H2)
  }
  if (fm == "pa") {
    e.values <- eigen(r, symmetric = symmetric)$values
    while (err > min.err) {
      eigens <- eigen(r.mat, symmetric = symmetric)
      if (nfactors > 1) {
        loadings <- eigens$vectors[, 1:nfactors] %*%
          diag(sqrt(eigens$values[1:nfactors]))
      }
      else {
        loadings <- eigens$vectors[, 1] * sqrt(eigens$values[1])
      }
      model <- loadings %*% t(loadings)
      new <- diag(model)
      comm1 <- sum(new)
      diag(r.mat) <- new
      err <- abs(comm - comm1)
      if (is.na(err)) {
        warning("imaginary eigen value condition encountered in fa\n Try again with SMC=FALSE \n exiting fa")
        break
      }
      comm <- comm1
      comm.list[[i]] <- comm1
      i <- i + 1
      if (i > max.iter) {
        if (warnings) {
          message("maximum iteration exceeded")
        }
        err <- 0
      }
    }
    eigens <- eigens$values
  }
  if (fm == "minrank") {
    mrfa <- MRFA(r, nfactors)
    loadings <- mrfa$loadings
    model <- loadings %*% t(loadings)
    e.values <- eigen(r)$values
    S <- r
    diag(S) <- diag(model)
    eigens <- eigen(S)$values
  }
  if ((fm == "wls") | (fm == "minres") | (fm ==
                                          "minchi") | (fm == "gls") | (fm == "uls") |
      (fm == "ml") | (fm == "mle") | (fm == "ols") |
      (fm == "old.min")) {
    uls <- fit(r, nfactors, fm, covar = covar)
    e.values <- eigen(r)$values
    result.res <- uls$res
    loadings <- uls$loadings
    model <- loadings %*% t(loadings)
    S <- r
    diag(S) <- diag(model)
    eigens <- eigen(S)$values
  }
  if (!is.double(loadings)) {
    warning("the matrix has produced imaginary results -- proceed with caution")
    loadings <- matrix(as.double(loadings), ncol = nfactors)
  }
  if (nfactors > 1) {
    sign.tot <- vector(mode = "numeric", length = nfactors)
    sign.tot <- sign(colSums(loadings))
    sign.tot[sign.tot == 0] <- 1
    loadings <- loadings %*% diag(sign.tot)
  }
  else {
    if (sum(loadings) < 0) {
      loadings <- -as.matrix(loadings)
    }
    else {
      loadings <- as.matrix(loadings)
    }
    colnames(loadings) <- "MR1"
  }
  switch(fm, alpha = {
    colnames(loadings) <- paste("alpha", 1:nfactors,
                                sep = "")
  }, wls = {
    colnames(loadings) <- paste("WLS", 1:nfactors,
                                sep = "")
  }, pa = {
    colnames(loadings) <- paste("PA", 1:nfactors, sep = "")
  }, gls = {
    colnames(loadings) <- paste("GLS", 1:nfactors,
                                sep = "")
  }, ml = {
    colnames(loadings) <- paste("ML", 1:nfactors, sep = "")
  }, minres = {
    colnames(loadings) <- paste("MR", 1:nfactors, sep = "")
  }, minrank = {
    colnames(loadings) <- paste("MRFA", 1:nfactors,
                                sep = "")
  }, minchi = {
    colnames(loadings) <- paste("MC", 1:nfactors, sep = "")
  })
  rownames(loadings) <- rownames(r)
  loadings[loadings == 0] <- 10^-15
  model <- loadings %*% t(loadings)
  f.loadings <- loadings
  rot.mat <- NULL
  if (rotate != "none") {
    if (nfactors > 1) {
      if (rotate == "varimax" | rotate == "Varimax" |
          rotate == "quartimax" | rotate == "bentlerT" |
          rotate == "geominT" | rotate == "targetT" |
          rotate == "bifactor" | rotate == "TargetT" |
          rotate == "equamax" | rotate == "varimin" |
          rotate == "Promax" |
          rotate == "promax" | rotate == "cluster" |
          rotate == "biquartimin" |
          rotate == "TargetQ") {
        Phi <- NULL
        switch(rotate, varimax = {
          rotated <- stats::varimax(loadings)
          loadings <- rotated$loadings
          rot.mat <- rotated$rotmat
        }, Varimax = {
          if (!requireNamespace("GPArotation")) {
            stop("I am sorry, to do this rotation requires the GPArotation package to be installed")
          }
          rotated <- GPArotation::Varimax(loadings, ...)
          loadings <- rotated$loadings
          rot.mat <- t(solve(rotated$Th))
        }, quartimax = {
          if (!requireNamespace("GPArotation")) {
            stop("I am sorry, to do this rotation requires the GPArotation package to be installed")
          }
          rotated <- GPArotation::quartimax(loadings,
                                            ...)
          loadings <- rotated$loadings
          rot.mat <- t(solve(rotated$Th))
        }, bentlerT = {
          if (!requireNamespace("GPArotation")) {
            stop("I am sorry, to do this rotation requires the GPArotation package to be installed")
          }
          rotated <- GPArotation::bentlerT(loadings,
                                           ...)
          loadings <- rotated$loadings
          rot.mat <- t(solve(rotated$Th))
        }, geominT = {
          if (!requireNamespace("GPArotation")) {
            stop("I am sorry, to do this rotation requires the GPArotation package to be installed")
          }
          rotated <- GPArotation::geominT(loadings, ...)
          loadings <- rotated$loadings
          rot.mat <- t(solve(rotated$Th))
        }, targetT = {
          if (!requireNamespace("GPArotation")) {
            stop("I am sorry, to do this rotation requires the GPArotation package to be installed")
          }
          rotated <- GPArotation::targetT(loadings, Tmat = diag(ncol(loadings)),
                                          ...)
          loadings <- rotated$loadings
          rot.mat <- t(solve(rotated$Th))
        }, bifactor = {
          rot <- psych::bifactor(loadings, ...)
          loadings <- rot$loadings
          rot.mat <- t(solve(rot$Th))
        }, TargetT = {
          if (!requireNamespace("GPArotation")) {
            stop("I am sorry, to do this rotation requires the GPArotation package to be installed")
          }
          rot <- GPArotation::targetT(loadings, Tmat = diag(ncol(loadings)),
                                      ...)
          loadings <- rot$loadings
          rot.mat <- t(solve(rot$Th))
        }, equamax = {
          rot <- psych::equamax(loadings, ...)
          loadings <- rot$loadings
          rot.mat <- t(solve(rot$Th))
        }, varimin = {
          rot <- psych::varimin(loadings, ...)
          loadings <- rot$loadings
          rot.mat <- t(solve(rot$Th))
        }, Promax = {
          pro <- psych::Promax(loadings, ...)
          loadings <- pro$loadings
          Phi <- pro$Phi
          rot.mat <- pro$rotmat
        }, promax = {
          pro <- psych::kaiser(loadings, rotate = "Promax",
                               ...)
          loadings <- pro$loadings
          rot.mat <- pro$rotmat
          Phi <- pro$Phi
        }, cluster = {
          loadings <- stats::varimax(loadings, ...)$loadings
          pro <- psych::target.rot(loadings)
          loadings <- pro$loadings
          Phi <- pro$Phi
          rot.mat <- pro$rotmat
        }, biquartimin = {
          ob <- psych::biquartimin(loadings, ...)
          loadings <- ob$loadings
          Phi <- ob$Phi
          rot.mat <- t(solve(ob$Th))
        }, TargetQ = {
          ob <- psych::TargetQ(loadings, ...)
          loadings <- ob$loadings
          Phi <- ob$Phi
          rot.mat <- t(solve(ob$Th))
        })
      }
      else {
        if (rotate == "oblimin" | rotate == "quartimin" |
            rotate == "simplimax" | rotate == "geominQ" |
            rotate == "bentlerQ" | rotate == "targetQ") {
          if (!requireNamespace("GPArotation")) {
            warning("I am sorry, to do these rotations requires the GPArotation package to be installed")
            Phi <- NULL
          }
          else {
            ob <- try(do.call(utils::getFromNamespace(rotate,
                                               "GPArotation"), list(loadings, ...)))
            if (inherits(ob, as.character("try-error"))) {
              warning("The requested transformaton failed, Promax was used instead as an oblique transformation")
              ob <- psych::Promax(loadings)
            }
            loadings <- ob$loadings
            Phi <- ob$Phi
            rot.mat <- t(solve(ob$Th))
          }
        }
        else {
          message("Specified rotation not found, rotate='none' used")
        }
      }
    }
  }
  signed <- sign(colSums(loadings))
  signed[signed == 0] <- 1
  loadings <- loadings %*% diag(signed)
  if (!is.null(Phi)) {
    Phi <- diag(signed) %*% Phi %*% diag(signed)
  }
  switch(fm, alpha = {
    colnames(loadings) <- paste("alpha", 1:nfactors,
                                sep = "")
  }, wls = {
    colnames(loadings) <- paste("WLS", 1:nfactors,
                                sep = "")
  }, pa = {
    colnames(loadings) <- paste("PA", 1:nfactors, sep = "")
  }, gls = {
    colnames(loadings) <- paste("GLS", 1:nfactors,
                                sep = "")
  }, ml = {
    colnames(loadings) <- paste("ML", 1:nfactors, sep = "")
  }, minres = {
    colnames(loadings) <- paste("MR", 1:nfactors, sep = "")
  }, minrank = {
    colnames(loadings) <- paste("MRFA", 1:nfactors,
                                sep = "")
  }, uls = {
    colnames(loadings) <- paste("ULS", 1:nfactors,
                                sep = "")
  }, old.min = {
    colnames(loadings) <- paste0("oldmin", 1:nfactors)
  }, minchi = {
    colnames(loadings) <- paste("MC", 1:nfactors, sep = "")
  })
  if (nfactors > 1) {
    ev.rotated <- diag(t(loadings) %*% loadings)
    ev.order <- order(ev.rotated, decreasing = TRUE)
    loadings <- loadings[, ev.order]
  }
  rownames(loadings) <- colnames(r)
  if (!is.null(Phi)) {
    Phi <- Phi[ev.order, ev.order]
  }
  class(loadings) <- "loadings"
  if (nfactors < 1)
    nfactors <- n
  result <- psych::factor.stats(r, loadings, Phi, n.obs = n.obs, np.obs = np.obs,
                                alpha = alpha)
  result$rotation <- rotate
  result$communality <- diag(model)
  if (max(result$communality > 1) && !covar)
    warning("An ultra-Heywood case was detected.  Examine the results carefully")
  if (fm == "minrank") {
    result$communalities <- mrfa$communality
  }
  else {
    if (fm == "pa" | fm == "alpha") {
      result$communalities <- comm1
    }
    else {
      result$communalities <- 1 - result.res$par
    }
  }
  result$uniquenesses <- diag(r - model)
  result$values <- eigens
  result$e.values <- e.values
  result$loadings <- loadings
  result$model <- model
  result$fm <- fm
  result$rot.mat <- rot.mat
  if (!is.null(Phi)) {
    colnames(Phi) <- rownames(Phi) <- colnames(loadings)
    result$Phi <- Phi
    Structure <- loadings %*% Phi
  }
  else {
    Structure <- loadings
  }
  class(Structure) <- "loadings"
  result$Structure <- Structure
  if (fm == "pa")
    result$communality.iterations <- unlist(comm.list)
  result$method = scores
  if (oblique.scores) {
    result$scores <- psych::factor.scores(x.matrix, f = loadings,
                                          Phi = NULL, method = scores)
  }
  else {
    result$scores <- psych::factor.scores(x.matrix, f = loadings,
                                          Phi = Phi, method = scores)
  }
  if (is.null(result$scores$R2))
    result$scores$R2 <- NA
  result$R2.scores <- result$scores$R2
  result$weights <- result$scores$weights
  result$scores <- result$scores$scores
  if (!is.null(result$scores))
    colnames(result$scores) <- colnames(loadings)
  result$factors <- nfactors
  result$r <- r
  result$np.obs <- np.obs
  result$fn <- "fa"
  result$fm <- fm
  if (is.null(Phi)) {
    if (nfactors > 1) {
      vx <- colSums(loadings^2)
    }
    else {
      vx <- sum(loadings^2)
    }
  }
  else {
    vx <- diag(Phi %*% t(loadings) %*% loadings)
  }
  vtotal <- sum(result$communality + result$uniquenesses)
  names(vx) <- colnames(loadings)
  varex <- rbind(`SS loadings` = vx)
  varex <- rbind(varex, `Proportion Var` = vx/vtotal)
  if (nfactors > 1) {
    varex <- rbind(varex, `Cumulative Var` = cumsum(vx/vtotal))
    varex <- rbind(varex, `Proportion Explained` = vx/sum(vx))
    varex <- rbind(varex, `Cumulative Proportion` = cumsum(vx/sum(vx)))
  }
  result$Vaccounted <- varex
  result$Call <- cl
  class(result) <- c("psych", "fa")
  return(result)
}
lctolg.helper <- function(K){
  LC <- GDINA::attributepattern(K)
  lc <- apply(LC, 1, paste, collapse = "")
  # lc2 <- as.vector(sapply(lc, function(x) paste0(which(as.numeric(unlist(strsplit(x, ""))) == 1), collapse = "")))
  lKq <- as.vector(unlist(sapply(1:(K-1), function(x) rep(2^x, factorial(K) / (factorial(x) * factorial(K - x))))))
  q1 <- unlist(sapply(1:length(lKq), function(x) rep(lc[-c(1, length(lc))][x], lKq[x])))
  q2 <- as.vector(sapply(q1, function(x) paste0(which(as.numeric(unlist(strsplit(x, ""))) == 1), collapse = "")))
  q3 <- unlist(sapply(1:length(lKq), function(x) apply(GDINA::attributepattern(log(lKq, 2)[x]), 1, paste, collapse = "")))
  dm <- matrix(NA, nrow = length(lc), ncol = length(q1), dimnames = list(lc, q1))
  for(C in 1:ncol(dm)){
    pos <- as.numeric(unlist(strsplit(q2[C], "")))
    key <- as.numeric(unlist(strsplit(q3[C], "")))
    if(length(pos) == 1){
      dm[LC[,pos] == key, C] <- 1
      dm[LC[,pos] != key, C] <- 0
    } else {
      dm[apply(apply(LC[,pos], 1, function(x) x == key), 2, all), C] <- 1
      dm[!apply(apply(LC[,pos], 1, function(x) x == key), 2, all), C] <- 0
    }
  }
  return(dm)
}
extract.Hull <- function(fit, what = "PVAF", R2.aux = NULL){
  if(!all(what %in% c("PVAF", "R2"))){stop("what = 'PVAF', 'R2'")}
  res <- list()
  if("PVAF" %in% what){res[["PVAF"]] <- GDINA::Qval(fit, digits = 8)$PVAF; colnames(res[["PVAF"]]) <- paste0("J", 1:ncol(res[["PVAF"]]))}
  if("R2" %in% what){
    dat <- fit$options$dat
    Q <- fit$options$Q
    N <- nrow(dat)
    K <- ncol(Q)
    J <- nrow(Q)
    qM <- GDINA::attributepattern(K)[-1,]
    qS <- apply(qM, 1, paste, collapse = "")
    np <- matrix(2^(rowSums(qM)), byrow = T, nrow = J, ncol = nrow(qM))

    LL.0j <- apply(dat, 2, function(x) sum((x * log(mean(x))) + ((1 - x) * log(1 - mean(x)))))
    post <- exp(GDINA::indlogPost(fit)) # Posterior probailities for examinees in each latent class
    lc.n <- colSums(post) # Expected number of examinees in each latent class
    lc.r <- t(sapply(1:J, function(x) colSums(post * dat[,x]))) # Expected number of correct responses for each item and latent class
    rownames(lc.r) <- paste0("J", 1:J)
    lc.pc <- apply(lc.r, 1, function(x) x / lc.n) # Expected probability of correct responses for each item and latent class
    if(is.null(R2.aux)){
      DM <- lctolg.helper(K)
    } else {
      DM <- R2.aux
    }
    lg.n <- as.vector(lc.n %*% DM)
    names(lg.n) <- colnames(DM)
    lg.r <- lc.r %*% DM
    lg.pc <- apply(lg.r, 1, function(x) x / lg.n)
    lg.post <- post %*% DM

    j.lg.index <- c(0, cumsum(as.vector(unlist(sapply(1:(K-1), function(x) rep(2^x, factorial(K) / (factorial(x) * factorial(K - x))))))))
    j.lg.index <- sapply(1:(length(j.lg.index) - 1), function(x) (j.lg.index[x] + 1):j.lg.index[x + 1])

    R2 <- matrix(NA, nrow = J, ncol = 2^K - 1, dimnames = list(1:J, qS))
    exp.cor <- array(NA, dim = c(N, 2^K - 1, J))

    for(j in 1:J){
      for(q in 1:(2^K - 2)){
        exp.cor[, q, j] <- lg.post[,j.lg.index[[q]], drop = FALSE] %*% lg.pc[j.lg.index[[q]], j, drop = FALSE]
        exp.cor[, q, j][round(exp.cor[, q, j], 14) == 1] <- 1 - 1e-15 # Avoid round 1
        LL.tmp <- sum((dat[,j] * log(exp.cor[, q, j])) + ((1 - dat[,j]) * log(1 - exp.cor[, q, j])))
        R2.McF <- 1 - (LL.tmp / LL.0j[j])
        R2[j, q] <- R2.McF
      }
      q <- q + 1
      exp.cor[, q, j] <- post %*% lc.pc[, j]
      exp.cor[, q, j][round(exp.cor[, q, j], 14) == 1] <- 1 - 1e-15 # Avoid round 1
      LL.tmp <- sum((dat[,j] * log(exp.cor[, q, j])) + ((1 - dat[,j]) * log(1 - exp.cor[, q, j])))
      R2.McF <- 1 - (LL.tmp / LL.0j[j])
      R2[j, q] <- R2.McF
    }

    res[["R2"]] <- t(R2)
  }
  return(res)
}
best.Kj <- function(M, best.pos = NULL){
  K <- log(nrow(M) + 1, 2)
  J <- ncol(M)
  Kj <- rowSums(GDINA::attributepattern(K)[-1,])
  if(is.null(best.pos)){
    maxloc <- matrix(NA, nrow = K, ncol = J, dimnames = list(paste0("K", 1:K), paste0("J", 1:J)))
    for(kj in 1:K){maxloc[kj,] <- apply(M[Kj == kj,, drop = FALSE], 2, which.max)}
    best.pos <- cumsum(c(0, table(Kj)[-K])) + maxloc
    best.M <- sapply(1:J, function(j) M[best.pos[,j],j])
    rownames(best.M) <- paste0("K", 1:K)
    colnames(best.M) <- paste0("J", 1:J)
  } else {
    best.M <- sapply(1:J, function(j) M[best.pos[,j],j])
    rownames(best.M) <- paste0("K", 1:K)
    colnames(best.M) <- paste0("J", 1:J)
  }
  return(list(best.pos = best.pos, best.M = best.M))
}
hull <- function(best.M, np, trim = T){
  if(is.vector(best.M)){best.M <- matrix(best.M, ncol = 1)}
  K <- nrow(best.M)
  J <- ncol(best.M)
  res <- best.M/NA
  for(j in 1:J){
    gf.j <- c(0, best.M[,j])
    np.j <- np
    if(trim){
      i <- 0
      while(1 < 2){
        i <- i + 1
        if((i + 2) > length(np.j)){break}
        np.dif1 <- np.j[i + 1] - np.j[i]; gf.dif1 <- gf.j[i + 1] - gf.j[i]
        np.dif2 <- np.j[i + 2] - np.j[i]; gf.dif2 <- gf.j[i + 2] - gf.j[i]
        at1 <- atan(gf.dif1/np.dif1)
        at2 <- atan(gf.dif2/np.dif2)
        if(at2 > at1){
          np.j <- np.j[-(i + 1)]
          gf.j <- gf.j[-(i + 1)]
          i <- 0
        } else {
          next
        }
      }
    }

    if(length(gf.j) == 2){
      res[,j] <- c(rep(NA, K - 1), 1)
    } else {
      st <- c()
      for(i in 2:(length(np.j) - 1)){
        st <- c(st, (((gf.j[i] - gf.j[i - 1]) / (np.j[i] - np.j[i - 1])) / ((gf.j[i + 1] - gf.j[i]) / (np.j[i + 1] - np.j[i]))))
      }
      res[names(st),j] <- st
    }
  }
  return(res)
}
select.q <- function(M, best.pos = NULL){
  bKj <- best.Kj(M, best.pos)
  K <- nrow(bKj$best.pos)
  J <- ncol(bKj$best.pos)
  H <- NULL
  np <- c(0, 2^(1:K))
  H.tmp <- H <- hull(bKj$best.M, np, TRUE)
  H.tmp[is.na(H.tmp)] <- 0
  bK <- sapply(1:J, function(j) which.max(stats::na.omit(abs(H.tmp[,j])))); names(bK) <- paste0("J", 1:J)
  bq <- sapply(1:J, function(j) bKj$best.pos[bK[j], j])
  sug.Q <- GDINA::attributepattern(K)[1 + bq,]
  return(list(sug.Q = sug.Q, bKj = bKj, H = H, bK = bK, bq = bq))
}
bootSE.parallel <- function(fit, bootsample = 50, type = "nonparametric", n.cores = 1, verbose = TRUE, seed = 12345){
  if(exists(".Random.seed", .GlobalEnv)){
    oldseed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- oldseed)
  } else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv))
  }
  set.seed(seed)
  Y <- GDINA::extract(fit, "dat")
  Q <- GDINA::extract(fit, "Q")
  N <- nrow(Y)
  J <- ncol(Y)
  K <- ncol(Q)
  no.mg <- GDINA::extract(fit, "ngroup")
  stopifnot(no.mg == 1)
  lambda <- delta <- itemprob <- vector("list", bootsample)
  GDINA.options <- formals(GDINA::GDINA)
  GDINA.options <- GDINA.options[-c(1, 2, length(GDINA.options))]
  tmp <- as.list(fit$extra$call)[-c(1:3)]
  GDINA.options[names(GDINA.options) %in% names(tmp)] <- tmp
  GDINA.options$verbose <- 0
  GDINA.options$model <- fit$model
  att <- GDINA::extract(fit, "attributepattern")

  cl <- parallel::makeCluster(n.cores, type = "SOCK")
  doSNOW::registerDoSNOW(cl)

  comb <- function(x, ...){lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))}

  if(verbose){
    cat("Bootstrapping Progress:", "\n")
    pb <- utils::txtProgressBar(max = bootsample, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {
    opts <- NULL
  }

  boot.parallel = foreach::foreach(r = 1:bootsample,
                                   .options.snow = opts,
                                   .combine = 'comb', .multicombine = TRUE, .packages = "GDINA",
                                   .init = list(lambda = list(), itemprob = list(), delta = list())) %dopar% {
                                     set.seed(r * seed)
                                     if(tolower(type) == "parametric"){
                                       simdat <- GDINA::simGDINA(N, Q, catprob.parm = fit$catprob.parm, attribute = att[sample(seq_len(nrow(att)), N, replace = TRUE, prob = fit$posterior.prob), ])$dat
                                     } else if (tolower(type) == "nonparametric"){
                                       constant <- T
                                       while(constant){
                                         simdat <- Y[sample(1:N, N, replace = TRUE),]
                                         constant <- any(apply(simdat, 2, stats::sd, na.rm = T) == 0)
                                       }
                                     }
                                     boot.out <- do.call(GDINA::GDINA, c(list(dat = simdat, Q = Q), GDINA.options))
                                     lambda <- boot.out$struc.parm
                                     itemprob <- boot.out$catprob.parm
                                     delta <- boot.out$delta.parm

                                     return(list(lambda = lambda, itemprob = itemprob, delta = delta))
                                   }
  parallel::stopCluster(cl)
  names(boot.parallel) <- c("lambda", "itemprob", "delta")

  se.ip <- lapply(do.call(Map, c(f = "rbind", boot.parallel$itemprob)), function(x) apply(x, 2, stats::sd))
  se.d <- lapply(do.call(Map, c(f = "rbind", boot.parallel$delta)), function(x) apply(x, 2, stats::sd))
  se.lambda <- apply(do.call(rbind, boot.parallel$lambda), 2, stats::sd)
  if(GDINA::extract(fit, "att.dist") == "higher.order"){
    se.jointAtt <- matrix(se.lambda, ncol = 2)
  } else {
    se.jointAtt <- se.lambda
  }

  res <- list(itemparm.se = se.ip, delta.se = se.d, lambda.se = se.jointAtt, boot.est = list(lambda = boot.parallel$lambda, itemprob = boot.parallel$itemprob, delta = boot.parallel$delta))
  return(res)
}
phi.ML <- function(phi, dist, J, posterior = FALSE, att.prior = NULL){
  if(is.null(att.prior)){att.prior <- rep(1 / nrow(dist), nrow(dist))}
  L.il <- t(phi^dist * (1 - phi)^(J - dist))
  mL.il <- t(t(L.il) * att.prior)
  if(posterior){
    pp.il <- t(apply(mL.il, 1, function(i) i / sum(i)))
    lp <- colMeans(pp.il)
    return(sum(log(apply(L.il %*% lp, 1, prod))))
  } else {
    return(sum(log(apply(mL.il, 1, sum))))
  }
}
cdmTools.AlphaNP <- function(Y, Q, gate = c("AND", "OR"), method = c("Hamming", "Weighted", "Penalized"), wg = 1, ws = 1){
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  gate <- match.arg(gate)
  method <- match.arg(method)
  nperson <- dim(Y)[1]
  nitem <- dim(Q)[1]
  natt <- dim(Q)[2]
  M <- 2^natt
  pattern <- cdmTools.AlphaPermute(natt)
  Ideal <- matrix(NA, M, nitem)
  for(m in 1:M){
    for(j in 1:nitem){
      if(gate == "AND"){
        u <- prod(pattern[m, ]^Q[j, ])
      } else if(gate == "OR"){
        u <- 1 - prod((1 - pattern[m, ])^Q[j, ])
      } else {
        return(warning("Gate specification not valid."))
      }
      Ideal[m, j] <- u
    }
  }
  if(method == "Hamming"){
    weight <- rep(1, nitem)
    ws <- wg <- 1
  } else if(method == "Weighted"){
    p.bar <- apply(Y, 2, function(x) mean(x, na.rm = TRUE))
    weight <- 1/(p.bar * (1 - p.bar))
    weight[weight > 1/(0.95 * 0.05)] <- 1/(0.95 * 0.05)
    ws <- wg <- 1
  } else if (method == "Penalized"){
    p.bar <- apply(Y, 2, function(x) mean(x, na.rm = TRUE))
    weight <- 1/(p.bar * (1 - p.bar))
    weight[weight > 1/(0.95 * 0.05)] <- 1/(0.95 * 0.05)
    if(ws == wg){warning("Penalzing weights for guess and slip are the same --> equivalent with the \"Weighted\" method.")}
  } else {
    return(warning("Method specification not valid."))
  }
  loss.matrix <- matrix(NA, nrow = M, ncol = nperson)
  est.class <- NULL
  est.pattern <- NULL
  n.tie <- rep(0, nperson)
  for(i in 1:nperson){
    valid <- !is.na(Y[i, ])
    Y.matrix <- matrix(rep(Y[i, ], M), M, nitem, byrow = TRUE)
    loss <- apply(matrix(rep(weight[valid], M), M, sum(valid), byrow = TRUE) * (wg * abs(Y.matrix[,valid] - Ideal[,valid]) * Y.matrix[,valid] + ws * abs(Y.matrix[,valid] - Ideal[,valid]) * (1 - Y.matrix[,valid])), 1, sum)
    loss.matrix[, i] <- loss
    min.loss <- which(loss == min(loss))
    if(length(min.loss) != 1){
      n.tie[i] <- length(min.loss)
      min.loss <- sample(min.loss, 1, prob = rep(1/length(min.loss), length(min.loss)))
    }
    est.class <- c(est.class, min.loss)
  }
  est.pattern <- pattern[est.class, ]
  est.ideal <- Ideal[est.class, ]
  output <- list(alpha.est = est.pattern, est.ideal = est.ideal,
                 est.class = est.class, n.tie = n.tie, pattern = pattern,
                 loss.matrix = loss.matrix, method = method, Q = Q, Y = Y)
  class(output) <- "AlphaNP"
  return(output)
}
cdmTools.aggregateCol <- function(mX, ind){
  uniq <- unique(ind)
  N <- nrow(mX)
  Lj <- length(uniq)
  res <- matrix(0, nrow = N, ncol = Lj)
  for(l in 1:Lj){
    loc <- which(ind == l)
    res[,l] <- rowSums(mX[, loc, drop = FALSE])
  }
  return(res)
}
cdmTools.matchMatrix <- function(A, B){
  return(t(t(match(apply(B, 1, paste, collapse = ""), apply(A, 1, paste, collapse = "")))))
}
cdmTools.AlphaPermute <- function(dim){
  alpha <- matrix(c(0, 1), 2, 1)
  for (i in 1:(dim - 1)) {
    alpha <- rbind(alpha, alpha)
    alpha <- cbind(alpha, c(rep(0, 2^i), rep(1, 2^i)))
  }
  return(alpha)
}
cdmTools.partial_order2 <- function(Kjj, AlphaPattern = NULL){
  if(is.null(AlphaPattern)){
    alp <- GDINA::attributepattern(Kjj)
  } else {
    alp <- AlphaPattern
  }
  alp <- cbind(c(1:nrow(alp)), rowSums(alp), alp)
  out <- NULL
  for(k in 1:max(alp[, 2])){
    for(i in 1:sum(alp[, 2] == k - 1)){
      alpk_1 <- alp[alp[, 2] == k - 1, , drop = FALSE]
      alpk <- alp[alp[, 2] == k, , drop = FALSE]
      out <- rbind(out, cbind(alpk[(apply(alpk, 1, function(x){all(x - alpk_1[i, ] >= 0)})), 1], alpk_1[i, 1]))
    }
  }
  colnames(out) <- c("l", "s")
  return(out)
}
cdmTools.m2l <- function(m, remove = NA){
  if(is.na(remove)){
    lapply(seq_len(nrow(m)), function(i) m[i, !is.na(m[i,
    ])])
  } else {
    lapply(seq_len(nrow(m)), function(i) m[i, m[i, ] != remove])
  }
}
cdmTools.model.table <- function(){
    data.frame(model.char=c("LOGGDINA","LOGITGDINA","UDF", "GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM", "MSDINA","BUGDINO","SISM"),
               model.num=c(-3:8),
               linkf.num = c(3,2,-1,1,1,1,1,2,3,1,1,1),
               linkf.char = c("log","logit","UDF","identity","identity","identity","identity","logit","log","identity","identity","identity"),
               rule = c(0,0,-1,0,1,2,3,3,3,4,5,6))
}
cdmTools.model2numeric <- function(model, J = 1){
  x <- cdmTools.model.table()
  if(is.numeric(model)){
    if (J != 1 && length(model) != J)
      model <- rep(model, J)
  } else {
    M <- x$model.char
    if (J != 1 && length(model) != J)
      model <- rep(model, J)
    model <- match(toupper(model), M) - 4
  }
  model
}
cdmTools.model2rule.j <- function(model.j){
  x <- cdmTools.model.table()
  if (is.character(model.j)) {
    x$rule[which(x$model.char == model.j)]
  }
  else if (is.numeric(model.j)) {
    x$rule[which(x$model.num == model.j)]
  }
}
cdmTools.LikNR <- function(mpar, mX, vlogPrior, vgroup, mloc, weights, simplify = TRUE){
  .Call("_GDINA_LikNR", PACKAGE = "GDINA", mpar, mX, vlogPrior, vgroup, mloc, weights, simplify)
}
cdmTools.designM <- function(Kj, rule, AlphaPattern = NULL){
  .Call("_GDINA_designM", PACKAGE = "GDINA", Kj, rule, AlphaPattern)
}
cdmTools.item_latent_group <- function(Q, AlphaPattern = NULL){
  .Call("_GDINA_item_latent_group", PACKAGE = "GDINA", Q, AlphaPattern)
}
est.polarity <- function(polarity, Q, polarity.initial = 1e-4, polarity.prior = NULL){
  J <- nrow(polarity)
  init.parm <- lapply(1:J, function(x) rep(NA, 4))
  item.prior <- plyr::llply(2^rowSums(Q), function(x) matrix(NA, nrow = x, ncol = 2))
  names(init.parm) <- names(item.prior) <- sapply(1:J, function(x) paste0("Item ", x))
  if(is.null(polarity.prior)){polarity.prior <- list(c(1, 1), c(1, 1), c(1, 1))}
  for(j in 1:J){
    alpha.j <- expand.grid(lapply(1:sum(Q[j,]), function(x) c(0,1)))
    tmp <- (polarity[j,1]*(2*alpha.j[,1]-1) + polarity[j,2]*(-2*alpha.j[,2]+1) + 2)/4
    init.parm[[j]][tmp == 0] <- polarity.initial
    init.parm[[j]][tmp == 0.5] <- 0.5
    init.parm[[j]][tmp == 1] <- 1 - polarity.initial
    names(init.parm[[j]]) <- paste0("P(", apply(alpha.j, 1, paste, collapse = ""), ")")
    item.prior[[j]][tmp == 0,] <- polarity.prior[[1]]
    item.prior[[j]][tmp == 0.5,] <- polarity.prior[[2]]
    item.prior[[j]][tmp == 1,] <- polarity.prior[[3]]
    rownames(item.prior[[j]]) <- apply(alpha.j, 1, paste, collapse = "")
    colnames(item.prior[[j]]) <- c("alfa", "beta")
  }
  return(list(init.parm = init.parm, item.prior = item.prior))
}
GDINA.MJ <- function(dat, Q, verbose = 0, item.prior = NULL, catprob.parm = NULL, mono.constr = FALSE, conv.crit = 1e-4, maxitr = 2000, bound = 1e-4, model = "GDINA"){

  options('nloptr.show.inequality.warning'=FALSE)
  dat <- as.matrix(dat)
  Q <- as.matrix(Q)

  K <- ncol(Q)

  Kj <- rowSums(Q > 0)
  Lj <- 2^Kj
  N <- nrow(dat)
  J <- ncol(dat)
  AlphaPattern <- GDINA::attributepattern(K)
  parloc <- GDINA::LC2LG(Q=Q)

  ncat <- J
  model <- cdmTools.model2numeric(model, ncat)
  rule <- sapply(model, cdmTools.model2rule.j)

  reduced.LG <- my.item_latent_group(Q)

  if(any(model == -1)) stop("design.matrix must be provided for user-defined models.",call. = FALSE)
  DesignMatrices <-  vector("list", ncat)
  for(j in seq_len(ncat)) {
    if(model[j] == 6){
      DesignMatrices[[j]] <- GDINA::designmatrix(model = model[j],Qj = Q[which(Q[,1]==j),-c(1:2),drop=FALSE])
    }else if(rule[j] >= 0 & rule[j]<= 3){
      DesignMatrices[[j]] <- my.designM(Kj[j], rule[j], reduced.LG[[j]])
    }
  }

  L <- nrow(AlphaPattern)  # The number of latent classes
  if(is.null(item.prior)){item.prior <- plyr::llply(Lj,function(x)matrix(1,nrow = x,ncol = 2))}

  prior <- rep(1/L, L)
  logprior <- log(prior)

  # initial values

  if(is.null(catprob.parm)){
    item.parm <- list()
    for (j in seq_len(J)) {
      if(model[j] == 0){
        if (Kj[j] == 1) {
          item.parm[[j]] <- c(0.2, 0.8)
        } else if (Kj[j] == 2) {
          item.parm[[j]] <- c(0.2, 0.5, 0.5, 0.8)
        } else if (Kj[j] == 3) {
          item.parm[[j]] <- c(0.2, 0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.8)
        } else if (Kj[j] > 3){
          item.parm[[j]] <- c(0.2, rep(0.5,2^Kj[j]-2), 0.8)
        }
      } else if(model[j] == 1){ # pongo runif(.1, .3) como en STAN para DINA
        if (Kj[j] == 1) {
          item.parm[[j]] <- c(stats::runif(1, .1, .3), 1 - stats::runif(1, .1, .3))
        } else if (Kj[j] > 1) {
          initg <- stats::runif(1, .1, .3)
          item.parm[[j]] <- c(initg, rep(initg,2^Kj[j]-2), 1 - stats::runif(1, .1, .3))
        }
      } else if(model[j] == 2){
        if (Kj[j] == 1) {
          item.parm[[j]] <- c(0.2, 0.8)
        } else if (Kj[j] > 1) {
          item.parm[[j]] <- c(0.2, rep(0.8,2^Kj[j]-2), 0.8)
        }
      }
    }
    item.parm <- t(sapply(item.parm, "length<-", value = max(sapply(item.parm, length))))
  } else {
    item.parm <- t(sapply(catprob.parm, "length<-", value = max(sapply(catprob.parm, length))))
  }

  ConstrMatrix <- vector("list",J)
  ConstrPairs <- vector("list",J)

  for(j in seq_len(J)) {

    ConstrPairs[[j]] <- cdmTools.partial_order2(Kj[j])
    nctj <- nrow(ConstrPairs[[j]])
    tmp <- matrix(0,nctj,2^Kj[j])
    tmp[matrix(c(seq_len(nctj),ConstrPairs[[j]][,1]),ncol = 2)] <- 1
    tmp[matrix(c(seq_len(nctj),ConstrPairs[[j]][,2]),ncol = 2)] <- -1
    ConstrMatrix[[j]] <- tmp

  }

  itr <- 0L

  parm0 <- c(item.parm)
  success <- TRUE
  while(itr < maxitr)  {
    estep <- my.LikNR(as.matrix(item.parm), as.matrix(dat), as.matrix(logprior), rep(1,N),
                      as.matrix(parloc), rep(1,N), TRUE)

    Rg = estep$Rg
    Ng = estep$Ng

    correction = c(.0005, .001)

    for(j in 1:J){ #for each item

      modelj <- model[j]
      designMj=DesignMatrices[[j]]

      Nj=Ng[j,1:2^Kj[j]]
      Rj=Rg[j,1:2^Kj[j]]

      if(modelj==1 | modelj==2){
        rNj <- c(rep(sum(Nj[designMj[,2]==0],na.rm = TRUE),sum(designMj[,2]==0)),
                 rep(sum(Nj[designMj[,2]==1],na.rm = TRUE),sum(designMj[,2]==1)))
        rRj <- c(rep(sum(Rj[designMj[,2]==0],na.rm = TRUE),sum(designMj[,2]==0)),
                 rep(sum(Rj[designMj[,2]==1],na.rm = TRUE),sum(designMj[,2]==1)))
        if (any(rNj<correction[2])){
          rNj[which(rNj<correction[2])] <- rNj[which(rNj<correction[2])] + correction[2]
          rRj[which(rNj<correction[2])] <- rRj[which(rNj<correction[2])] + correction[1]
        }

        Nj <- rNj
        Rj <- rRj
      }

      r1 <- c(item.prior[[j]][,1]) - 1
      r2 <- c(item.prior[[j]][,2]) - 1
      n <- r1 + r2

      #EM and BM estimates
      Pj <- (Rj + r1)/(Nj + n)
      Pj[Pj <= bound] <- bound
      Pj[Pj >= 1 - bound] <- 1 - bound
      Pj[is.na(Pj)] <- bound

      if(verbose==1)
        cat("\nitem",j,"min(Nj)=",min(Nj + n))

      if(mono.constr&&any(c(ConstrMatrix[[j]]%*%Pj)<0)){

        obj <- function(x0){
          -1*sum(Rj*log(x0)+(Nj-Rj)*log(1-x0))-sum(r1*log(x0)+r2*log(1-x0))
        }
        dev <- function(x0){
          -1*(Rj + r1 - (Nj + n) * x0)/(x0-x0^2)
        }
        ineq <- function(x0){
          c(ConstrMatrix[[j]] %*% x0)
        }
        ineq.jac <- function(x0){
          ConstrMatrix[[j]]
        }
        x00 <- item.parm[j,1:(2^Kj[j])]

        x00[x00<bound] <- bound
        x00[x00>1-bound] <- 1-bound

        optims <- try(nloptr::slsqp(x0 = x00,fn=obj,
                                    gr = dev,
                                    hin=ineq,
                                    hinjac = ineq.jac,
                                    lower = rep(bound,length(x00)),
                                    upper = rep(1-bound,length(x00))),silent = TRUE)

        if(inherits(optims,"try-error")){
          warning(paste("Optimization failed for item",j),call. = FALSE)
          if(verbose==1)
            print(optims)
          return(list(success=FALSE))
        }

        item.parm[j,1:(2^Kj[j])] <- optims$par
      }else{
        item.parm[j,1:(2^Kj[j])] <- Pj
      }

    }
    prior <- c(exp(estep$logprior))
    prior <- prior/sum(prior)
    logprior <- log(prior)


    parm1 <- c(item.parm)

    maxchg = max(abs(parm1-parm0),na.rm = TRUE)

    parm0 <- parm1
    itr <- itr + 1

    if(is.infinite(estep$LL))
      stop("-2LL is not finite.",call. = FALSE)
    if(verbose==1L) {
      cat('\rIter =',itr,' Max. abs. change =',formatC(maxchg,digits = 5, format = "f"),
          ' Deviance  =',formatC(-2 * estep$LL,digits = 3, format = "f"),'                                                                                 ')
    }else if (verbose==2L) {
      cat('Iter =',itr,' Max. abs. change =',formatC(maxchg,digits = 5, format = "f"),
          ' Deviance  =',formatC(-2 * estep$LL,digits = 3, format = "f"),'                                                                                \n')
    }

    if(maxchg < conv.crit) break
  }
  estep <- my.LikNR(as.matrix(item.parm),
                    as.matrix(dat),
                    as.matrix(logprior),
                    rep(1,N),
                    as.matrix(parloc),
                    rep(1,N),
                    FALSE)

  EAP <- 1*((exp(estep$logpost) %*% AlphaPattern) > 0.5000)
  MAP <- AlphaPattern[max.col(estep$logpost),]
  options('nloptr.show.inequality.warning'=TRUE)
  list(catprob.parm = cdmTools.m2l(item.parm), posterior.prob = exp(estep$logprior),EAP=EAP, MAP=MAP,success=success,item.prior=item.prior, logpost = estep$logpost, loglik = estep$loglik, Q = Q)
}
relfit.GDINA.MJ <- function(fit, item.prior = NULL){
  N <- nrow(fit$EAP)
  lik.il <- exp(fit$loglik)
  lik.il <- t(sapply(1:N, function(i) lik.il[i,] * fit$posterior.prob[,1]))
  Deviance <- -2*sum(log(rowSums(lik.il)))
  Q <- fit$Q
  popu.parm <- 2^(ncol(Q)) - 1
  free.l <- c()
  for(j in 1:nrow(Q)){
    l <- sum(2^sum(Q[j,]))
    fixed.l <- 0
    if(!is.null(item.prior)){fixed.l <- sum(apply(item.prior[[j]] == 1, 1, all))}
    free.l <- c(free.l, l - fixed.l)
  }
  item.parm <- sum(free.l)
  npar <- popu.parm + item.parm
  AIC <- Deviance + 2 * npar
  BIC <- Deviance + log(N) * npar
  CAIC <- Deviance + (log(N) + 1) * npar
  SABIC <- Deviance + log((N + 2)/24) * npar
  return(list(Deviance = Deviance, npar = npar, AIC = AIC, BIC = BIC, CAIC = CAIC, SABIC = SABIC))
}
