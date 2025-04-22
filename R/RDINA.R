#' Restricted DINA model
#'
#' @description Estimation of the \emph{restricted deterministic input, noisy "and" gate} model (R-DINA; Nájera et al., 2023).
#' In addition to the non-compensatory (i.e., conjunctive) condensation rule of the DINA model, the compensatory (i.e., disjunctive) rule of the DINO model can be also applied (i.e., R-DINO model).
#' The R-DINA/R-DINO model should be only considered for applications involving very small sample sizes (N < 100; Nájera et al., 2023), and model fit evaluation and comparison with competing models (e.g., DINA/DINO, G-DINA) is highly recommended.
#'
#' @param dat A \emph{N} individuals x \emph{J} items (\code{matrix} or \code{data.frame}). Missing values need to be coded as \code{NA}. Caution is advised if missing data are present.
#' @param Q A \emph{J} items x \emph{K} attributes Q-matrix (\code{matrix} or \code{data.frame}).
#' @param gate Either a conjunctive (\code{"AND"}) or disjunctive (\code{"OR"}) condensation rule to estimate the RDINA or RDINO model, respectively. Default is \code{"AND"}.
#' @param att.prior A 2^\emph{K} attributes vector containing the prior distribution for each latent class. The sum of all elements does not have to be equal to 1, since the vector will be normalized. Default is \code{NULL}, which is a uniform prior distribution.
#' @param est Use the Brent's method (\code{"Brent"}) or the expectation-maximization algorithm (\code{"EM"}) to estimate the model? Default is \code{"Brent"}, since it is faster and both algorithms are virtually equivalent for the RDINA/RDINO model.
#' @param tau.alpha Attribute profile estimator (either \code{"MAP"}, \code{"EAP"}, or \code{"MLE"}) used to calculate the estimated classification accuracy as done with the \code{CA} function of the \code{GDINA} package (Ma & de la Torre, 2020).
#' @param emp.bayes Use empirical Bayes estimation for structural parameters. Default is \code{FALSE}.
#' @param boot Use bootstrapping to increase robustness in posterior probabilities estimation. Default is \code{FALSE}.
#' @param n.boots Number of bootstrapping samples. Default is 500.
#' @param n.cores Number of CPU processors to speed up computation when bootstrapping is used. Default is 1.
#' @param maxitr Maximum number of iterations. Default is 1000.
#' @param conv.crit Convergence criterion regarding the maximum absolute change in either the phi parameter estimate or the marginal posterior probabilities of attribute mastery. Default is 0.0001.
#' @param init.phi Initial value for the phi parameter. Default is 0.2.
#' @param bound.p Lowest value for probability estimates (highest would be 1 - bound.p). Default is 1e-06.
#' @param verbose Print information after each iteration. Default is \code{TRUE}.
#' @param seed Random number generation seed (e.g., to solve ties in case they occur with MLE or MAP estimation). Default is \code{NULL}, which means that no specific seed is used.
#'
#' @return \code{RDINA} returns an object of class \code{RDINA}.
#' \describe{
#' \item{\code{MLE}}{Estimated attribute profiles with the MLE estimator (\code{matrix}).}
#' \item{\code{MAP}}{Estimated attribute profiles with the MAP estimator (\code{matrix}).}
#' \item{\code{EAP}}{Estimated attribute profiles with the EAP estimator (\code{matrix}).}
#' \item{\code{phi}}{Phi parameter estimate (\code{numeric}).}
#' \item{\code{post.probs}}{A (\code{list}) containing the estimates of the posterior probability of each examinee in each latent class (\code{pp}), marginal posterior probabilities of attribute mastery (\code{mp}), and posterior probability of each latent class (\code{lp}).}
#' \item{\code{likelihood}}{A (\code{list}) containing the likelihood of each examinee in each latent class (\code{lik_il}) and the model log-likelihood (\code{logLik}).}
#' \item{\code{test.fit}}{Relative model fit indices (\code{list}).}
#' \item{\code{class.accu}}{A (\code{list}) containing the classification accuracy estimates at the test-level (\code{tau}), latent class-level (\code{tau_l}), and attribute-level (\code{tau_k}).}
#' \item{\code{specifications}}{Function call specifications (\code{list}).}
#' }
#'
#' @author {Pablo Nájera, Universidad Pontificia Comillas}
#'
#' @references
#' Ma, W., & de la Torre, J. (2020). GDINA: An R package for cognitive diagnosis modeling. \emph{Journal of Statistical Software}, \emph{93}(14). https://doi.org/10.18637/jss.v093.i14
#'
#' Nájera, P., Abad, F. J., Chiu, C.-Y., & Sorrel, M. A. (2023). The Restricted DINA model: A Comprehensive Cognitive Diagnostic Model for Classroom-Level Assessments. \emph{Journal of Educational and Behavioral Statistics}.
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
#' GS <- data.frame(guessing = rep(0.2, J), slip = rep(0.2, J))
#' sim <- simGDINA(20, Q, GS, model = "DINA")
#' simdat <- sim$dat # Simulated data
#' simatt <- sim$attribute # Generating attributes
#' fit.RDINA <- RDINA(simdat, Q) # Apply the GNPC method
#' ClassRate(fit.RDINA$EAP, simatt) # Check classification accuracy
#' }
RDINA <- function(dat, Q, gate = "AND", att.prior = NULL, est = "Brent", tau.alpha = "MAP", emp.bayes = FALSE, boot = FALSE, n.boots = 500,
                  n.cores = 1, maxitr = 1000, conv.crit = 1e-04, init.phi = 0.2, bound.p = 1e-06, verbose = TRUE, seed = NULL){

  #-------------------
  # Arguments control
  #-------------------

  if(!gate %in% c("AND", "OR")){stop("Error in RDINA: gate must be 'AND' or 'OR'")}
  if(!est %in% c("Brent", "EM")){stop("Error in RDINA: est must be 'Brent' or 'EM'")}
  if(!tau.alpha %in% c("MLE", "MAP", "EAP")){stop("Error in RDINA: tau.alpha must be 'MLE', 'MAP', or 'EAP'")}
  if(!is.null(att.prior)){if(!is.vector(att.prior)){stop("Error in RDINA: att.prior must be either NULL or a vector")}}
  if(!is.logical(emp.bayes)){stop("Error in RDINA: emp.bayes must be logical")}
  if(!is.logical(boot)){stop("Error in RDINA: boot must be logical")}
  if(!is.numeric(n.boots)){stop("Error in RDINA: n.boots must be numeric")}
  if(n.boots <= 0){stop("Error in RDINA: n.boots must be higher than 0")}
  if(!is.numeric(n.cores)){stop("Error in RDINA: n.cores must be numeric")}
  if(n.cores <= 0){stop("Error in RDINA: n.cores must be higher than 0")}
  if(n.cores > parallel::detectCores()){stop(paste0("Error in RDINA: n.cores cannot suprass the number of maximum CPU cores: ", parallel::detectCores()))}
  if(!boot & n.cores > 1){warning("Warning in RDINA: n.cores equal to 1 will be used as no bootstrapping is used")}
  if(!is.numeric(maxitr)){stop("Error in RDINA: maxitr must be numeric")}
  if(maxitr < 1){stop("Error in RDINA: maxitr must be equal to or higher than 1")}
  if(!is.numeric(bound.p)){stop("Error in RDINA: conv.crit must be numeric")}
  if(bound.p <= 0){stop("Error in RDINA: conv.crit must be higher than 0")}
  if(!is.numeric(init.phi)){stop("Error in RDINA: init.phi must be numeric")}
  if(init.phi <= 0){stop("Error in RDINA: init.phi must be higher than 0")}
  if(!is.numeric(bound.p)){stop("Error in RDINA: bound.p must be numeric")}
  if(bound.p <= 0){stop("Error in RDINA: bound.p must be higher than 0")}
  if(!is.logical(verbose)){stop("Error in RDINA: verbose must be logical")}
  if(!is.null(seed)){if(!is.vector(seed)){stop("Error in RDINA: seed must be either NULL or a vector")}}

  #-----------------------
  # Information gathering
  #-----------------------

  if(is.null(seed)){seed <- sample(1:100000, size = 1)}
  set.seed(seed)
  RDINA.call <- match.call()
  N <- nrow(dat)
  J <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K
  lclass <- GDINA::attributepattern(K)
  lgroups <- apply(lclass, 1, paste, collapse = "")
  if(is.null(att.prior)){att.prior <- rep(1/L, L)}
  att.prior <- att.prior/sum(att.prior)

  #-------
  # Brent
  #-------

  if(est == "Brent"){

    # Without bootstrapping

    NPC <- cdmTools.AlphaNP(dat, Q, gate, method = "Hamming")
    match_lclass <- match(apply(lclass, 1, paste, collapse = ""), apply(cdmTools.AlphaPermute(K), 1, paste, collapse = ""))
    alpha.est <- NPC$alpha.est
    dist.li <- NPC$loss.matrix[match_lclass, ]

    if(!boot){
      phi <- stats::optimize(f = phi.ML, dist = dist.li, J = J, posterior = FALSE, att.prior = att.prior, interval = c(bound.p, 0.5 - bound.p), maximum = TRUE)$maximum
      lik.il <- t(phi^dist.li * (1 - phi)^(J - dist.li))
      marg.lik.il <- t(t(lik.il) * att.prior)
      pp <- t(apply(marg.lik.il, 1, function(i) i/sum(i)))
      lp <- colMeans(pp)
      lp[lp == 0] <- bound.p; lp[lp == 1] <- 1 - bound.p; lp <- lp / sum(lp)
      r <- 2
      if(emp.bayes){
        while(r <= maxitr){
          phi_back <- phi
          pp_back <- pp
          phi <- stats::optimize(f = phi.ML, dist = dist.li, J = J, posterior = FALSE, att.prior = lp, interval = c(bound.p, 0.5 - bound.p), maximum = TRUE)$maximum
          lik.il <- t(phi^dist.li * (1 - phi)^(J - dist.li))
          marg.lik.il <- t(t(lik.il) * lp)
          pp <- t(apply(marg.lik.il, 1, function(i) i/sum(i)))
          lp <- colMeans(pp)
          lp[lp == 0] <- bound.p; lp[lp == 1] <- 1 - bound.p; lp <- lp / sum(lp)
          phi_change_tmp <- phi - phi_back
          mp_change_tmp <- pp - pp_back
          max_abs_change <- max(abs(c(phi_change_tmp, mp_change_tmp)))
          if(verbose){cat("\r", "Iter =", r, " Max. abs. change =", round(max_abs_change, 5), "        ")}
          if(max_abs_change < conv.crit){break}
          r <- r + 1
        }
      }
    }

    # With bootstrapping

    if(boot){

      if(n.cores == 1){
        b.pp <- array(NA, dim = c(N, L, n.boots))
        for(b in 1:n.boots){
          set.seed(b + seed)
          dat.b <- dat[sample(1:N, size = N, replace = TRUE),]
          NPC.b <- cdmTools.AlphaNP(dat.b, Q, gate, method = "Hamming")
          match_lclass.b <- match(apply(lclass, 1, paste, collapse = ""), apply(cdmTools.AlphaPermute(K), 1, paste, collapse = ""))
          alpha.est.b <- NPC.b$alpha.est
          dist.li.b <- NPC.b$loss.matrix[match_lclass, ]
          phi <- stats::optimize(f = phi.ML, dist = dist.li.b, J = J, posterior = FALSE, att.prior = att.prior, interval = c(bound.p, 0.5 - bound.p), maximum = TRUE)$maximum
          lik.il <- t(phi^dist.li.b * (1 - phi)^(J - dist.li.b))
          marg.lik.il <- t(t(lik.il) * att.prior)
          pp <- t(apply(marg.lik.il, 1, function(i) i/sum(i)))
          lp <- colMeans(pp)
          lp[lp == 0] <- bound.p; lp[lp == 1] <- 1 - bound.p; lp <- lp / sum(lp)
          r <- 2
          if(emp.bayes){
            while(r <= maxitr){
              phi_back <- phi
              pp_back <- pp
              phi <- stats::optimize(f = phi.ML, dist = dist.li.b, J = J, posterior = FALSE, att.prior = lp, interval = c(bound.p, 0.5 - bound.p), maximum = TRUE)$maximum
              lik.il <- t(phi^dist.li.b * (1 - phi)^(J - dist.li.b))
              marg.lik.il <- t(t(lik.il) * lp)
              pp <- t(apply(marg.lik.il, 1, function(i) i/sum(i)))
              lp <- colMeans(pp)
              lp[lp == 0] <- bound.p; lp[lp == 1] <- 1 - bound.p; lp <- lp / sum(lp)
              phi_change_tmp <- phi - phi_back
              mp_change_tmp <- pp - pp_back
              max_abs_change <- max(abs(c(phi_change_tmp, mp_change_tmp)))
              if(verbose){cat("\r", "Bootstrap sample =", b, " Iter =", r, " Max. abs. change =", round(max_abs_change, 5), "        ")}
              if(max_abs_change < conv.crit){break}
              r <- r + 1
            }
          }
          b.catprob <- list()
          for(j in 1:J){
            if(gate == "AND"){b.catprob[[j]] <- c(rep(phi, 2^sum(Q[j,]) - 1), 1 - phi)}
            if(gate == "OR"){b.catprob[[j]] <- c(phi, rep(1 - phi, 2^sum(Q[j,]) - 1))}
          }
          b.rdina <- GDINA::GDINA(dat = dat, Q = Q, att.prior = lp, catprob.parm = b.catprob, control = list(maxitr = 0))
          b.pp[,,b] <- as.matrix(exp(b.rdina$technicals$logposterior.i))
        }
        pp <- apply(b.pp, 1:2, mean)
        lp <- colMeans(pp)
        lp[lp == 0] <- bound.p; lp[lp == 1] <- 1 - bound.p; lp <- lp / sum(lp)
        phi <- stats::optimize(f = phi.ML, dist = dist.li, J = J, posterior = FALSE, att.prior = lp, interval = c(bound.p, 0.5 - bound.p), maximum = TRUE)$maximum
        lik.il <- t(phi^dist.li * (1 - phi)^(J - dist.li))
        marg.lik.il <- t(t(lik.il) * lp)
      }

      if(n.cores > 1){
        cl <- parallel::makeCluster(n.cores, type = "SOCK")
        doSNOW::registerDoSNOW(cl)
        if(verbose){
          cat("Bootstrapping Progress:", "\n")
          pb <- utils::txtProgressBar(max = n.boots, style = 3)
          progress <- function(n) utils::setTxtProgressBar(pb, n)
          opts <- list(progress = progress)
        } else {
          opts <- NULL
        }
        b.pp = foreach::foreach(b = 1:n.boots,
                                .options.snow = opts,
                                .combine = 'list', .multicombine = TRUE, .packages = "GDINA") %dopar% {
                                  set.seed(b + seed)
                                  dat.b <- dat[sample(1:N, size = N, replace = TRUE),]
                                  NPC.b <- cdmTools.AlphaNP(dat.b, Q, gate, method = "Hamming")
                                  match_lclass.b <- match(apply(lclass, 1, paste, collapse = ""), apply(cdmTools.AlphaPermute(K), 1, paste, collapse = ""))
                                  alpha.est.b <- NPC.b$alpha.est
                                  dist.li.b <- NPC.b$loss.matrix[match_lclass, ]
                                  phi <- stats::optimize(f = phi.ML, dist = dist.li.b, J = J, posterior = FALSE, att.prior = att.prior, interval = c(bound.p, 0.5 - bound.p), maximum = TRUE)$maximum
                                  lik.il <- t(phi^dist.li.b * (1 - phi)^(J - dist.li.b))
                                  marg.lik.il <- t(t(lik.il) * att.prior)
                                  pp <- t(apply(marg.lik.il, 1, function(i) i/sum(i)))
                                  lp <- colMeans(pp)
                                  lp[lp == 0] <- bound.p; lp[lp == 1] <- 1 - bound.p; lp <- lp / sum(lp)
                                  r <- 2
                                  if(emp.bayes){
                                    while(r <= maxitr){
                                      phi_back <- phi
                                      pp_back <- pp
                                      phi <- stats::optimize(f = phi.ML, dist = dist.li.b, J = J, posterior = FALSE, att.prior = lp, interval = c(bound.p, 0.5 - bound.p), maximum = TRUE)$maximum
                                      lik.il <- t(phi^dist.li.b * (1 - phi)^(J - dist.li.b))
                                      marg.lik.il <- t(t(lik.il) * lp)
                                      pp <- t(apply(marg.lik.il, 1, function(i) i/sum(i)))
                                      lp <- colMeans(pp)
                                      lp[lp == 0] <- bound.p; lp[lp == 1] <- 1 - bound.p; lp <- lp / sum(lp)
                                      phi_change_tmp <- phi - phi_back
                                      mp_change_tmp <- pp - pp_back
                                      max_abs_change <- max(abs(c(phi_change_tmp, mp_change_tmp)))
                                      if(max_abs_change < conv.crit){break}
                                      r <- r + 1
                                    }
                                  }
                                  b.catprob <- list()
                                  for(j in 1:J){
                                    if(gate == "AND"){b.catprob[[j]] <- c(rep(phi, 2^sum(Q[j,]) - 1), 1 - phi)}
                                    if(gate == "OR"){b.catprob[[j]] <- c(phi, rep(1 - phi, 2^sum(Q[j,]) - 1))}
                                  }
                                  b.rdina <- GDINA::GDINA(dat = dat, Q = Q, att.prior = lp, catprob.parm = b.catprob, control = list(maxitr = 0))
                                  b.pp <- as.matrix(exp(b.rdina$technicals$logposterior.i))
                                  return(b.pp)
                                }
        parallel::stopCluster(cl)
        b.pp <- array(unlist(b.pp), dim = c(N, L, n.boots))
        b.pp <- lapply(1:n.boots, function(i) b.pp[,,i])
        pp <- Reduce("+", b.pp) / length(b.pp)
        lp <- colMeans(pp)
        lp[lp == 0] <- bound.p; lp[lp == 1] <- 1 - bound.p; lp <- lp / sum(lp)
        phi <- stats::optimize(f = phi.ML, dist = dist.li, J = J, posterior = FALSE, att.prior = lp, interval = c(bound.p, 0.5 - bound.p), maximum = TRUE)$maximum
        lik.il <- t(phi^dist.li * (1 - phi)^(J - dist.li))
        marg.lik.il <- t(t(lik.il) * lp)
      }

    }

  }

  #----
  # EM
  #----

  if(est == "EM"){
    phi <- init.phi
    mp_change <- c()
    r <- 0
    lp <- att.prior
    lp[lp == 0] <- bound.p; lp[lp == 1] <- 1 - bound.p; lp <- lp / sum(lp)
    for(R in 1:maxitr){
      r <- r + 1
      eta.lj <- switch(gate,
                       AND = t(apply(lclass, 1, function(a) apply(Q, 1, function(q) prod(a^q)))),
                       OR = t(apply(lclass, 1, function(a) apply(Q, 1, function(q) 1 - prod((1 - a)^q)))))
      jl <- t(sapply(1:J, function(j) phi^(1 - eta.lj[, j]) * (1 - phi)^(eta.lj[, j])))
      colnames(jl) <- paste0("L", 1:L)
      rownames(jl) <- paste0("J", 1:J)
      lik.il <- matrix(NA, nrow = N, ncol = L, dimnames = list(paste0("N", 1:N), paste0("L", 1:L)))
      for(l in 1:L){
        jl.l <- matrix(rep(jl[, l], N), nrow = N, byrow = T)
        lik.il[, l] <- apply(jl.l^dat * (1 - jl.l)^(1 - dat), 1, prod)
      }
      mp_change_tmp <- c()
      if(r > 1){pp_back <- pp}
      if(emp.bayes){
        marg.lik.il <- t(t(lik.il) * lp)
      } else {
        marg.lik.il <- t(t(lik.il) * att.prior)
      }
      pp <- t(apply(marg.lik.il, 1, function(i) i/sum(i)))
      lp <- colMeans(pp)
      lp[lp == 0] <- bound.p; lp[lp == 1] <- 1 - bound.p; lp <- lp / sum(lp)
      if(r > 1){mp_change_tmp <- c(mp_change_tmp, pp - pp_back)}
      l_pp <- colSums(pp)
      IR <- matrix(data = NA, nrow = J, ncol = 4, dimnames = list(paste("J", 1:J, sep = ""), c("I0", "R0", "I1", "R1")))
      for(j in 1:J){
        i0 <- i1 <- r0 <- r1 <- 0
        for(l in 1:L){
          if(t(eta.lj)[j, l] == 0){
            i0 <- i0 + l_pp[l]
          } else {
            i1 <- i1 + l_pp[l]
          }
          r0 <- pp[c(as.double(which(dat[, j] == 1))), c(as.double(which(t(eta.lj)[j, ] == 0)))]
          r1 <- pp[c(as.double(which(dat[, j] == 1))), c(as.double(which(t(eta.lj)[j, ] == 1)))]
        }
        IR[j, 1] <- i0
        IR[j, 2] <- sum(r0)
        IR[j, 3] <- i1
        IR[j, 4] <- sum(r1)
      }
      phi_change_tmp <- phi - sum(IR[, 2], (IR[, 3] - IR[, 4]))/sum(IR[, 1], IR[, 3])
      phi <- sum(IR[, 2], (IR[, 3] - IR[, 4]))/sum(IR[, 1], IR[, 3])
      max_abs_change <- max(abs(c(phi_change_tmp, mp_change_tmp)))
      if(verbose){cat("\r", "Iter =", r, " Max. abs. change =", round(max_abs_change, 5), "        ")}
      if(max_abs_change < conv.crit){break}
    }
  }

  #---------------------
  # Extract information
  #---------------------

  log.cond.lik <- sum(log(apply(lik.il, 1, max)))
  log.marg.lik <- sum(log(apply(marg.lik.il, 1, sum)))
  mp <- pp %*% lclass
  lp <- colMeans(pp)
  lp[lp == 0] <- bound.p; lp[lp == 1] <- 1 - bound.p; lp <- lp / sum(lp)
  colnames(lik.il) <- colnames(pp) <- apply(lclass, 1, paste, collapse = "")
  colnames(mp) <- paste0("A", 1:K)
  names(lp) <- apply(lclass, 1, paste, collapse = "")
  MLE.est <- cbind(as.data.frame(lclass[sapply(apply(marg.lik.il, 1, function(x) which(x == max(x))), function(y) ifelse(length(y) > 1, sample(y, size = 1), y)), ]), MM = ifelse(rowSums(t(apply(marg.lik.il, 1, function(i) i == max(i)))) == 1, FALSE, TRUE))
  MAP.est <- cbind(as.data.frame(lclass[sapply(apply(pp, 1, function(x) which(x == max(x))), function(y) ifelse(length(y) > 1, sample(y, size = 1), y)), ]), MM = ifelse(rowSums(t(apply(pp, 1, function(i) i == max(i)))) == 1, FALSE, TRUE))
  EAP.est <- ifelse(mp >= 0.5, 1, 0)
  tau.est <- as.matrix(switch(tau.alpha, MLE = MLE.est, MAP = MAP.est, EAP = EAP.est)[, 1:K])
  Deviance <- -2 * log.marg.lik
  npar <- ifelse(emp.bayes, 1 + (L - 1), 1)
  AIC <- Deviance + 2 * npar
  BIC <- Deviance + log(N) * npar
  CAIC <- Deviance + (log(N) + 1) * npar
  SABIC <- Deviance + log((N + 2)/24) * npar
  gr <- cdmTools.matchMatrix(lclass, tau.est)
  pseudo.gr <- setdiff(seq(nrow(lclass)), unique(gr))
  gr <- c(gr, pseudo.gr)
  lab <- apply(lclass, 1, paste0, collapse = "")
  post <- cbind(t(pp), matrix(0, nrow(lclass), length(pseudo.gr)))
  CCM <- cdmTools.aggregateCol(post, gr)/c(nrow(tau.est) * lp)
  tau_c <- diag(CCM)
  tau <- sum(tau_c * c(lp))
  tau_k <- colMeans(tau.est * mp + (1 - tau.est) * (1 - mp))
  names(tau_c) <- rownames(CCM) <- colnames(CCM) <- lab
  if(is.null(colnames(Q))){colnames(Q) <- paste0("A", 1:K)}
  if(is.null(rownames(Q))){rownames(Q) <- paste0("J", 1:J)}

  #-----------------
  # Return outcomes
  #-----------------

  specs <- list(dat = dat, Q = Q, gate = gate, att.prior = att.prior, est = est, tau.alpha = tau.alpha,
                maxitr = maxitr, conv.crit = conv.crit, init.phi = init.phi, verbose = verbose, seed = seed)
  res <- list(MLE = MLE.est, MAP = MAP.est, EAP = EAP.est,
              phi = phi, post.probs = list(pp = pp, mp = mp, lp = lp),
              likelihood = list(lik_il = lik.il, logLik = log.marg.lik),
              test.fit = list(Deviance = Deviance, npar = npar, AIC = AIC, BIC = BIC, CAIC = CAIC, SABIC = SABIC),
              class.accu = list(tau = tau, tau_l = tau_c, tau_k = tau_k, CCM = CCM), specifications = specs,
              call = RDINA.call)
  class(res) <- "RDINA"
  return(res)
}
