#' Forced-choice data simulation based on the G-DINA model
#'
#' @description Simulate forced-choice (FC) responses based on the G-DINA model (de la Torre, 2011) and the FC-DCM (Huang, 2023).
#' This function accommodates FC responses to the \code{simGDINA} function from the \code{GDINA} package (Ma & de la Torre, 2020).
#'
#' @param N A \code{numeric} value indicating the sample size.
#' @param Q.items A binary \code{matrix} of dimensions \code{J} statements x \code{K} attributes indicating what statements measure what attributes.
#' @param n.blocks A \code{numeric} value indicating the number of forced-choice blocks.
#' @param polarity A \code{matrix} of dimensions \code{n.blocks} x 2 indicating whether each statament in each block is direct (1) or inverse (-1). Default is a matrix full of 1 (i.e., all statements are direct).
#' @param att A \code{matrix} of dimensions \code{N} individuals x \code{K} attributes indicating the attribute profile of each individual. Default is \code{NULL}, meaning that attribute profiles will be simulated based on the specifications listed on \code{GDINA.args} or \code{FCDCM.args}.
#' @param model Use the G-DINA model (\code{"GDINA"}) or the FC-DCM (\code{"FCDCM"}) as the generating model. Default is \code{"GDINA"}.
#' @param GDINA.args A \code{list} of arguments used if \code{model = "GDINA"}.
#' \describe{
#' \item{\code{GS}}{A \code{J} statements x 2 matrix indicating the guessing and slip parameter of each statement. Default is \code{NULL}.}
#' \item{\code{GS.items}}{Only used if \code{GDINA.args$GS = NULL}. A \code{vector} of length 2 indicating the minimum and maximum value for the random generation of guessing and slip parameters for each statement. Default is \code{c(1/3, 1/3).}}
#' \item{\code{AC}}{A \code{numeric} value indicating the attribute correlations in line with the multivariate normal threshold model (Chiu et al., 2009). Default is 0.}
#' \item{\code{AT}}{A \code{numeric} value indicating the attribute thresholds in line with the multivariate normal threshold model (Chiu et al., 2009). Default is 0.}
#' }
#' @param FCDCM.args A \code{list} of arguments used if \code{model = "FCDCM"}.
#' \describe{
#' \item{\code{d0}}{A \code{vector} of length 2 indicating the minimum and maximum value for the baseline probability for each FC block (see Huang, 2023). Default is \code{c(0.2, 0.2).}}
#' \item{\code{sd}}{A \code{vector} of length 2 indicating the minimum and maximum value for the statement utility parameters (see Huang, 2023). Default is \code{c(0.15, 0.15).}}
#' \item{\code{a}}{A \code{numeric} value indicating the minimum and maximum discrimination parameter for the higher-order model. Default is \code{c(0, 0).}}}
#' \item{\code{b}}{A \code{numeric} value indicating the location parameter for the higher-order model. Default is 0.}}}
#' }
#' @param seed Random number generation seed. Default is \code{NULL}.
#'
#' @return \code{simFCGDINA} returns an object of class \code{simFCGDINA}.
#' \describe{
#' \item{\code{dat}}{Generated FC responses (\code{matrix}).}
#' \item{\code{att}}{Generated attribute profiles (\code{matrix}).}
#' \item{\code{Q}}{Generated Q-matrix of FC blocks (\code{matrix}).}
#' \item{\code{LCprob}}{Generated block response probabilities for each latent class (\code{matrix}).}
#' \item{\code{item.pairs}}{Statements used in each FC block (\code{matrix}).}
#' \item{\code{q_att}}{Attribute measured by each statement as used by Huang (2023) (\code{matrix}).}
#' \item{\code{q_sta}}{Relative position of each statement as used by Huang (2023) (\code{matrix}).}
#' \item{\code{simGDINA}}{Object of class \code{simGDINA} (\code{list}).}
#' \item{\code{polarity}}{Polarity matrix indicating the direction of each statement in each block (\code{matrix}).}
#' \item{\code{GS}}{Generated guessing and slip parameter for each statement (\code{matrix}).}
#' }
#'
#' @author {Pablo Nájera, Universidad Pontificia Comillas}
#'
#' @references
#' Huang, H.-Y. (2023). Diagnostic Classification Model for Forced-Choice Items and Noncognitive Tests. \emph{Educational and Psychological Measurement}, \emph{83}(1), 146-180. https://doi.org/10.1177/00131644211069906
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R package for cognitive diagnosis modeling. \emph{Journal of Statistical Software}, \emph{93}(14). https://doi.org/10.18637/jss.v093.i14
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(GDINA)
#' set.seed(123)
#' Q.items <- do.call("rbind", replicate(5, diag(5), simplify = FALSE)) # Q-matrix for the unidimensional statements
#' GS <- cbind(runif(n = nrow(Q.items), min = 0.1, max = 0.3), runif(n = nrow(Q.items), min = 0.1, max = 0.3)) # Guessing and slip parameter for each statement
#' n.blocks <- 30 # Number of forced-choice blocks
#' polarity <- matrix(1, nrow = n.blocks, ncol = 2) # Block polarity (1 = direct statement; -1 = indirect statement)
#' sim <- simFCGDINA(N = 1000, Q.items, n.blocks = n.blocks, polarity = polarity, model = "GDINA", GDINA.args = list(GS = GS), seed = 123)
#' }
simFCGDINA <- function(N, Q.items, n.blocks = NULL, polarity = NULL, att = NULL, model = "GDINA",
                       GDINA.args = list(GS = NULL, GS.items = c(1/3, 1/3), AC = 0, AT = 0),
                       FCCDM.args = list(d0 = c(0.2, 0.2), sd = c(0.15, 0.15), a = c(0, 0), b = 0),
                       seed = NULL){

  library(GDINA)
  if(!is.null(seed)){set.seed(seed)}
  J.items <- nrow(Q.items)
  K <- ncol(Q.items)
  if(is.null(n.blocks)){n.blocks <- J.items / 2}
  J <- n.blocks
  if(is.null(GDINA.args$GS.items)){GDINA.args$GS.items <- c(1/3, 1/3)}
  if(is.null(GDINA.args$AC)){GDINA.args$AC <- 0}
  if(is.null(GDINA.args$AT)){GDINA.args$AT <- 0}
  if(is.null(FCCDM.args$d0)){FCCDM.args$d0 <- c(0.2, 0.2)}
  if(is.null(FCCDM.args$sd)){FCCDM.args$sd <- c(0.15, 0.15)}
  if(is.null(FCCDM.args$a)){FCCDM.args$a <- c(0, 0)}
  if(is.null(FCCDM.args$b)){FCCDM.args$b <- 0}
  if(model == "FCCDM" & !is.null(polarity)){
    polarity <- matrix(1, nrow = J, ncol = 2)
    warning("FCCDM only supports homopolar blocks. Polarity has been changed accordingly.")
  }
  if(is.null(polarity)){polarity <- matrix(1, nrow = J, ncol = 2)}

  kj <- apply(Q.items, 1, function(x) which(x == 1))
  item.pairs <- matrix(NA, nrow = J, ncol = 2)
  Q <- matrix(0, nrow = J, ncol = K)
  q_att <- q_sta <- matrix(NA, nrow = J, ncol = 2)
  item.pool <- 1:J.items
  while(!all(rowSums(Q) == 2)){
    j <- which(rowSums(Q) == 0)[1]
    j1 <- sample(item.pool, size = 1)
    k1 <- kj[j1]
    pre.j2 <- (1:J.items)[kj != k1]
    if(length(intersect(item.pool, pre.j2)) == 0){
      if(length(item.pool <= 2)){
        item.pool <- 1:J.items
        next
      } else {
        next
      }
    }
    if(j1 %in% item.pairs){
      prev.j2 <- which(j1 == item.pairs, arr.ind = TRUE)[1,]
      prev.j2 <- setdiff(item.pairs[prev.j2[1],], j1)
      pre.j2 <- setdiff(pre.j2, prev.j2)
    }
    pre.j2 <- intersect(item.pool, pre.j2)
    if(length(pre.j2) == 0){
      if(length(item.pool <= 3)){
        item.pool <- 1:J.items
        next
      } else {
        next
      }
    }
    if(length(pre.j2) > 1){pre.j2 <- sample(pre.j2, size = 1)}
    j2 <- pre.j2
    k2 <- kj[j2]
    item.pairs[j,] <- c(j1, j2)
    Q[j, c(k1, k2)] <- 1
    q_sta[j,] <- c(which(j1 == which(Q.items[,which(Q.items[j1,] == 1)] == 1)), which(j2 == which(Q.items[,which(Q.items[j2,] == 1)] == 1)))
    q_att[j,] <- c(k1, k2)
    item.pool <- setdiff(item.pool, c(j1, j2))
    q_sta[j,] <- q_sta[j,][order(q_att[j,])]
    item.pairs[j,] <- item.pairs[j,][order(q_att[j,])]
    q_att[j,] <- q_att[j,][order(q_att[j,])]
    if(length(item.pool) <= 1){item.pool <- 1:J.items}
  }

  if(model == "GDINA"){
    if(is.null(GDINA.args$GS)){
      GS.items <- cbind(runif(J.items, GDINA.args$GS.items[1], GDINA.args$GS.items[2]), runif(J.items, GDINA.args$GS.items[1], GDINA.args$GS.items[2]))
    } else {
      GS.items <- GDINA.args$GS
    }
    catprob.parm <- list()
    for(j in 1:J){
      p0_1 <- GS.items[item.pairs[j, 1], 1]
      p1_1 <- 1 - GS.items[item.pairs[j, 1], 2]
      p0_2 <- GS.items[item.pairs[j, 2], 1]
      p1_2 <- 1 - GS.items[item.pairs[j, 2], 2]
      if(polarity[j, 1] == -1){p0_1 <- 1 - p0_1; p1_1 <- 1 - p1_1}
      if(polarity[j, 2] == -1){p0_2 <- 1 - p0_2; p1_2 <- 1 - p1_2}
      Py_00 <- ((1-p0_2)*p0_1)/(((1-p0_2)*p0_1)+(p0_2*(1-p0_1)))
      Py_10 <- ((1-p0_2)*p1_1)/(((1-p0_2)*p1_1)+(p0_2*(1-p1_1)))
      Py_01 <- ((1-p1_2)*p0_1)/(((1-p1_2)*p0_1)+(p1_2*(1-p0_1)))
      Py_11 <- ((1-p1_2)*p1_1)/(((1-p1_2)*p1_1)+(p1_2*(1-p1_1)))
      catprob.parm[[j]] <- c(Py_00, Py_10, Py_01, Py_11)
      names(catprob.parm[[j]]) <- c("P(00)", "P(10)", "P(01)", "P(11)")
    }
    if(is.null(att)){
      cutoffs <- seq(-GDINA.args$AT, GDINA.args$AT, length.out = K)
      m <- rep(0, K)
      vcov <- matrix(GDINA.args$AC, K, K); diag(vcov) <- 1
      sim <- simGDINA(N, Q, catprob.parm = catprob.parm, att.dist = "mvnorm", mvnorm.parm = list(mean = m, sigma = vcov, cutoffs = cutoffs))
    } else {
      sim <- simGDINA(N, Q, catprob.parm = catprob.parm, attribute = att)
    }
    dat <- sim$dat
    att <- sim$attribute
    return(list(dat = dat, att = att, Q = Q, LCprob = sim$LCprob.parm, item.pairs = item.pairs, q_att = q_att, q_sta = q_sta, simGDINA = sim, polarity = polarity, GS = GS.items))
  }

  if(model == "FCCDM"){
    sd <- runif(J.items, FCCDM.args$sd[1], FCCDM.args$sd[2])
    d0 <- runif(J, FCCDM.args$d0[1], FCCDM.args$d0[2])
    theta <- rnorm(N, 0, 1)
    a <- rlnorm(K, meanlog = FCCDM.args$a[1], sdlog = FCCDM.args$a[2])
    b <- seq(-FCCDM.args$b, FCCDM.args$b, length.out = K)
    if(is.null(att)){
      mp <- att <- matrix(NA, nrow = N, ncol = K)
      for(n in 1:N){
        mp[n,] <- round(exp(a * (theta[n] - b))/(1 + exp(a * (theta[n] - b))), 5)
        att[n,] <- rbinom(K, 1, prob = mp[n,])
        if(anyNA(att[n,])){
          if(sum(mp[n,] >= 1) > 0){att[n, which(mp[n,] > 1)] <- 1}
          if(sum(mp[n,] <= 0) > 0){att[n, which(mp[n,] < 0)] <- 0}
        }
      }
    }
    L <- 2^K
    LCprob.parm <- matrix(NA, nrow = L, ncol = J)
    for(l in 1:L){
      for(j in 1:J){
        LCprob.parm[l, j] <- d0[j] +
          (0.5 - d0[j]) * (attributepattern(K)[l, q_att[j, 1]] >= attributepattern(K)[l, q_att[j, 2]]) + #GDINA::attributepatern
          (sd[item.pairs[j, 1]] * attributepattern(K)[l, q_att[j, 1]] + sd[item.pairs[j, 2]] * (1 - attributepattern(K)[l, q_att[j, 2]])) * (attributepattern(K)[l, q_att[j, 1]] > attributepattern(K)[l, q_att[j, 2]])
      }
    }
    dat <- matrix(NA, nrow = N, ncol = J)
    for(n in 1:N){
      lc <- which(sapply(1:L, function(l) all(attributepattern(K)[l,] == att[n,])))
      dat[n,] <- rbinom(J, 1, prob = LCprob.parm[lc,])
    }
    return(list(dat = dat, att = att, Q = Q, LCprob = LCprob.parm, d0 = d0, sd = sd, theta = theta, a = a, b = b, item.pairs = item.pairs, q_att = q_att, q_sta = q_sta, polarity = polarity))
  }
}
