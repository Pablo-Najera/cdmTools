% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FCGDINA.R
\name{FCGDINA}
\alias{FCGDINA}
\title{G-DINA model for forced-choice blocks}
\usage{
FCGDINA(
  dat,
  Q,
  polarity = NULL,
  polarity.initial = 1e-04,
  att.dist = "saturated",
  att.prior = NULL,
  verbose = 1,
  higher.order = list(),
  catprob.parm = NULL,
  control = list()
)
}
\arguments{
\item{dat}{A \emph{N} individuals x \emph{J} items (\code{matrix} or \code{data.frame}). Missing values need to be coded as \code{NA}. Caution is advised if missing data are present.}

\item{Q}{A \emph{F} blocks x \emph{K} attributes Q-matrix (\code{matrix} or \code{data.frame}). Each q-vector must measure two attributes, reflecting the attributes measured by its statements.}

\item{polarity}{A \emph{F} blocks x 2 (\code{matrix} or \code{data.frame}). Each row reflects the direction of the first and second statement, where 1 and -1 corresponds to direct and inverse statements, respectively. Default is \code{NULL}, denoting that all statements are direct.}

\item{polarity.initial}{A \code{numeric} value that indicates the initial value for the estimation of the probability of endorsement for the latent group whose ideal response is equal to 0. The initial value for the latent group whose ideal response is equal to 1 will be 1 - \code{polarity.initial}. The initial value for latent groups without a clear ideal response is always equal to 0.5. This argument is ignored if \code{catprob.parm != NULL}. Default is \code{1e-4}.}

\item{att.dist}{How is the joint attribute distribution estimated? It can be \code{"saturated"}, \code{"higher.order"}, \code{"fixed"}, \code{"independent"}, and \code{"loglinear"}. Only considered if EM estimation is used. Default is \code{"saturated"}. See the \code{GDINA} package documentation for more information.}

\item{att.prior}{A \code{vector} of length 2^K to speficy attribute prior distribution for the latent classes. Only considered if EM estimation is used. Default is \code{NULL}. See the \code{GDINA} package documentation for more information.}

\item{verbose}{How to print calibration information after each EM iteration? Can be 0, 1 or 2, indicating to print no information, information for current iteration, or information for all iterations.}

\item{higher.order}{A \code{list} specifying the higher-order joint attribute distribution with the following components. Only considered if EM estimation is used. See the \code{GDINA} package documentation for more information.}

\item{catprob.parm}{A \code{list} of initial values for probabilities of endorsement for each nonzero category. Default is \code{NULL}. See the \code{GDINA} package documentation for more information.}

\item{control}{A \code{list} of control parameters. Only considered if EM estimation is used. See the \code{GDINA} package documentation for more information.}
}
\value{
\code{FCGDINA} returns an object of class \code{FCGDINA}.
\describe{
\item{\code{GDINA.obj}}{Estimation output from the \code{GDINA} function of the \code{GDINA.MJ} (Ma & Jiang, 2021) function, depending on whether EM or BM estimation has been used (\code{list}).}
\item{\code{technical}}{Information about initial values (\code{list}).}
\item{\code{specifications}}{Function call specifications (\code{list}).}
}
}
\description{
Estimation of the G-DINA model for forced-choice responses according to Nájera et al. (2024).
Block polarity (i.e., statement direction) and initial values for parameters can be specified to determine the design of the forced-choice blocks.
The \code{GDINA} package (Ma & de la Torre, 2020) is used to estimate the model via expectation maximumation (EM) algorithm if no priors are used.
To estimate the \emph{forced-choice diagnostic classification model} (FC-DCM; Huang, 2023) using Bayes modal estimation, please check the codes provided in https://osf.io/h6x9e/.
Only unidimensional statements (i.e., bidimensional blocks) are currently supported.
}
\examples{
\donttest{
library(GDINA)
set.seed(123)
# Q-matrix for the unidimensional statements
Q.items <- do.call("rbind", replicate(5, diag(5), simplify = FALSE))
# Guessing and slip
GS <- cbind(runif(n = nrow(Q.items), min = 0.1, max = 0.3),
            runif(n = nrow(Q.items), min = 0.1, max = 0.3))
n.blocks <- 30 # Number of forced-choice blocks

#----------------------------------------------------------------------------------------
# Illustration with simulated data using only direct statements (i.e., homopolar blocks)
#----------------------------------------------------------------------------------------

# Block polarity (1 = direct statement; -1 = indirect statement)
polarity <- matrix(1, nrow = n.blocks, ncol = 2)
sim <- simFCGDINA(N = 1000, Q.items, n.blocks = n.blocks, polarity = polarity,
                  model = "GDINA", GDINA.args = list(GS = GS), seed = 123)
Q <- sim$Q # Generated Q-matrix of forced-choice blocks
dat <- sim$dat # Generated responses
att <- sim$att # Generated attribute profiles

fit <- FCGDINA(dat = dat, Q = Q, polarity = polarity) # Fit the G-DINA model with EM estimation
ClassRate(personparm(fit$GDINA.obj), att) # Classification accuracy

#-------------------------------------------------------------------------------------------
# Illustration with simulated data using some inverse stataments (i.e., heteropolar blocks)
#-------------------------------------------------------------------------------------------

polarity <- matrix(1, nrow = n.blocks, ncol = 2)
# Including 15 inverse statements
polarity[sample(x = 1:(2*n.blocks), size = 15, replace = FALSE)] <- -1
sim <- simFCGDINA(N = 1000, Q.items, n.blocks = n.blocks, polarity = polarity,
                  model = "GDINA", GDINA.args = list(GS = GS), seed = 123)
Q <- sim$Q
dat <- sim$dat
att <- sim$att

fit <- FCGDINA(dat = dat, Q = Q, polarity = polarity)
ClassRate(personparm(fit$GDINA.obj), att)
}
}
\references{
Huang, H.-Y. (2023). Diagnostic Classification Model for Forced-Choice Items and Noncognitive Tests. \emph{Educational and Psychological Measurement}, \emph{83}(1), 146-180. https://doi.org/10.1177/00131644211069906

Ma, W., & de la Torre, J. (2020). GDINA: An R package for cognitive diagnosis modeling. \emph{Journal of Statistical Software}, \emph{93}(14). https://doi.org/10.18637/jss.v093.i14

Nájera, P., Kreitchmann, R. S., Escudero, S., Abad, F. J., de la Torre, J., & Sorrel, M. A. (2025). A General Diagnostic Modeling Framework for Forced-Choice Assessments. \emph{British Journal of Mathematical and Statistical Psychology}.
}
\author{
{Pablo Nájera, Universidad Pontificia Comillas}
}
