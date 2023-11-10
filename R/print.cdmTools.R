#' @export
print.estQ <- function(x, ...){
  cat(paste0("Estimated Q-matrix using the DFL method using ", x$specifications$criterion, " criterion:", "\n", "\n"))
  print(as.data.frame(x$est.Q))
}
#' @export
print.genQ <- function(x, ...){
  cat("Generated Q-matrix:\n\n")
  print(as.data.frame(x$gen.Q))
}
#' @export
print.GNPC <- function(x, ...){
  conv <- ifelse(x$n.ite < x$specifications$maxitr, TRUE, FALSE)
  cat("GNPC method\n\n")
  cat(paste0("Number of examinees:  ", nrow(x$specifications$dat), "\n"))
  cat(paste0("Number of items:      ", nrow(x$specifications$Q), "\n"))
  cat(paste0("Number of attributes: ", ncol(x$specifications$Q), "\n", "\n"))
  cat(paste0("Number of iterations: ", x$n.ite, "\n"))
  cat(paste0("Max. change in alpha: ", round(max(x$hist.change), 3), "\n"))
  cat(paste0("Mean change in alpha: ", round(mean(x$hist.change), 3), "\n"))
  cat(paste0("Last change in alpha: ", round(x$hist.change[x$n.ite], 3), "\n"))
  cat(paste0("Convergence: ", conv, "\n"))
}
#' @export
print.is.Qid <- function(x, ...){
  cat("Model =", x$specifications$model, "\n\n")
  if(length(x$conditions) == 3){
    cat("Identifiability conditions:", "\n")
    cat("A) Completeness         =", x$conditions$completeness, "\n")
    cat("B) Distinctiveness      =", x$conditions$distinctiveness, "\n")
    cat("C) Repetition           =", x$conditions$repetition, "\n\n")
  } else if(length(x$conditions) == 5){
    cat("Identifiability conditions:", "\n")
    cat("A) Completeness         =", x$conditions$completeness, "\n")
    cat("B) Distinctiveness      =", x$conditions$distinctiveness, "\n")
    cat("C) Repetition           =", x$conditions$repetition, "\n")
    cat("D) Generic completeness =", x$conditions$generic.completeness, "\n")
    cat("E) Generic repetition   =", x$conditions$generic.repetition, "\n\n")
  }
  cat("Strict identifiabilty   =", x$strict, "\n")
  cat("Generic identifiabilty  =", x$generic, "\n")
}
#' @export
print.missQ <- function(x, ...){
  mQ <- as.data.frame(x$miss.Q)
  tQ <- as.data.frame(x$Q)
  mQ[mQ != tQ] <- paste0(mQ[mQ != tQ], "*")
  cat("Misspecified Q-matrix:\n\n")
  print(mQ, right = FALSE)
  cat("Note: * denotes a misspecified element.\n")
}
#' @export
print.modelcompK <- function(x, ...){
  if(is.null(x$specifications$Qs)){
    txt <- paste0("Suggested number of attributes based on model fit:")
    for(l in 1:length(table(x$sug.K))){
      txt <- paste0(txt, "\n", names(table(x$sug.K))[l], " attributes: ", paste(names(x$sug.K)[x$sug.K == as.numeric(names(table(x$sug.K))[l])], collapse = ", "))
    }
  } else {
    txt <- paste0("Selected Q-matrix based on model fit:")
    for(l in 1:length(table(x$sel.Q))){
      txt <- paste0(txt, "\n", names(table(x$sel.Q))[l], ": ", paste(names(x$sel.Q)[x$sel.Q == names(table(x$sel.Q))[l]], collapse = ", "))
    }
  }
  cat(txt)
}
#' @export
print.orderQ <- function(x, ...){
  config.MAD <- as.vector(x$configs[x$configs[,2] == min(x$configs[,2]), 1])
  min.MAD <- min(x$configs[,2])
  config.CC <- as.vector(x$configs[x$configs[,3] == max(x$configs[,3]), 1])
  max.CC <- max(x$configs[,3])
  cat(paste0("Best attribute ordering(s) according to MAD: ", paste(config.MAD, collapse = ", "), " (MAD = ", min.MAD, ")", "\n",
             "Best attribute ordering(s) according to CC:  ", paste(config.CC, collapse = ", "), " (CC = ", max.CC, ")"))
}
#' @export
print.paK <- function(x, ...){
  txt <- paste0("Suggested number of attributes based on parallel analysis:")
  for(l in 1:length(table(x$sug.K))){
    txt <- paste0(txt, "\n", names(table(x$sug.K))[l], " attributes: ", paste(names(x$sug.K)[x$sug.K == as.numeric(names(table(x$sug.K))[l])], collapse = ", "))
  }
  cat(txt)
  print(x$plot)
}
#' @export
print.valQ <- function(x, ...){
  vQ <- as.data.frame(x$sug.Q)
  oQ <- as.data.frame(x$Q)
  vQ[vQ != oQ] <- paste0(vQ[vQ != oQ], "*")
  implement <- switch(x$specifications$iterative,
                      "none" = "non-iterative implementation",
                      "test" = "test-level iterative implementation",
                      "test.att" = "attribute-test-level iterative implementation",
                      "item" = "item-level iterative implementation")
  cat(paste0("Q-matrix validation based on the Hull method", "\n", "(using ", x$specifications$index, " index and ", implement, ")", "\n", "\n"))
  cat("Suggested Q-matrix:\n\n")
  print(vQ, right = FALSE)
  cat("Note: * denotes a modified element.\n")
}
#' @export
print.RDINA <- function(x, ...){
  if(x$specifications$gate %in% c("AND", "OR")){
    Q <- x$specifications$Q
    pckg <- paste0("cdmTools ", utils::packageDescription("cdmTools")$Version,
                   ", NPCD ", utils::packageDescription("NPCD")$Version,
                   ", GDINA ", utils::packageDescription("GDINA")$Version)
  } else {
    Q <- x$specifications$Q
    pckg <- paste0("cdmTools ", utils::packageDescription("cdmTools")$Version,
                   ", GDINA ", utils::packageDescription("GDINA")$Version)
  }
  cat.phi.true <- ifelse(is.null(x$phi.true), "", paste0(" (True \U1D711 = ", sprintf("%.*f", 4, x$phi.true), ")"))
  cat(
    paste0(
      "============================================================", "\n",
      "                         R-DINA model                       ", "\n", "\n",
      "Packages used: ", pckg, "\n", "\n",
      "Model estimated using the '", x$specifications$gate, "' gate and ", x$specifications$est, " optimization,", "\n",
      nrow(x$specifications$dat), " individuals, ", nrow(Q), " items, and ", ncol(Q), " attributes.", "\n", "\n",
      "Estimated \U1D711 = ", sprintf("%.*f", 4, x$phi), cat.phi.true, "\n", "\n",
      "Relative model fit:", "\n",
      "  -2LL = ", sprintf("%.*f", 2, x$test.fit$Deviance), " | npar  = ", x$test.fit$npar, "\n",
      "  AIC  = ", sprintf("%.*f", 2, x$test.fit$AIC), " | BIC   = ", sprintf("%.*f", 2, x$test.fit$BIC), "\n",
      "  CAIC = ", sprintf("%.*f", 2, x$test.fit$CAIC), " | SABIC = ", sprintf("%.*f", 2, x$test.fit$SABIC), "\n", "\n",
      "Estimated classification accuracy (\U1D70F):", "\n",
      "  Test level      = ",  sprintf("%.*f", 3, x$class.accu$tau), "\n",
      "  Profile level   = from ", sprintf("%.*f", 3, min(x$class.accu$tau_l)), " {", names(which.min(x$class.accu$tau_l)), "} to ", sprintf("%.*f", 3, max(x$class.accu$tau_l)), " {", names(which.max(x$class.accu$tau_l)), "}", "\n",
      "  Attribute level = from ", sprintf("%.*f", 3, min(x$class.accu$tau_k)), " (", names(which.min(x$class.accu$tau_k)), ") to ", sprintf("%.*f", 3, max(x$class.accu$tau_k)), " (", names(which.max(x$class.accu$tau_k)), ")", "\n",
      "============================================================="
    )
  )
}
