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
  cat(x$message)
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
