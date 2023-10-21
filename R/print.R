#' Print the summary of a fitted model represented by a \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}
#' @param digits The number of significant digits for printing out the summary; defaults to 2
#'
#' @export
print.coevfit <- function(object, digits = 2) {
  # get summary of parameters from cmdstanr
  s <- as.data.frame(object$fit$summary(variables = c("alpha","sigma")))
  # print variables and response distributions
  cat("Variables: ")
  for (j in 1:length(object$variables)) {
    if (j == 1) {
      cat(names(object$variables)[j], "=", as.character(object$variables)[j], "\n")
    } else {
      cat("          ",
          names(object$variables)[j], "=", as.character(object$variables)[j], "\n")
    }
  }
  # print data name and number of observations
  cat(paste0("     Data: ", object$data_name, " (Number of observations: ", nrow(object$data), ")\n"))
  # print mcmc settings
  cat(paste0("    Draws: ", object$fit$num_chains(), " chains, each with iter = ",
             object$fit$metadata()$iter_sampling, "; warmup = ",
             object$fit$metadata()$iter_warmup, "; thin = ",
             object$fit$metadata()$thin, "\n           total post-warmup draws = ",
             object$fit$metadata()$iter_sampling * object$fit$num_chains(), "\n\n"))
  # print autoregressive effects
  cat("Autoregressive selection effects:\n")
  auto <- s[grepl("alpha", s$variable) & (substring(s$variable, 7, 7) == substring(s$variable, 9, 9)),]
  rownames(auto) <- names(object$variables)[as.integer(substring(auto$variable, 7, 7))]
  print(format(round(auto[,2:ncol(auto)], digits), nsmall = digits))
  cat("\n")
  # print cross effects
  cat("Cross selection effects:\n")
  cross <- s[grepl("alpha", s$variable) & (substring(s$variable, 7, 7) != substring(s$variable, 9, 9)),]
  rownames(cross) <-
    paste0(
      names(object$variables)[as.integer(substring(cross$variable, 9, 9))],
      " \U27F6 ",
      names(object$variables)[as.integer(substring(cross$variable, 7, 7))]
    )
  print(format(round(cross[,2:ncol(cross)], digits), nsmall = digits))
  cat("\n")
  # print drift
  cat("Drift scale parameters:\n")
  drift <- s[grepl("sigma", s$variable),]
  rownames(drift) <- names(object$variables)[as.integer(substring(drift$variable, 7, 7))]
  print(format(round(drift[,2:ncol(drift)], digits), nsmall = digits))
}
