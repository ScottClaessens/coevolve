#' Plot flowfield of expected evolutionary change from a fitted \code{coevfit}
#' object
#'
#' @param object An object of class \code{coevfit}
#' @param var1 A character string equal to one of the coevolving variables in
#'   the model
#' @param var2 A character string equal to one of the coevolving variables in
#'   the model
#' @param nullclines Logical (defaults to FALSE); whether to show coloured
#'   nullclines
#' to indicate where each variable is at equilibrium, depending on the state of
#'   the other
#'
#' @return A flowfield plot drawn directly to the device
#' @export
#'
#' @examples
#' \dontrun{
#' # fit dynamic coevolutionary model
#' fit <- coev_fit(
#'   data = authority$data,
#'   variables = list(
#'     political_authority = "ordered_logistic",
#'     religious_authority = "ordered_logistic"
#'   ),
#'   id = "language",
#'   tree = authority$phylogeny,
#'   # additional arguments for cmdstanr::sample()
#'   chains = 4,
#'   parallel_chains = 4,
#'   seed = 1
#'   )
#'
#' # plot flow field
#' coev_plot_flowfield(
#'   object = fit,
#'   var1 = "political_authority",
#'   var2 = "religious_authority"
#'   )
#' }
coev_plot_flowfield <- function(object, var1, var2, nullclines = FALSE) {
  # stop if object is not of class coevfit
  if (!methods::is(object, "coevfit")) {
    stop2("Argument 'object' must be a fitted coevolutionary model of class coevfit.")
  }
  if (!is.character(var1) | length(var1) != 1) {
    # stop if var1 not character string of length one
    stop2("Argument 'var1' must be a character string of length one.")
  } else if (!(var1 %in% names(object$variables))) {
    # stop if var1 not included in model
    stop2("Argument 'var1' must be a variable included in the fitted model.")
  }
  if (!is.character(var2) | length(var2) != 1) {
    # stop if var2 not character string of length one
    stop2("Argument 'var2' must be a character string of length one.")
  } else if (!(var2 %in% names(object$variables))) {
    # stop if var2 not included in model
    stop2("Argument 'var2' must be a variable included in the fitted model.")
  }
  # stop if var1 and var2 are the same variable
  if (var1 == var2) {
    stop2("Argument 'var1' and 'var2' must refer to different variables.")
  }
  # stop if nullclines not logical
  if (!is.logical(nullclines)) {
    stop2("Argument 'nullclines' must be logical.")
  }
  # get IDs for variables
  id_var1 <- which(names(object$variables) == var1)
  id_var2 <- which(names(object$variables) == var2)
  # get posterior draws
  draws <- posterior::as_draws_rvars(object$fit)
  # medians and median absolute deviations for all variables
  eta  <- apply(
    draws$eta[,1:object$stan_data$N_tips,], 3, posterior::rvar_median
    )
  meds <- unlist(lapply(eta, stats::median))
  mads <- unlist(lapply(eta, stats::mad))
  lowers <- meds - 2.5*mads
  uppers <- meds + 2.5*mads
  # get median parameter values for A and b
  A <- stats::median(draws$A)
  b <- stats::median(draws$b)
  # function for flow field diagram
  OU <- function(t, y, parameters) {
    dy <- numeric(2)
    # variable 1
    dy[1] <- b[id_var1]
    for (j in 1:length(names(object$variables))) {
      if (j == id_var1) {
        # autoregressive effect
        dy[1] <- dy[1] + A[id_var1,j]*y[1]
      } else if (j == id_var2) {
        # cross-lagged effect of predictor
        dy[1] <- dy[1] + A[id_var1,j]*y[2]
      } else {
        # cross-lagged effects of other variables held at their median values
        dy[1] <- dy[1] + A[id_var1,j]*meds[j]
      }
    }
    # variable 2
    dy[2] <- b[id_var2]
    for (j in 1:length(names(object$variables))) {
      if (j == id_var2) {
        # autoregressive effect
        dy[2] <- dy[2] + A[id_var2,j]*y[2]
      } else if (j == id_var1) {
        # cross-lagged effect of predictor
        dy[2] <- dy[2] + A[id_var2,j]*y[1]
      } else {
        # cross-lagged effects of other variables held at their median values
        dy[2] <- dy[2] + A[id_var2,j]*meds[j]
      }
    }
    return(list(dy))
  }
  # create flow field diagram
  suppressWarnings({
    OU.flowField <-
      phaseR::flowField(
        OU,
        xlim = c(lowers[id_var1], uppers[id_var1]),
        ylim = c(lowers[id_var2], uppers[id_var2]),
        parameters = NA,
        add = FALSE,
        xlab = "",
        ylab = "",
        points = 12,
        col = "grey",
        xaxt = 'n',
        yaxt = 'n',
        arrow.type = "proportional",
        frac = 1.5,
        xaxs = "i",
        yaxs = "i",
        axes = FALSE,
        lwd = 2
      )
  })
  # var 1 label
  graphics::mtext(
    side = 1,
    text = paste0(var1, " (z-score)"),
    at = meds[id_var1],
    line = 2.5,
    cex = 1.3
  )
  # var 2 label
  graphics::mtext(
    side = 2,
    text = paste0(var2, " (z-score)"),
    at = meds[id_var2],
    line = 2.5,
    cex = 1.3
  )
  # add nullclines to phase plane
  suppressWarnings({
    if (nullclines) {
      nc <-
        phaseR::nullclines(
          OU,
          xlim = c(lowers[id_var1], uppers[id_var1]),
          ylim = c(lowers[id_var2], uppers[id_var2]),
          parameters = NA,
          points = 20,
          axes = FALSE,
          col = c("#c55852","#5387b6"),
          add.legend = FALSE,
          lwd = 2
        )
    }
  })
  # add axes
  graphics::axis(
    side = 1,
    at = c(lowers[id_var1], meds[id_var1], uppers[id_var1]),
    labels =
      (c(lowers[id_var1], meds[id_var1], uppers[id_var1]) - meds[id_var1]) /
      mads[id_var1]
  )
  graphics::axis(
    side = 2,
    at = c(lowers[id_var2], meds[id_var2], uppers[id_var2]),
    labels =
      (c(lowers[id_var2], meds[id_var2], uppers[id_var2]) - meds[id_var2]) /
      mads[id_var2]
  )
}
