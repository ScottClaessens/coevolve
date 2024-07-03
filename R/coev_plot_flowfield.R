#' Plot flowfield of expected evolutionary change from a fitted \code{coevfit} object
#'
#' @param object An object of class \code{coevfit}
#' @param var1 A character string equal to one of the coevolving variables in the model
#' @param var2 A character string equal to one of the coevolving variables in the model
#' @param nullclines Logical (defaults to FALSE); whether to show coloured nullclines
#' to indicate where each variable is at equilibrium, depending on the state of the other
#'
#' @return A flowfield plot drawn directly to the device
#' @export
#'
#' @examples
#' \dontrun{
#' # simulate data
#' n <- 20
#' tree <- ape::rcoal(n)
#' d <- data.frame(
#'   id = tree$tip.label,
#'   x = rbinom(n, size = 1, prob = 0.5),
#'   y = ordered(sample(1:4, size = n, replace = TRUE))
#' )
#' # fit dynamic coevolutionary model
#' m <- coev_fit(
#'   data = d,
#'   variables = list(
#'     x = "bernoulli_logit",
#'     y = "ordered_logistic"
#'   ),
#'   id = "id",
#'   tree = tree,
#'   # additional arguments for cmdstanr::sample()
#'   chains = 4,
#'   parallel_chains = 4,
#'   iter_warmup = 500,
#'   seed = 1
#' )
#' # plot flowfield
#' coev_plot_flowfield(m)
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
  # medians and median absolute deviations for both variables
  eta_var1 <- posterior::rvar_median(draws$eta[1:nrow(object$data), id_var1])
  eta_var2 <- posterior::rvar_median(draws$eta[1:nrow(object$data), id_var2])
  med_var1 <- stats::median(eta_var1)
  med_var2 <- stats::median(eta_var2)
  mad_var1 <- stats::mad(eta_var1)
  mad_var2 <- stats::mad(eta_var2)
  lower_var1 <- med_var1 - 2.5*mad_var1
  lower_var2 <- med_var2 - 2.5*mad_var2
  upper_var1 <- med_var1 + 2.5*mad_var1
  upper_var2 <- med_var2 + 2.5*mad_var2
  # get median parameter values for A and b
  A <- stats::median(draws$A)
  b <- stats::median(draws$b)
  # function for flow field diagram
  OU <- function(t, y, parameters) {
    dy <- numeric(2)
    dy[1] <- y[1]*A[id_var1,id_var1] + y[2]*A[id_var1,id_var2] + b[id_var1]
    dy[2] <- y[2]*A[id_var2,id_var2] + y[1]*A[id_var2,id_var1] + b[id_var2]
    list(dy)
  }
  # create flow field diagram
  suppressWarnings({
    OU.flowField <-
      phaseR::flowField(
        OU,
        xlim = c(lower_var1, upper_var1),
        ylim = c(lower_var2, upper_var2),
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
    at = med_var1,
    line = 2.5,
    cex = 1.3
  )
  # var 2 label
  graphics::mtext(
    side = 2,
    text = paste0(var2, " (z-score)"),
    at = med_var2,
    line = 2.5,
    cex = 1.3
  )
  # add nullclines to phase plane
  suppressWarnings({
    if (nullclines) {
      nc <-
        phaseR::nullclines(
          OU,
          xlim = c(lower_var1, upper_var1),
          ylim = c(lower_var2, upper_var2),
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
    at = c(lower_var1, med_var1, upper_var1),
    labels = (c(lower_var1, med_var1, upper_var1) - med_var1) / mad_var1
  )
  graphics::axis(
    side = 2,
    at = c(lower_var2, med_var2, upper_var2),
    labels = (c(lower_var2, med_var2, upper_var2) - med_var2) / mad_var2
  )
}
