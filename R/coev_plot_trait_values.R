#' Plot estimated trait values from a fitted \code{coevfit} object
#'
#' Produce a pairs plot of the estimated trait values for taxa from a fitted
#' \code{coevfit} object. The plot includes scatterplot(s) of median estimated
#' trait values, heatmap(s) of median estimated trait values, and density plots
#' of the marginal distributions for each variable with associated posterior
#' uncertainty.
#'
#' @srrstats {G1.3, G1.4, G2.1a} Function documentation begins here, with
#'   expected data types and definitions of statistical terminology and inputs
#'
#' @param object An object of class \code{coevfit}
#' @param variables If NULL (default), the function returns a pairs plot
#'   including all coevolving variables from the model. Otherwise, a character
#'   vector declaring the variables to be included.
#' @param ndraws An integer indicating the number of draws to return in the
#'   density plots on the diagonal. The default and maximum number of draws is
#'   the size of the posterior sample.
#' @param tree_id An integer indicating the tree ID to use when making
#'   posterior predictions. Set to \code{NULL} by default, which will use draws
#'   from every tree, integrating phylogenetic uncertainty.
#' @param xlim Limits for the x-axis. If \code{NULL} (default), limits are set
#'   to the minimum and maximum estimated trait values.
#' @param ylim Limits for the y-axis. If \code{NULL} (default), limits are set
#'   to the minimum and maximum estimated trait values.
#'
#' @returns A patchwork of \code{ggplot} objects
#'
#' @author Scott Claessens \email{scott.claessens@@gmail.com}, Erik Ringen
#'   \email{erikjacob.ringen@@uzh.ch}
#'
#' @references
#' Ringen, E., Martin, J. S., & Jaeggi, A. (2021). Novel phylogenetic methods
#' reveal that resource-use intensification drives the evolution of "complex"
#' societies. \emph{EcoEvoRXiv}. \code{doi:10.32942/osf.io/wfp95}
#'
#' Sheehan, O., Watts, J., Gray, R. D., Bulbulia, J., Claessens, S., Ringen,
#' E. J., & Atkinson, Q. D. (2023). Coevolution of religious and political
#' authority in Austronesian societies. \emph{Nature Human Behaviour},
#' \emph{7}(1), 38-45. \code{10.1038/s41562-022-01471-y}
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
#' # pairs plot of trait values
#' coev_plot_trait_values(fit)
#' }
#'
#' @export
coev_plot_trait_values <- function(object, variables = NULL, ndraws = 50,
                                   tree_id = NULL, xlim = NULL, ylim = NULL) {
  # run checks
  run_checks_plot_trait_values(object, variables, ndraws, tree_id, xlim, ylim)
  # get posterior trait values
  # Handle both cmdstanr and nutpie
  if (inherits(object$fit, "nutpie_fit")) {
    # For nutpie, extract draws first
    draws_obj <- object$fit$draws()
    draws_rvars <- posterior::as_draws_rvars(draws_obj)
    eta <- draws_rvars$eta
  } else {
    # For cmdstanr, use as_draws_rvars method
    draws_rvars <- posterior::as_draws_rvars(object$fit)
    eta <- draws_rvars$eta
  }
  # initial rvars dimensions: [trees, nodes, variables]
  # if particular tree defined, subset to that tree
  if (!is.null(tree_id)) eta <- eta[tree_id, , ]
  # average over tree(s)
  eta <- apply(
    eta, c(2, 3), function(x) {
      posterior::rvar(as.vector(posterior::draws_of(x)))
    }
  )
  # new rvars dimensions: [nodes, variables]
  # subset to taxa (remove internal nodes)
  eta <- eta[1:object$stan_data$N_tips, ]
  # new rvars dimensions: [tips, variables]
  # create dataset for plotting
  d <- tibble::tibble(tip = seq_len(nrow(eta)))
  for (j in seq_along(object$variables)) {
    d[[names(object$variables)[j]]] <- eta[, j]
  }
  # if particular variables defined, keep only those
  if (!is.null(variables)) {
    d <- d[, variables]
  } else {
    d <- d[, -1] # else just remove tip column
  }
  # get median trait values
  d_med <-
    d |>
    dplyr::rowwise() |>
    dplyr::mutate(
      dplyr::across(dplyr::where(is.list), function(x) stats::median(x[[1]]))
    ) |>
    tidyr::unnest(dplyr::everything())
  # get model summary
  s <- summary(object)
  # get draws
  all_draw_ids <- 1:(s$iter * s$chains * ifelse(is.null(tree_id), s$ntrees, 1))
  if (!is.null(ndraws)) {
    draw_ids <- sample(all_draw_ids, size = ndraws)
  } else {
    draw_ids <- all_draw_ids
  }
  # function for plotting in upper triangle
  plot_upper_tri <- function(i, j) {
    ggplot2::ggplot(
      data = d_med,
      mapping = ggplot2::aes(
        x = !!dplyr::sym(colnames(d)[j]),
        y = !!dplyr::sym(colnames(d)[i])
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::scale_x_continuous(
        expand = c(0, 0),
        limits = c(
          ifelse(is.null(xlim), min(d_med), xlim[1]),
          ifelse(is.null(xlim), max(d_med), xlim[2])
        )
      ) +
      ggplot2::scale_y_continuous(
        expand = c(0, 0),
        limits = c(
          ifelse(is.null(ylim), min(d_med), ylim[1]),
          ifelse(is.null(ylim), max(d_med), ylim[2])
        )
      ) +
      ggplot2::theme_bw()
  }
  # function for plotting in lower triangle
  plot_lower_tri <- function(i, j) {
    ggplot2::ggplot(
      data = d_med,
      mapping = ggplot2::aes(
        x = !!dplyr::sym(colnames(d)[j]),
        y = !!dplyr::sym(colnames(d)[i])
      )
    ) +
      ggplot2::stat_density_2d(
        geom = "polygon",
        contour = TRUE,
        ggplot2::aes(fill = ggplot2::after_stat(.data$level)),
        bins = 10
      ) +
      ggplot2::scale_fill_distiller(
        palette = "Reds",
        direction = 1
      ) +
      ggplot2::scale_x_continuous(
        expand = c(0, 0),
        limits = c(
          ifelse(is.null(xlim), min(d_med), xlim[1]),
          ifelse(is.null(xlim), max(d_med), xlim[2])
        )
      ) +
      ggplot2::scale_y_continuous(
        expand = c(0, 0),
        limits = c(
          ifelse(is.null(ylim), min(d_med), ylim[1]),
          ifelse(is.null(ylim), max(d_med), ylim[2])
        )
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")
  }
  # function for plotting on diagonal
  plot_diag <- function(i) {
    d |>
      tidyr::unnest(dplyr::everything()) |>
      tibble::rownames_to_column("taxa") |>
      dplyr::rowwise() |>
      dplyr::mutate(
        dplyr::across(
          !.data$taxa, function(x) {
            list(posterior::draws_of(x)[, 1])
          }
        )
      ) |>
      tidyr::unnest(!.data$taxa) |>
      dplyr::mutate(
        iter = rep(all_draw_ids, times = max(as.numeric(.data$taxa)))
      ) |>
      dplyr::filter(.data$iter %in% draw_ids) |>
      ggplot2::ggplot() +
      ggplot2::geom_density(
        ggplot2::aes(
          x = !!dplyr::sym(colnames(d)[i]),
          group = .data$iter
        ),
        linewidth = 0.1,
        colour = "lightblue"
      ) +
      ggplot2::geom_density(
        data = d_med,
        ggplot2::aes(x = !!dplyr::sym(colnames(d)[i])),
        colour = "black"
      ) +
      ggplot2::scale_x_continuous(
        expand = c(0, 0),
        limits = c(
          ifelse(is.null(xlim), min(d_med), xlim[1]),
          ifelse(is.null(xlim), max(d_med), xlim[2])
        )
      ) +
      ggplot2::scale_y_continuous(
        name = colnames(d)[i],
        expand = c(0, 0)
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  }
  # list of plots
  plot_list <- as.list(rep(NA, times = ncol(d) * ncol(d)))
  # loop over plots
  count <- 0
  for (i in seq_len(ncol(d))) {
    for (j in seq_len(ncol(d))) {
      # increment count
      count <- count + 1
      # diagonal plots
      if (i == j) {
        p <- plot_diag(i)
        if (i != 1) {
          p <- p + ggplot2::theme(axis.title.y = ggplot2::element_blank())
        }
        if (i != ncol(d)) {
          p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                  axis.text.x = ggplot2::element_blank())
        }
        plot_list[[count]] <- p
      } else if (i < j) {
        # upper triangle plots
        p <- plot_upper_tri(i, j)
        if (i != ncol(d)) {
          p <- p + ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                                  axis.text.y = ggplot2::element_blank())
        }
        if (j != 1) {
          p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                  axis.text.x = ggplot2::element_blank())
        }
        plot_list[[count]] <- p
      } else if (i > j) {
        # lower triangle plots
        p <- plot_lower_tri(i, j)
        if (i != ncol(d)) {
          p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                  axis.text.x = ggplot2::element_blank())
        }
        if (j != 1) {
          p <- p + ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                                  axis.text.y = ggplot2::element_blank())
        }
        plot_list[[count]] <- p
      }
    }
  }
  # plot grid
  patchwork::wrap_plots(
    plot_list,
    nrow = ncol(d),
    ncol = ncol(d)
  )
}

#' Internal helper function for checking coev_plot_trait_values() arguments
#'
#' @srrstats {G1.4a} Non-exported function documented here
#'
#' @description Checks arguments for coev_plot_trait_values()
#'
#' @returns Error message if any of the checks fail
#'
#' @noRd
run_checks_plot_trait_values <- function(object, variables, ndraws, tree_id,
                                         xlim, ylim) {
  #' @srrstats {G5.2, G5.2a} Unique error messages for each input
  # stop if object is not of class coevfit
  #' @srrstats {G2.1} Assertion on type of input
  if (!methods::is(object, "coevfit")) {
    stop2(
      paste0(
        "Argument 'object' must be a fitted coevolutionary model ",
        "of class coevfit."
      )
    )
  }
  if (!is.null(variables)) {
    if (!is.character(variables) || !(length(variables) >= 2)) {
      # stop if variables is not a character string
      #' @srrstats {G2.0, G2.1} Assertion on length and type of input
      stop2(
        paste0(
          "Argument 'variables' must be a character vector ",
          "of at least length 2."
        )
      )
    } else if (!all(variables %in% names(object$variables))) {
      # stop if variables not included in the fitted model
      stop2(
        paste0(
          "Argument 'variables' contains variable names that are not ",
          "included in the fitted model."
        )
      )
    }
  }
  # stop if ndraws is not a single integer between 1 and the total num draws
  if (!is.null(ndraws)) {
    if (!is.numeric(ndraws)) {
      #' @srrstats {G2.1} Assertion on type of input
      stop2("Argument 'ndraws' must be numeric.")
    } else if (!all(as.integer(ndraws) == ndraws) || length(ndraws) != 1) {
      #' @srrstats {G2.0, G2.1, G2.2, G2.4, G2.4a} Assertion on length and type
      #' of input, convert to integer
      stop2("Argument 'ndraws' must be a single integer.")
    } else {
      # Get number of draws from stored nsamples in coevfit object
      n_draws_total <- object$nsamples
      if (ndraws < 1 || ndraws > n_draws_total) {
        stop2(
          "Argument 'ndraws' must be between 1 and the total number of draws."
        )
      }
    }
  }
  # stop if tree_id is not a single integer between 1 and the total num trees
  if (!is.null(tree_id)) {
    if (!is.numeric(tree_id)) {
      #' @srrstats {G2.1} Assertion on type of input
      stop2("Argument 'tree_id' must be numeric.")
    } else if (!all(as.integer(tree_id) == tree_id) || length(tree_id) != 1) {
      #' @srrstats {G2.0, G2.1, G2.2, G2.4, G2.4a} Assertion on length and type
      #' of input, convert to integer
      stop2("Argument 'tree_id' must be a single integer.")
    } else if (tree_id < 1 || tree_id > object$stan_data$N_tree) {
      stop2(
        "Argument 'tree_id' must be between 1 and the total number of trees."
      )
    }
  }
  # stop if xlim is not a numeric vector of length 2
  #' @srrstats {G2.0, G2.1} Assertion on length and type of input
  if (!is.null(xlim)) {
    if (!(is.numeric(xlim) && is.vector(xlim) && length(xlim) == 2)) {
      stop2("Argument 'xlim' must be a numeric vector of length 2.")
    }
  }
  # stop if ylim is not a numeric vector of length 2
  #' @srrstats {G2.0, G2.1} Assertion on length and type of input
  if (!is.null(ylim)) {
    if (!(is.numeric(ylim) && is.vector(ylim) && length(ylim) == 2)) {
      stop2("Argument 'ylim' must be a numeric vector of length 2.")
    }
  }
}
