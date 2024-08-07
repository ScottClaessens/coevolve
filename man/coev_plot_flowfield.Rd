% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coev_plot_flowfield.R
\name{coev_plot_flowfield}
\alias{coev_plot_flowfield}
\title{Plot flowfield of expected evolutionary change from a fitted \code{coevfit}
object}
\usage{
coev_plot_flowfield(object, var1, var2, nullclines = FALSE)
}
\arguments{
\item{object}{An object of class \code{coevfit}}

\item{var1}{A character string equal to one of the coevolving variables in
the model}

\item{var2}{A character string equal to one of the coevolving variables in
the model}

\item{nullclines}{Logical (defaults to FALSE); whether to show coloured
nullclines
to indicate where each variable is at equilibrium, depending on the state of
the other}
}
\value{
A flowfield plot drawn directly to the device
}
\description{
Plot flowfield of expected evolutionary change from a fitted \code{coevfit}
object
}
\examples{
\dontrun{
# fit dynamic coevolutionary model
fit <- coev_fit(
  data = authority$data,
  variables = list(
    political_authority = "ordered_logistic",
    religious_authority = "ordered_logistic"
  ),
  id = "language",
  tree = authority$phylogeny,
  # additional arguments for cmdstanr::sample()
  chains = 4,
  parallel_chains = 4,
  seed = 1
  )

# plot flow field
coev_plot_flowfield(
  object = fit,
  var1 = "political_authority",
  var2 = "religious_authority"
  )
}
}
