% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coev_calculate_theta.R
\name{coev_calculate_theta}
\alias{coev_calculate_theta}
\title{Calculate optimal trait values (theta) for a fitted \code{coevfit} object}
\usage{
coev_calculate_theta(object, intervention_values)
}
\arguments{
\item{object}{An object of class \code{coevfit}}

\item{intervention_values}{A named list of variables and associated
intervention values for calculating optimal trait values. All coevolving
variables must be declared separately in the named list without repetition.
If the intervention value for a particular variable is set to NA, this
variable is treated as a free variable. Otherwise, if the intervention
value for a particular variable is specified, the variable is held
constant at this trait value in the calculation. At least one variable must
be declared as a free variable and at least one variable must be held
constant (e.g., \code{list(var1 = NA, var2 = 0)}).}
}
\value{
Posterior samples in matrix format
}
\description{
Calculate optimal trait values (theta) for a fitted \code{coevfit} object
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

# calculate theta given intervention values
coev_calculate_theta(
  object = fit,
  intervention_values = list(
    political_authority = NA,
    religious_authority = 0
    )
  )
}
}
