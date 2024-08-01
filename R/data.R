#' Political and religious authority in Austronesian societies
#'
#' @description This dataset and associated phylogeny describes the states of
#'   political and religious authority in 97 Austronesian societies. These data
#'   were originally compiled in Sheehan et al. (2023). Political authority is
#'   defined as a right to manage interactions between living human beings,
#'   whereas religious authority is defined as a right to manage interactions
#'   between living human beings and supernatural agents or powers. Authority
#'   is coded as being absent, sublocal (smaller than the local community),
#'   local (coextensive with the local community), or supralocal (consisting of
#'   more than one local community).
#'
#' @format A list containing a dataset and an associated language phylogeny. The
#'   phylogeny is a pruned maximum clade credibility tree for 97 Austronesian
#'   languages. The dataset is a data frame with 97 observations and three
#'   variables:
#' \describe{
#'  \item{language}{The name of the Austronesian language linked to each
#'    society}
#'  \item{political_authority}{An ordered factor indicating whether political
#'    authority is absent, sublocal, local, or supralocal in each society}
#'  \item{religious_authority}{An ordered factor indicating whether religious
#'    authority is absent, sublocal, local, or supralocal in each society}
#' }
#'
#' @source Sheehan, O., Watts, J., Gray, R. D., Bulbulia, J., Claessens, S.,
#'   Ringen, E. J., & Atkinson, Q. D. (2023). Coevolution of religious and
#'   political authority in Austronesian societies. \emph{Nature Human
#'   Behaviour}, \emph{7}(1), 38-45.
#'
#' @examples
#' \dontrun{
#' # fit model to authority data
#' m <-
#'   coev_fit(
#'     data = authority$data,
#'     variables = list(
#'       political_authority = "ordered_logistic",
#'       religious_authority = "ordered_logistic"
#'     ),
#'     id = "language",
#'     tree = authority$phylogeny,
#'     # arguments for cmdstanr::sample()
#'     parallel_chains = 4,
#'     seed = 1,
#'     # set prior manually
#'     prior = list(A_offdiag = "normal(0, 2)")
#'   )
#' # print model summary
#' summary(m)
#' # plot delta theta
#' coev_plot_delta_theta(m)
#' }
#'
"authority"

#' Brain weight, body weight, diet, and sociality in primates
#'
#' @description This dataset and associated phylogeny describes brain weights,
#'   body weights, dietary categories, social systems, mating systems, and group
#'   sizes for 143 primate species. These data were originally compiled in
#'   DeCasien et al. (2017). The phylogeny is a 50% majority rule consensus
#'   tree from the 10k Trees website.
#'
#' @format A list containing a dataset and an associated phylogeny. The
#'   phylogeny is a pruned consensus tree for 143 species of primate (excluding
#'   humans). The dataset is a data frame with 143 observations and 11
#'   variables:
#' \describe{
#'  \item{species}{The name of the primate species}
#'  \item{brain_weight}{Average brain weight, in grams}
#'  \item{body_weight}{Average body weight, in grams}
#'  \item{diet}{A factor with four levels indicating whether the species'
#'    dietary category is folivore, frugivore, frugivore/folivore, or omnivore}
#'  \item{percent_fruit}{Percentage of fruit in diet}
#'  \item{social_system}{A factor with four levels indicating whether the
#'    species' social system is solitary, pair-living, harem polygyny, or
#'    polygynandry}
#'  \item{mating_system}{A factor with five levels indicating whether the
#'    species' mating system is spatial polygyny, monogamy, polyandry, harem
#'    polygyny, or polygynandry}
#'  \item{group_size}{Average group size}
#' }
#'
#' @source DeCasien, A. R., Williams, S. A., & Higham, J. P. (2017). Primate
#'   brain size is predicted by diet but not sociality. \emph{Nature Ecology &
#'   Evolution}, \emph{1}(5), 0112.
#' @source https://10ktrees.nunn-lab.org/Primates/index.html
#'
#' @examples
#' \dontrun{
#' # fit model to primates data
#' m <-
#'   coev_fit(
#'     data = primates$data,
#'     variables = list(
#'       brain_weight = "lognormal",
#'       body_weight = "lognormal"
#'     ),
#'     id = "species",
#'     tree = primates$phylogeny,
#'     # arguments to cmdstanr::sample()
#'     parallel_chains = 4,
#'     seed = 1
#'   )
#' # print model summary
#' summary(m)
#' # plot delta theta
#' coev_plot_delta_theta(m)
#' }
#'
"primates"

#' Example dataset with repeated observations
#'
#' @description This example dataset and associated phylogeny are used as an
#'   example of repeated observations in phylogenetic modelling. The data are
#'   adapted from de Villemeruil & Nakagawa (2014) to include only the first
#'   twenty species.
#'
#' @format A list containing a dataset and an associated phylogeny. The dataset
#'   is a data frame with 100 observations and 3 variables:
#' \describe{
#'  \item{species}{The name of the species}
#'  \item{y}{An example continuous variable}
#'  \item{x}{An example continuous variable}
#' }
#'
#' @source de Villemeruil P. & Nakagawa, S. (2014). General quantitative genetic
#'   methods for comparative biology. In L. Garamszegi (Ed.), *Modern
#'   phylogenetic comparative methods and their application in evolutionary
#'   biology: concepts and practice* (pp. 287-303). Springer, New York.
#'
#' @examples
#' \dontrun{
#' # fit model to primates data
#' m <-
#'   coev_fit(
#'     data = repeated$data,
#'     variables = list(
#'       y = "normal",
#'       x = "normal"
#'     ),
#'     id = "species",
#'     tree = repeated$phylogeny,
#'     # arguments to cmdstanr::sample()
#'     parallel_chains = 4,
#'     seed = 1
#'   )
#' # print model summary
#' summary(m)
#' # plot delta theta
#' coev_plot_delta_theta(m)
#' }
#'
"repeated"
