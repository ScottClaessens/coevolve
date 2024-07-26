library(ape)
library(readxl)
library(tidyverse)

# download xls file to temporary location
temp <-
  tempfile(fileext = ".xls")

file <-
  utils::download.file(
    paste0(
      "https://static-content.springer.com/esm/",
      "art%3A10.1038%2Fs41559-017-0112/MediaObjects/",
      "41559_2017_BFs415590170112_MOESM250_ESM.xls"
    ),
    destfile = temp,
    mode = "wb"
  )

# load xls sheet "Brain Final"
brain_data <-
  read_xls(
    path = temp,
    sheet = 3
  ) %>%
  # rename columns
  transmute(
    species = Taxon,
    brain_weight = `Final Brain Weight (g)`
  )

# load xls sheet "Body Data"
body_data <-
  read_xls(
    path = temp,
    sheet = 4
  ) %>%
  # rename columns
  transmute(
    species = Taxon,
    body_weight = `Final Body Weight (g)`
  )

# load xls sheet "Diet Data"
diet_data <-
  read_xls(
    path = temp,
    sheet = 5
  ) %>%
  # rename columns
  transmute(
    species = Taxon,
    diet = factor(
      c(
        "Fol"      = "Folivore",
        "Frug"     = "Frugivore",
        "Frug/Fol" = "Frugivore/Folivore",
        "Om"       = "Omnivore"
        )
      [`Diet Category`]
      ),
    percent_fruit = parse_number(`% Fruit`, na = "n/a")
  )

# load xls sheet "System Data"
system_data <-
  read_xls(
    path = temp,
    sheet = 6
  ) %>%
  # rename columns
  transmute(
    species = Taxon,
    social_system = factor(`Social System`),
    mating_system = factor(`Mating System`)
  )

# load xls sheet "Group Size Pivot"
group_size_data <-
  read_xls(
    path = temp,
    sheet = 8,
    skip = 5
  ) %>%
  # rename columns
  transmute(
    species = `Row Labels`,
    group_size = `Average of Group Size`
  )

# put datasets together
data <-
  brain_data %>%
  left_join(body_data, by = "species") %>%
  left_join(diet_data, by = "species") %>%
  left_join(system_data, by = "species") %>%
  left_join(group_size_data, by = "species") %>%
  # remove one species not present in phylogeny
  filter(species != "Hoolock_hoolock") %>%
  # as data frame
  as.data.frame()

# load consensus phylogeny of primates from 10k trees website
phylogeny <-
  read.nexus(
    file = paste0(
      "https://10ktrees.nunn-lab.org/Primates/downloads/version3/",
      "10kTrees_consensusTree_version3.nex"
    )
  )

# prune phylogeny to species in dataset
phylogeny <-
  keep.tip(
    phylogeny$`'con50majrule_noCladeCredibilityValues'`,
    data$species
  )

# put together
primates <- list(
  data = data,
  phylogeny = phylogeny
)

# save
usethis::use_data(primates, overwrite = TRUE)
