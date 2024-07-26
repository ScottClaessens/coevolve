library(tidyverse)

# load dataset for political authority from github
political_authority <-
  read_tsv(
    file = paste0(
      "https://raw.githubusercontent.com/ScottClaessens/phyloAuthority/main/",
      "data/political_authority_20_09_21.txt"
      ),
    col_types = "cd"
  )

# load dataset for religious authority from github
religious_authority <-
  read_tsv(
    file = paste0(
      "https://raw.githubusercontent.com/ScottClaessens/phyloAuthority/main/",
      "data/religious_authority_20_09_21.txt"
    ),
    col_types = "cd"
  )

# get ordered authority levels
authority_levels <- c("Absent", "Sublocal", "Local", "Supralocal")
authority_levels <-
  factor(
    authority_levels,
    levels = authority_levels,
    ordered = TRUE
  )

# join datasets
data <-
  left_join(
    political_authority,
    religious_authority,
    by = "Language"
  ) %>%
  # rename and mutate variables in dataset
  transmute(
    language = Language,
    political_authority = authority_levels[`Political Authority` + 1],
    religious_authority = authority_levels[`Religious Authority` + 1]
  ) %>%
  # get dataset as data frame
  as.data.frame()

# load phylogenetic samples from github
trees <-
  ape::read.nexus(
    file = paste0(
      "https://raw.githubusercontent.com/ScottClaessens/phyloAuthority/main/",
      "data/authority_97_20_09_21.trees"
    )
  )

# get maximum clade credibility tree
phylogeny <- phangorn::mcc(trees)

# put together list
authority <- list(
  data = data,
  phylogeny = phylogeny
  )

# save
usethis::use_data(authority, overwrite = TRUE)
