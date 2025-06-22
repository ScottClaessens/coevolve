library(ape)
library(geosphere)
library(phangorn)
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
  read.nexus(
    file = paste0(
      "https://raw.githubusercontent.com/ScottClaessens/phyloAuthority/main/",
      "data/authority_97_20_09_21.trees"
    )
  )

# get maximum clade credibility tree
phylogeny <- mcc(trees)

# load longitude and latitude values
lon_lat <-
  read_tsv(
    file = paste0(
      "https://raw.githubusercontent.com/ScottClaessens/phyloAuthority/main/",
      "data/authority_coordinates.txt"
    ),
    col_types = "cd"
  ) %>%
  dplyr::select(Language, Longitude, Latitude) %>%
  column_to_rownames(var = "Language") %>%
  as.matrix()

# convert to geographic distance matrix
dist_mat <- distm(lon_lat)
rownames(dist_mat) <- rownames(lon_lat)
colnames(dist_mat) <- rownames(lon_lat)

# put together list
authority <- list(
  data = data,
  phylogeny = phylogeny,
  distance_matrix = dist_mat
)

# save
usethis::use_data(authority, overwrite = TRUE)
