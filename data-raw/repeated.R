library(ape)
library(tidyverse)

# load data
data <-
  read.table(
    file = paste0(
      "https://raw.githubusercontent.com/MPCMEvolution/MPCMArchive/master/",
      "online_practical_material/data_files/ch11/data_repeat.txt"
    ),
    header = TRUE
  ) %>%
  # rename variables
  transmute(
    species = species,
    x = cofactor,
    y = phen
  ) %>%
  # only include species 1-20
  filter(species %in% paste0("sp_", 1:20))

# load phylogeny
phylogeny <-
  read.nexus(
    file = paste0(
      "https://raw.githubusercontent.com/MPCMEvolution/MPCMArchive/master/",
      "online_practical_material/data_files/ch11/phylo.nex"
    )
  )

# prune phylogeny to species in dataset
phylogeny <- keep.tip(phylogeny, unique(data$species))

# put together list
repeated <- list(
  data = data,
  phylogeny = phylogeny
)

# save
usethis::use_data(repeated, overwrite = TRUE)
