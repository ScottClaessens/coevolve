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
  ) |>
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
  ) |>
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
  ) |>
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
  ) |>
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
  ) |>
  # rename columns
  transmute(
    species = `Row Labels`,
    group_size = `Average of Group Size`
  )

# put datasets together
data <-
  brain_data |>
  left_join(body_data, by = "species") |>
  left_join(diet_data, by = "species") |>
  left_join(system_data, by = "species") |>
  left_join(group_size_data, by = "species") |>
  # remove one species not present in phylogeny
  filter(species != "Hoolock_hoolock")

# add clade information according to groupings here:
# https://10ktrees.nunn-lab.org/Primates/downloadTrees.php
clade <-
  tribble(
    ~species,                       ~clade,
    "Alouatta_caraya",              "Atelidae",
    "Alouatta_palliata",            "Atelidae",
    "Alouatta_pigra",               "Atelidae",
    "Alouatta_seniculus",           "Atelidae",
    "Aotus_azarai",                 "Cebidae",
    "Aotus_lemurinus",              "Cebidae",
    "Aotus_trivirgatus",            "Cebidae",
    "Arctocebus_calabarensis",      "Loridae",
    "Ateles_belzebuth",             "Atelidae",
    "Ateles_fusciceps",             "Atelidae",
    "Ateles_geoffroyi",             "Atelidae",
    "Ateles_paniscus",              "Atelidae",
    "Avahi_laniger",                "Lemuriformes",
    "Avahi_occidentalis",           "Lemuriformes",
    "Brachyteles_arachnoides",      "Atelidae",
    "Cacajao_calvus",               "Pitheciidae",
    "Callicebus_moloch",            "Pitheciidae",
    "Callicebus_torquatus",         "Pitheciidae",
    "Callimico_goeldii",            "Cebidae",
    "Callithrix_argentata",         "Cebidae",
    "Callithrix_geoffroyi",         "Cebidae",
    "Callithrix_jacchus",           "Cebidae",
    "Callithrix_penicillata",       "Cebidae",
    "Callithrix_pygmaea",           "Cebidae",
    "Cebus_albifrons",              "Cebidae",
    "Cebus_apella",                 "Cebidae",
    "Cebus_capucinus",              "Cebidae",
    "Cebus_olivaceus",              "Cebidae",
    "Cercocebus_agilis",            "Papionini",
    "Cercocebus_torquatus_atys",    "Papionini",
    "Cercocebus_galeritus",         "Papionini",
    "Cercocebus_torquatus",         "Papionini",
    "Cercopithecus_ascanius",       "Cercopithecini",
    "Cercopithecus_campbelli",      "Cercopithecini",
    "Cercopithecus_cephus",         "Cercopithecini",
    "Cercopithecus_diana",          "Cercopithecini",
    "Cercopithecus_lhoesti",        "Cercopithecini",
    "Cercopithecus_mitis",          "Cercopithecini",
    "Cercopithecus_mona",           "Cercopithecini",
    "Cercopithecus_neglectus",      "Cercopithecini",
    "Cercopithecus_nictitans",      "Cercopithecini",
    "Cercopithecus_pogonias",       "Cercopithecini",
    "Cheirogaleus_major",           "Lemuriformes",
    "Cheirogaleus_medius",          "Lemuriformes",
    "Chiropotes_satanas",           "Pitheciidae",
    "Chlorocebus_aethiops",         "Cercopithecini",
    "Chlorocebus_pygerythrus",      "Cercopithecini",
    "Colobus_angolensis",           "Colobinae",
    "Colobus_guereza",              "Colobinae",
    "Colobus_polykomos",            "Colobinae",
    "Colobus_satanas",              "Colobinae",
    "Daubentonia_madagascariensis", "Lemuriformes",
    "Erythrocebus_patas",           "Cercopithecini",
    "Eulemur_coronatus",            "Lemuriformes",
    "Eulemur_fulvus_fulvus",        "Lemuriformes",
    "Eulemur_macaco_macaco",        "Lemuriformes",
    "Eulemur_mongoz",               "Lemuriformes",
    "Eulemur_rubriventer",          "Lemuriformes",
    "Euoticus_elegantulus",         "Galagonidae",
    "Galago_alleni",                "Galagonidae",
    "Galago_moholi",                "Galagonidae",
    "Galago_senegalensis",          "Galagonidae",
    "Galagoides_demidoff",          "Galagonidae",
    "Galagoides_zanzibaricus",      "Galagonidae",
    "Gorilla_beringei",             "Hominoidea",
    "Gorilla_gorilla_gorilla",      "Hominoidea",
    "Hapalemur_griseus",            "Lemuriformes",
    "Hylobates_agilis",             "Hominoidea",
    "Hylobates_klossii",            "Hominoidea",
    "Hylobates_lar",                "Hominoidea",
    "Hylobates_moloch",             "Hominoidea",
    "Hylobates_muelleri",           "Hominoidea",
    "Hylobates_pileatus",           "Hominoidea",
    "Indri_indri",                  "Lemuriformes",
    "Lagothrix_lagotricha",         "Atelidae",
    "Lemur_catta",                  "Lemuriformes",
    "Leontopithecus_chrysomelas",   "Cebidae",
    "Leontopithecus_rosalia",       "Cebidae",
    "Lepilemur_mustelinus",         "Lemuriformes",
    "Lepilemur_ruficaudatus",       "Lemuriformes",
    "Lophocebus_albigena",          "Papionini",
    "Lophocebus_aterrimus",         "Papionini",
    "Loris_tardigradus",            "Loridae",
    "Macaca_arctoides",             "Papionini",
    "Macaca_fascicularis",          "Papionini",
    "Macaca_fuscata",               "Papionini",
    "Macaca_mulatta",               "Papionini",
    "Macaca_nemestrina",            "Papionini",
    "Macaca_nigra",                 "Papionini",
    "Macaca_radiata",               "Papionini",
    "Macaca_silenus",               "Papionini",
    "Macaca_sinica",                "Papionini",
    "Macaca_sylvanus",              "Papionini",
    "Mandrillus_leucophaeus",       "Papionini",
    "Mandrillus_sphinx",            "Papionini",
    "Microcebus_murinus",           "Lemuriformes",
    "Microcebus_rufus",             "Lemuriformes",
    "Miopithecus_talapoin",         "Cercopithecini",
    "Mirza_coquereli",              "Lemuriformes",
    "Nasalis_larvatus",             "Colobinae",
    "Nomascus_concolor",            "Hominoidea",
    "Nycticebus_coucang",           "Loridae",
    "Nycticebus_pygmaeus",          "Loridae",
    "Otolemur_crassicaudatus",      "Galagonidae",
    "Otolemur_garnettii",           "Galagonidae",
    "Pan_paniscus",                 "Hominoidea",
    "Pan_troglodytes_troglodytes",  "Hominoidea",
    "Papio_anubis",                 "Papionini",
    "Papio_cynocephalus",           "Papionini",
    "Papio_hamadryas",              "Papionini",
    "Papio_papio",                  "Papionini",
    "Papio_ursinus",                "Papionini",
    "Perodicticus_potto",           "Loridae",
    "Piliocolobus_badius",          "Colobinae",
    "Pithecia_pithecia",            "Pitheciidae",
    "Pongo_abelii",                 "Hominoidea",
    "Pongo_pygmaeus",               "Hominoidea",
    "Presbytis_comata",             "Colobinae",
    "Presbytis_melalophos",         "Colobinae",
    "Procolobus_verus",             "Colobinae",
    "Propithecus_diadema",          "Lemuriformes",
    "Propithecus_verreauxi",        "Lemuriformes",
    "Pygathrix_nemaeus",            "Colobinae",
    "Saguinus_fuscicollis",         "Cebidae",
    "Saguinus_geoffroyi",           "Cebidae",
    "Saguinus_midas",               "Cebidae",
    "Saguinus_mystax",              "Cebidae",
    "Saguinus_oedipus",             "Cebidae",
    "Saimiri_boliviensis",          "Cebidae",
    "Saimiri_oerstedii",            "Cebidae",
    "Saimiri_sciureus",             "Cebidae",
    "Semnopithecus_entellus",       "Colobinae",
    "Symphalangus_syndactylus",     "Hominoidea",
    "Tarsius_bancanus",             "Tarsiidae",
    "Tarsius_syrichta",             "Tarsiidae",
    "Theropithecus_gelada",         "Papionini",
    "Trachypithecus_cristatus",     "Colobinae",
    "Trachypithecus_johnii",        "Colobinae",
    "Trachypithecus_obscurus",      "Colobinae",
    "Trachypithecus_phayrei",       "Colobinae",
    "Trachypithecus_pileatus",      "Colobinae",
    "Trachypithecus_vetulus",       "Colobinae",
    "Varecia_variegata_variegata",  "Lemuriformes"
  )

# link clade data
data <-
  data |>
  left_join(clade, by = "species") |>
  dplyr::select(species, clade, everything()) |>
  as.data.frame()

# load consensus phylogeny of primates from 10k trees website
phylogeny <- read.nexus(file = "data-raw/primates.nex")

# prune phylogeny to species in dataset
phylogeny <- keep.tip(phylogeny, data$species)

# put together
primates <- list(
  data = data,
  phylogeny = phylogeny
)

# save
usethis::use_data(primates, overwrite = TRUE)
