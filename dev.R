library(dplyr)
library(purrr)
library(speedy)

# read eDNA results

system("git clone -b data --depth=1 https://github.com/iobis/edna-results.git")

occurrence_files <- list.files("edna-results/data", "*Occurrence*", full.names = TRUE)

occurrence <- map(occurrence_files, read.table, sep = "\t", quote = "", header = TRUE) %>%
  bind_rows() %>%
  mutate_if(is.character, na_if, "") %>%
  filter(taxonRank == "species" & scientificName != "Homo sapiens")

# run speedy

options(timeout = 120)
st <- storr::storr_rds("storr")

aphiaids <- unique(as.numeric(stringr::str_replace(occurrence$scientificNameID, "urn:lsid:marinespecies.org:taxname:", "")))

x <- map(sample(aphiaids), function(aphiaid) {
  message(aphiaid)
  if (!st$exists(aphiaid)) {
    d <- get_dist(aphiaid = aphiaid)
    st$set(aphiaid, d)
  }
  NULL
}, .progress = TRUE)

# check

d <- st$get("216683")
plot_dist(d)
