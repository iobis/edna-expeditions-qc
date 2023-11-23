library(dplyr)
library(purrr)
library(speedy)
library(stringr)
library(ggplot2)
library(sf)
library(glue)

st <- storr::storr_rds("storr")
st_qc <- storr::storr_rds("storr_qc")

# read eDNA results

dna_files <- list.files("edna-results/data", "*DNADerivedData*", full.names = TRUE)
occurrence_files <- list.files("edna-results/data", "*Occurrence*", full.names = TRUE)

dna <- map(dna_files, read.table, sep = "\t", quote = "", header = TRUE) %>%
  bind_rows() %>%
  mutate_if(is.character, na_if, "")

occurrence <- map(occurrence_files, read.table, sep = "\t", quote = "", header = TRUE) %>%
  bind_rows() %>%
  mutate_if(is.character, na_if, "") %>%
  mutate(
    aphiaid = str_replace(scientificNameID, "urn:lsid:marinespecies.org:taxname:", ""),
    species = ifelse(taxonRank == "species", scientificName, NA)
  ) %>%
  left_join(dna, by = "occurrenceID")

# analysis

detections <- occurrence %>%
  filter(!is.na(species) & !is.na(decimalLongitude)) %>%
  group_by(higherGeography, decimalLongitude, decimalLatitude, species, aphiaid) %>%
  summarize() %>%
  ungroup()

check_occurrence <- function(aphiaid, decimalLongitude, decimalLatitude) {
  # key <- glue("{aphiaid}_{decimalLongitude}_{decimalLatitude}")
  key <- paste0(aphiaid, "_", decimalLongitude, "_", decimalLatitude)
  message(key)
  if (!st_qc$exists(key)) {
    gc()
    gc()
    if (!st$exists(aphiaid)) {
      dist <- get_dist(aphiaid = as.numeric(aphiaid))
      st$set(aphiaid, dist)
    } else {
      dist <- st$get(aphiaid)
    }
    # TODO: cleanup below
    if (nrow(dist$gbif) == 0) {
      dist$gbif <- NULL
    } else {
      dist$gbif <- st_set_crs(dist$gbif, 4326)
    }
    if (nrow(dist$obis) == 0) {
      dist$obis <- NULL
    } else {
      dist$obis <- st_set_crs(dist$obis, 4326)
    }
    dist$worms <- st_sf(st_sfc(), crs = 4326)
    if (is.null(dist$envelope$envelope)) {
      dist$envelope$envelope <- st_sf(st_sfc(), crs = 4326)
    }
    point <- st_sfc(st_point(c(decimalLongitude, decimalLatitude)), crs = 4326)
    db <- bind_rows(dist$obis, dist$gbif)
    if (nrow(db) > 0) {
      d <- as.numeric(min(st_distance(point, db)))
    } else {
      d <- Inf
    }
    worms_intersect <- st_intersection(dist$worms, point)
    envelope_intersect <- st_intersection(dist$envelope$envelope, point)
    result <- list(distance = d, worms = as.logical(nrow(worms_intersect)), thermal = as.logical(nrow(envelope_intersect)))
    st_qc$set(key, result)
  }
}


for (i in 1:nrow(detections)) {
  check_occurrence(detections$aphiaid[i], detections$decimalLongitude[i], detections$decimalLatitude[i])
}







# i <- 578
# aphiaid <- detections$aphiaid[i]
# decimalLongitude <- detections$decimalLongitude[i]
# decimalLatitude <- detections$decimalLatitude[i]
#
# plot_dist(dist) +
#   geom_sf(data = point, shape = 21, size = 2, color = "red", stroke = 2)


#######!!!!!!!!!!!!!!!!!
  # fix unaccepted https://www.marinespecies.org/aphia.php?p=taxdetails&id=400783&from=rss
# in resolve_taxonomy
# add in data files / species lists as well?
# also: list vsearch consensus


