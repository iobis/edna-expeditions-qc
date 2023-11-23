#' @title eDNA species list quality control
#'
#' @description eDNA species list quality control
#'
#' @docType package
#' @name ednaqc
#' @import dplyr glue worrms speedy stringr purrr memoise
#' @author Pieter Provoost, \email{p.provoost@unesco.org}
NULL

find_plausibility <- function(aphiaid, decimalLongitude, decimalLatitude) {
  message(glue("Calculating plausibility for {aphiaid} at POINT({decimalLongitude} {decimalLatitude})"))
  pl <- get_dist(aphiaid = aphiaid) %>%
    calculate_plausibility()
  pl_score <- pl %>%
    st_intersection(st_sfc(st_point(c(decimalLongitude, decimalLatitude)), crs = 4326)) %>%
    pull(plausibility)
  if (length(pl_score) == 0) {
    pl_score <- NA
  }
  pl_score
}

find_plausibility_cached <- memoise(find_plausibility)

analyze <- function(occurrence) {

  occurrence <- occurrence %>%
    mutate(aphiaid = as.numeric(str_extract(scientificNameID, "[0-9]+")))

  # aggregate

  occurrence <- occurrence %>%
    group_by(aphiaid, decimalLongitude, decimalLatitude, phylum, class, order, family, genus, scientificName) %>%
    summarize()

  # add WoRMS flags

  aphiaids <- unique(occurrence$aphiaid)
  aphiaid_batches <- split(aphiaids, as.integer((seq_along(aphiaids) - 1) / 50))
  aphia_records <- map(aphiaid_batches, wm_record) %>%
    bind_rows() %>%
    select(aphiaid = AphiaID, isMarine, isBrackish, isTerrestrial, isFreshwater, isExtinct)

  occurrence <- occurrence %>%
    left_join(aphia_records, by = "aphiaid")

  # check plausibility

  plausibilities <- map(1:nrow(occurrence), function(i) {
    suppressMessages(find_plausibility_cached(occurrence$aphiaid[i], occurrence$decimalLongitude[i], occurrence$decimalLatitude[i]))
  }, .progress = TRUE)




}
