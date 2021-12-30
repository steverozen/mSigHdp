GetPossibleICAMSCatalog <- function(input.catalog) {
  if (mode(input.catalog) == "character") {
      input.catalog <- ICAMS::ReadCatalog(input.catalog, stop.on.error = FALSE)
      if (is.na(input.catalog[1,1])) {
        input.catalog <- data.table::fread(input.catalog)
      }
  }
  return(input.catalog)
}
