IsICAMSCatalog <- function(matrix) {
  # Fragile, to be replaced by an ICAMS function
  # add in more conditions here because it breaks when input.catalog
  # is a subset of an ICAMS catalog
  # e.g. input.catalog <- ICAMS.catalog[1:10,]
  if (!is.matrix(matrix)) return(FALSE)
  return(any(grepl("Catalog", class(matrix)))&&
           nrow(matrix) %in% c(96,192,1536,78,144,136,83,166,1697))

}
