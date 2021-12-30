GetInputCatalogAsMatrix <- function(input.catalog) {
  if (mode(input.catalog) == "character") {
    input.catalog.data <- data.table::fread(input.catalog)
    # If it is an ICAMS catalog on disk it will have rownames,
    # which we have to discard.
    classes <- unlist(lapply(input.catalog.data, class))
    rowname.cols <- which(classes == "character") # Assume these are for rownames
    if (rowname.cols != 1:length(rowname.cols)) {
      stop("The character columns in ", input.catalog,
           " do not appear to be intended as row names")
    }
    input.catalog.data <- input.catalog.data[ , -rowname.cols]
    classes <- classes[-rowname.cols]
    if (any(classes != "integer")) {
      stop("The matrix portion of ", input.catalog,
           " must be all integers, got ",
           paste(unique(classes), collapse = ", "))
    }
    input.catalog <- as.matrix(input.catalog.data)
  }
  if (any(input.catalog < 0)) {
    stop("Elements < 0 found in input.catalog")
  }
  non.pos <- which(colSums(input.catalog) <= 0)
  if (length(non.pos) > 0) {
    warning("removing columns with sums <= 0: ", paste(non.pos, collapse = ", "))
    input.catalog <- input.catalog[ , -non.pos]
  }
  return(input.catalog)
}
