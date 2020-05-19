
eload <- function(file) {
  env <- new.env()
  load(file, envir = env)
  return(env)
}
