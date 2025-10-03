#' Path to packaged extdata
#' @param ... Path components within inst/extdata
#' @return A file path to an extdata resource
#' @export
data_file <- function(...) {
  pkg <- utils::packageName()
  p <- system.file("extdata", ..., package = pkg)
  if (nzchar(p)) return(p)
  # fallback when running from source
  file.path(here::here("inst", "extdata", ...))
}
