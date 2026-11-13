# Capitalize island labels; this helper was absent from the downloaded repository.
firstup <- function(x) {
  x <- as.character(x)
  substr(x, 1L, 1L) <- toupper(substr(x, 1L, 1L))
  x
}
