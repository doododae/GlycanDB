parse_formula <- function(f) {
  tokens <- regmatches(f, gregexpr("[A-Z][a-z]?[0-9]*", f))[[1]]
  el <- sub("[0-9]+$", "", tokens)
  n  <- as.integer(sub("^[A-Za-z]+", "", tokens))
  n[is.na(n)] <- 1L                # "N" with no number means 1
  tapply(n, el, sum)               # also merges repeated elements
}

add_formulas <- function(...) {
  counts <- unlist(lapply(c(...), function(f) as.list(parse_formula(f))))
  tot <- tapply(counts, names(counts), sum)
  
  # Hill order: C, H, then the rest alphabetically (no C: all alphabetical)
  el <- names(tot)
  ord <- if ("C" %in% el) {
    c("C", if ("H" %in% el) "H", sort(setdiff(el, c("C", "H"))))
  } else sort(el)
  
  n <- as.vector(tot[ord])
  paste0(ord, ifelse(n == 1, "", n), collapse = "")
}

sub_formulas <- function(a, b) {
  pa <- parse_formula(a); pb <- parse_formula(b)
  el <- union(names(pa), names(pb))
  res <- setNames(rep(0L, length(el)), el)
  res[names(pa)] <- res[names(pa)] + pa
  res[names(pb)] <- res[names(pb)] - pb
  res[res != 0]
}