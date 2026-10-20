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
  
  format_formula(tot)
}

format_formula <- function(counts) {
  counts <- counts[counts != 0]
  el <- names(counts)
  ord <- if ("C" %in% el) {
    c("C", if ("H" %in% el) "H", sort(setdiff(el, c("C", "H"))))
  } else sort(el)
  n <- as.vector(counts[ord])
  paste0(ord, ifelse(n == 1, "", n), collapse = "")
}

sub_formulas <- function(a, b) {
  pa <- parse_formula(a); pb <- parse_formula(b)
  el <- union(names(pa), names(pb))
  res <- setNames(rep(0L, length(el)), el)
  res[names(pa)] <- res[names(pa)] + as.vector(pa)
  res[names(pb)] <- res[names(pb)] - as.vector(pb)
  
  if (any(res < 0)) {
    warning("'", b, "' is not a subset of '", a, "'")
    return(NA_character_)
  }
  format_formula(res)
}

format_db <- function(db, reducing_end, non_end) {
  og_data <- db
  
  #add any reducing/non-reducing ends to the backbone database
  chem_mod <- data.frame(
    name = c("pNP", "DUA", "1,6-Anhydro"),
    mass = c(121.0164, -18.0106, -18.0106),
    formula = c("C6H3NO2", "H2O", "H20")
  )
  
  if('pnp' %in% reducing_end) {
    mod_formula <- chem_mod %>% filter(name == "pNP") %>% pull(formula)
    tag_mass <- chem_mod %>% filter(name == "pNP") %>% pull(mass)
    db$name <- paste0(db$name, "_pNP")
    db$formula <- mapply(add_formulas, db$formula, mod_formula)
    db$neutral_mass <- db$neutral_mass + tag_mass
  } else if('1-6-anhydro' %in% reducing_end) {
    mod_formula <- chem_mod %>% filter(name == "1-6-anhydro") %>% pull(formula)
    tag_mass <- chem_mod %>% filter(name == "1-6-anhydro") %>% pull(mass)
    db$name <- paste0(db$name, "_1,6-Anhydro")
    db$formula <- mapply(sub_formulas, db$formula, mod_formula)
    db$neutral_mass <- db$neutral_mass + tag_mass
  }
  
  if('normal' %in% reducing_end) {
    temp <- db
    db <- union(db, og_data)
    #save reducing end only data in case needed
    og_data <- temp
  }
  
  if('dua' %in% non_end) {
    mod_formula <- chem_mod %>% filter(name == "DUA") %>% pull(formula)
    tag_mass <- chem_mod %>% filter(name == "DUA") %>% pull(mass)
    db$name <- sub("HS_", "ΔUA_", db$name, fixed = TRUE)
    db$formula <- mapply(sub_formulas, db$formula, mod_formula)
    db$neutral_mass <- db$neutral_mass + tag_mass
  }
  
  if('normal' %in% non_end) {
    db <- union(db, og_data)
  }
  
  return(db)
}