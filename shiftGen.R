getShifts <- function(na, nh, mn, fa) {
  #create an adduct shifts combination matrix with user input
  
  adducts <- data.frame(
    name = c("NH3", "Na", "Mn", "FA"),
    mass = c(17.0265, 22.9898, 54.9380, 46.0055),
    h_replace = c(0, 1, 2, 0)
  )
  
  max_depth <- data.frame(
    NH3 = nh,
    Na = na,
    Mn = mn,
    FA = fa
  )
  
  counts <- expand.grid(
    setNames(lapply(adducts$name, function(n) 0:max_depth[[n]]), adducts$name)
  )
  
  per_unit <- adducts$mass - 1.0078 * adducts$h_replace
  
  counts$shift <- as.vector(as.matrix(counts) %*% per_unit)
  
  counts$combo <- apply(counts[adducts$name], 1, function(r) {
    keep <- r >= 0
    if (!any(keep)) "none" else paste0(adducts$name[keep], "x", r[keep], collapse = " ")
  })
  
  counts <- counts[counts$combo != "none", ]
  
  return(counts)
}