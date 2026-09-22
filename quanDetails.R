getDetails <- function(structure, data) {
  outp <- filter(data, name == structure) |>
    select(name, NH3, Na, Mn, FA, charge, Exep.isotopic.mz, mono_mw, ppm, time, abundance, score)
  return(outp)
}