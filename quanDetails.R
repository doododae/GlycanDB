getDetails <- function(structure, data) {
  outp <- filter(data, name == structure) |>
    select(name, NH3, Na, Mn, FA, mono_mw, exp_mass, Exep.isotopic.mz, charge, ppm, time, abundance, score)
  return(outp)
}