getDetails <- function(structure, data) {
  outp <- filter(data, name == structure) |>
    select(name, NH3, Na, Mn, FA, neutral_mass, exp_mass, charge, Exep.isotopic.mz, mono_mw, ppm, time, abundance, score)
  return(outp)
}