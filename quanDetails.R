getDetails <- function(structure, data) {
  outp <- filter(data, name == structure) |>
          select(name, floating_Na, floating_NH3, charge, Exep.isotopic.mz, mono_mw, time, abundance, score)
  return(outp)
}