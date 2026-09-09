getDetails <- function(structure, data) {
  if(hasName(data, "floating_Na")) {
    outp <- filter(data, name == structure) |>
      select(name, floating_Na, floating_NH3, charge, Exep.isotopic.mz, mono_mw, time, abundance, score)
  } else if(hasName(data, "floating_Mn")) {
    outp <- filter(data, name == structure) |>
      select(name, floating_Mn, floating_NH3, charge, Exep.isotopic.mz, mono_mw, time, abundance, score)
  } else {
    outp <- filter(data, name == structure) |>
      select(name, floating_NH3, charge, Exep.isotopic.mz, mono_mw, time, abundance, score)
  }
  return(outp)
}