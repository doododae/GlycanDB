# ALL COLUMNS (name, HexA, HexN, Ac, S, formula, neutral_mass, floating_Na, floating_NH3)
quanSearch <- function(iso_mass, ppm, db, shift) {
  
  outp <- filter(db, ppm >= abs((neutral_mass - iso_mass) / iso_mass * 10^6)) |>
          mutate(ppm = abs(round(((neutral_mass - iso_mass) / iso_mass * 10^6), 2)))
  
  return(outp)
}