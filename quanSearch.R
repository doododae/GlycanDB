# ALL COLUMNS (name, HexA, HexN, Ac, S, formula, neutral_mass, floating_Na, floating_NH3)
quanSearch <- function(iso_mass, ppm, db, shift) {
  outp <- db
  result <- data.frame(stringsAsFactors = FALSE) 
  
  data <- filter(db, abs(iso_mass + shift - neutral_mass) < 0.5 & ppm > abs((neutral_mass - iso_mass - shift) / (iso_mass - shift) * 10^6)) |>
          mutate(ppm = abs(round(((neutral_mass - iso_mass + shift) / (iso_mass + shift) * 10^6), 2)))
  
  if (count(data) > 0) {
    result <- bind_rows(result, data)
    outp <- result
  } 
  else {
    result <- result
  }
  
  outp <- result
  return(outp)
}