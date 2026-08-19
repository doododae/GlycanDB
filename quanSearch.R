# ALL COLUMNS (name, HexA, HexN, Ac, S, formula, neutral_mass, floating_Na, floating_NH3)
quanSearch <- function(iso_mass, ppm, db) {
  outp <- db
  result <- data.frame(stringsAsFactors = FALSE) 
  
  calc_mass = (10^6 * iso_mass) / (ppm + 10^6)
  n_calc_mass = (10^6 * iso_mass) / (-ppm + 10^6)
  
  data <- filter(db, abs(neutral_mass - iso_mass) < 0.5 & between(neutral_mass, calc_mass, n_calc_mass)) |>
          mutate(ppm = abs(round(((iso_mass - neutral_mass) / neutral_mass * 10^6), 2)))
  
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