# ALL COLUMNS (name, HexA, HexN, Ac, S, formula, neutral_mass, floating_Na, floating_NH3)
quanSearch <- function(iso, ppm, db, shift) {
  
  iso$exp_mass <- iso$mono_mw - shift
  outp <- filter(db, ppm >= abs((neutral_mass - iso_mass) / iso_mass * 10^6)) |>
          mutate(ppm = abs(round(((neutral_mass - iso_mass) / iso_mass * 10^6), 2))) |>
  mutate(
    peak_no = iso$peak.No[i], 
    charge = iso$charge[i], 
    mz = iso$mz[i],
    mono_mw = iso$mono_mw[i],
    exp_mass = exp_mass,
    abundance = iso$abundance[i],
    scan_range = iso$scan_range[i],
    scan_count = iso$scan_count[i],
    time = iso$time[i],
    NH3 = shifts$NH3[x],
    Na = shifts$Na[x],
    Mn = shifts$Mn[x],
    FA = shifts$FA[x]
  )
  
  return(outp)
}