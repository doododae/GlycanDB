qualSearch <- function(mz, charge, ppm, iso_peak, db) {
  outp_db <- db
  
  if(!is.na(mz) && !is.na(charge) && !is.na(ppm)) {
    #measured mass formula
    mea_mass = round((mz * charge) + (charge * 1.0078), 5)
    #calc_mass for lower bound
    calc_mass = (10^6 * mea_mass) / (ppm + 10^6)
    #calc_mass for upper bound
    n_calc_mass = (10^6 * mea_mass) / (-ppm + 10^6)
    
    if(iso_peak == 'no') {
      result <- data.frame(stringsAsFactors = FALSE) 
      for(i in c(1:5)) {
        peak_mass = round(mea_mass - (i-1) * 1.00335, 5)
        pcalc_mass = (10^6 * peak_mass) / (ppm + 10^6)
        p_ncalc_mass = (10^6 * peak_mass) / (-ppm + 10^6)
        
        data <- filter(db, between(neutral_mass, pcalc_mass, p_ncalc_mass)) |>
                mutate(experimental_mass = peak_mass) |>
                mutate(isotopic_peak = i) |>
                mutate(ppm = abs(round(((peak_mass - neutral_mass) / neutral_mass * 10^6), 2)))
        
        if(count(data) > 0) {
          result <- bind_rows(result, data)
          outp_db <- result
        }
        else {
          result <- result
        }
      }
      outp_db <- result
    }
    else {
      outp_db <- filter(db, between(neutral_mass, calc_mass, n_calc_mass)) |>
                 mutate(experimental_mass = mea_mass) |>
                 mutate(ppm = abs(round(((mea_mass - neutral_mass) / neutral_mass * 10^6), 2)))
    }
  }
  
  outp_db <- reactable(outp_db, 
    # ALL COLUMNS (name, HexA, HexN, Ac, S, formula, neutral_mass, floating_Na, floating_NH3)
    columns = list(
      name = colDef(minWidth = 120),
      neutral_mass = colDef(minWidth = 120),
      experimental_mass = colDef(minWidth = 120),
      DP = colDef(show = F),
      formula = colDef(show = F),
      floating_Na = colDef(show = F),
      floating_NH3 = colDef(show = F)
    ),
    defaultColDef = colDef(
      show = T, 
      minWidth = 60
    ), 
    details = colDef(
      name = "More",
      details = JS("function(rowInfo) {
        return `Details for row: ${rowInfo.index}` +
          `<pre>${JSON.stringify(rowInfo.values, null, 2)}</pre>`
      }"),
      html = TRUE,
      width = 60
    )
  )
  
  return(outp_db)
}