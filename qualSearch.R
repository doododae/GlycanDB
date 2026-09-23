qualSearch <- function(mz, charge, ppm, iso_peak, path, na, nh, mn, fa) {

  db_path = path
  
  db = read.table(file = db_path , sep = '\t', header = TRUE)
  
  outp_db <- db
  
  shifts <- getShifts(na, nh, mn, fa)
  
  if(!is.na(mz) && !is.na(charge) && !is.na(ppm)) {
    result <- data.frame(stringsAsFactors = FALSE) 
    for(x in c(1:nrow(shifts))) {
      #experimental mass formula --- (mz * charge + charge * 1.0078) - shift
      exp_mass = round((mz * charge) + (charge * 1.0078), 5) - shifts$shift[x]
      
      if(iso_peak == 'no') {
        for(i in c(1:5)) {
          peak_mass = round(exp_mass - (i-1) * 1.00335, 5)
          
          data <- filter(db, ppm >= abs((neutral_mass - peak_mass) / peak_mass * 10^6)) |>
                  mutate("Exp Mass" = peak_mass + shifts$shift[x]) |>
                  mutate("Iso Peak" = i) |>
                  mutate("Na" = shifts$Na[x]) |>
                  mutate("NH3" = shifts$NH3[x]) |>
                  mutate("Mn" = shifts$Mn[x]) |>
                  mutate("FA" = shifts$FA[x]) |>
                  mutate(ppm = abs(round(((neutral_mass - peak_mass + shifts$shift[x]) / (peak_mass + shifts$shift[x]) * 10^6), 2)))
          
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
        data <- filter(db, ppm >= abs((neutral_mass - exp_mass) / exp_mass * 10^6)) |>
                    mutate(exp_mass = exp_mass + shifts$shift[x]) |>
                    mutate("Na" = shifts$Na[x]) |>
                    mutate("NH3" = shifts$NH3[x]) |>
                    mutate("Mn" = shifts$Mn[x]) |>
                    mutate("FA" = shifts$FA[x]) |>
                    mutate(ppm = abs(round(((neutral_mass - exp_mass + shifts$shift[x]) / (exp_mass + shifts$shift[x]) * 10^6), 2)))
        if(count(data) > 0) {
          result <- bind_rows(result, data)
          out_db <- result
        }
        else {
          result <- result
        }
        outp_db <- result
      }
    }
  }

  outp_db <- reactable(outp_db, 
    columns = list(
      name = colDef(minWidth = 120),
      neutral_mass = colDef(minWidth = 120),
      DP = colDef(show = FALSE),
      formula = colDef(show = FALSE)
    ),
    defaultColDef = colDef(
      show = TRUE, 
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