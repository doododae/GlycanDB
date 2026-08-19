quanSummary <- function(data) {
  if(nrow(data) > 1) {
    outp <- data %>%
      group_by(name) %>%
      summarise(
        abundance = round(sum(abundance), 0),
        score = round(max(score), 2)
      )  
  }
  else {
    outp <- data
  }
  return(outp)
}