source('quanSearch.R')
source('gaussianSmooth.R')
library(data.table)

# scans.csv, isos.csv, ppm filter, database, 
# scan interval, scan start time, scan end time 
# dp lower bound, dp upper bound 
# Adducts - NA, NH3, Mn, FA
# smoothing window, findpeaks nups, findpeaks ndowns, findpeaks threshold, TIC grouping filter
hepQuan <- function(scan, iso, ppm, db, 
                    minscan, start, end, 
                    dp_lwr, dp_upr, 
                    na, nh, mn, fa, 
                    window = 100, nups = 3, ndowns = 3, threshold = 2, ticfilter = .95,
                    reducing_end, non_end) {
  
  ### LOAD & PREPARE DB ###
  incProgress(0.10, detail = "Prepping database...")
  db_path = db
  
  data = read.table(file = db_path , sep = '\t', header = TRUE)
  
  # Filter db table based on DP range
  data <- filter(data, between(data$DP, dp_lwr, dp_upr))
  
  reducing_end <- if (is.null(reducing_end)) 0 else reducing_end
  
  non_end <- if (is.null(non_end)) 0 else non_end
  
  data <- format_db(data, reducing_end, non_end)
  ### END OF DB PREP ###
  
  
  # Step 1: TIC grouping
  incProgress(0.05, detail = "TIC Grouping & Peak Smoothing...")
  scan_starting = scan$scan_num[1] - 1
  scan$scan_num <- scan$scan_num - scan_starting
  iso$scan_num <- iso$scan_num - scan_starting
  
  
  #peak smoothing, if scan num <100, skip TIC grouping step
  if(nrow(scan) > 100) {
    peaks <- gaussianSmooth(scan$tic, window)
    peaks[which(is.na(peaks))] <- 0
    peaks <- as.numeric(peaks)
    #TIC grouping
    tic_group <- findpeaks(peaks, nups = nups, ndowns = ndowns, threshold = threshold)
    tic_group <- as.data.frame(tic_group)
    
    #TIC grouping filter
    tic_group <- tic_group[which(scan$tic[tic_group$V2] * ticfilter > scan$tic[tic_group$V3]),]
  }
  else {
    tic_group <-as.data.frame(matrix(0, 1, 1))
    tic$V2 <- min(scan$scan_num)
    tic_group$V3 <- max(scan$scan_num)
  }
  
  #add peak group
  tic_group$group <- c(1:nrow(tic_group))
  
  #add peak group information of each scan numbers
  iso$peak.No <- 0
  for (i in c(1:nrow(tic_group))) {
    iso$peak.No[which(iso$scan_num>=tic_group$V3[i]&iso$scan_num<=tic_group$V4[i])] <- tic_group$group[i] 
  }
  
  #delete scan numbers which are not in peak group
  #iso<- iso[which(iso$peak.No!=0),]
  
  #Step2 mass grouping
  incProgress(0.05, detail = "Mass Grouping...")
  #round monoMW
  iso$round_mw <- round(iso$monoisotopic_mw, 2) 
  
  #add each MW counts
  iso$mw_count <- 1
  
  #add bpi and time of each scan number
  index <- iso$scan_num
  iso$time <- scan$scan_time[index]
  
  #peak from start elution time to end elution time for analysis
  iso <- iso[which(iso$time>=start),]
  iso <- iso[which(iso$time<=end),]
  
  #delete unnecessary information
  iso <- iso[,c("scan_num","charge","mz","round_mw","mw_count","abundance","peak.No","monoisotopic_mw","time")]
  
  iso_raw <- iso
  
  #combine data based mw, peak group, charge
  iso <- iso %>%
    group_by(round_mw,peak.No,charge) %>%
    summarise(
      mz = mean(mz),
      mono_mw = mean(monoisotopic_mw),
      scan_range = max(scan_num)-min(scan_num) + 1,
      abundance = sum(abundance),
      scan_count = sum(mw_count),
      time = mean(time)
    )
  
  #delete peaks with scan number
  iso <- filter(iso, scan_range >= minscan)
  
  print(paste("number of iso rows:", nrow(iso)))
  
  #returns shift combo matrix
  incProgress(0.10, detail = "Applying Adduct Shifts...")
  shifts <- getShifts(na = na, nh = nh, mn = mn, fa = fa)
  
  print(paste("number of shift rows:", nrow(shifts)))

  raw_outp <- crossing(iso, shifts) %>%
    mutate(
      NH3 = NH3,
      Na = Na,
      Mn = Mn,
      FA = FA,
      shift = shift,
      mono_mw = mono_mw,
      exp_mass = mono_mw - shift
    )
  write_csv(raw_outp, file = "shift_output.csv")
  
  # combo of all iso mono_mw and shifts
  combos <- crossing(iso, shifts) %>%
    mutate(
      adducts = Mn + NH3 + Na + FA,
      exp_mass = mono_mw - shift
    )
  
  data_idx <- data %>%
    mutate(match_row = row_number())
  
  incProgress(0.3, detail = "Matching Masses...")
  matches <- combos %>%
    rowwise() %>%
    mutate(
      match_row = list(
        which(
          ppm >= abs(data$neutral_mass - exp_mass) / exp_mass * 1e6
        )
      )
    ) %>%
    ungroup() %>%                          # ungroup here, once
    filter(lengths(match_row) > 0) %>%     # now counts matches per row correctly
    unnest(match_row) %>%                  # one row per match
    inner_join(data_idx, by = "match_row") %>%
    mutate(ppm = abs(neutral_mass - exp_mass) / exp_mass * 1e6)
  
  adduct_cols <- c("NH3", "Na", "Mn", "FA")
  
  # TRUE if the distinct values are exactly 0, 1, 2, ... max (with max > 0)
  is_seq <- function(x) {
    u <- sort(unique(x))
    max(u) > 0 && u[1] == 0 && all(diff(u) == 1)
  }
  
  matches <- matches %>%
    group_by(name) %>%
    filter(any(sapply(across(all_of(adduct_cols)), is_seq))) %>%
    ungroup()
  
  # delete high adductives
  # Adducts (NH3+Na+Mn+FA) < S + HexA + charge – 2
  result <- matches %>%
    filter(adducts < S + HexA + charge - 2) %>%
    filter(neutral_mass != 0) 
  
  print(result)
  
  #delete not matched peaks
  if(nrow(result) > 0) {

    incProgress(0.45, detail = "Finalizing...")
    
    result$gaussian <- 0.5
    
    #calculate Gaussian similarity
    length <- length(result$neutral_mass)
    
    for (i in c(1:length)) {
      data <- filter(iso_raw, round_mw == round(result$mono_mw[i], 2))
      
      if(nrow(data) > 2) {
        td <- data$time
        d <- data$abundance
        mu <- data$time[data$abundance == max(data$abundance)]
        num_peak_pts <- length(data$abundance)
        
        sigma <- max(data$time) - min(data$time)
        h <- max(data$abundance)
        
        fit <- try(nls(d ~ SSgauss(td, mu, sigma, h)), silent = TRUE)
        
        if(class(fit) != "try-error") {
          gaussPts <- as.matrix(fitted(fit))
          gaussPts_std <- (gaussPts-mean(gaussPts)) / sd(gaussPts)
          gaussPts_scale <- gaussPts_std / norm(gaussPts_std, type="F")
          
          d <- as.matrix(d)
          peak_intensity_std <- (d-mean(d)) / sd(d)
          peak_intensity_scale <- peak_intensity_std / norm(peak_intensity_std, type="F")
          
          gauss_similarity <- sum(gaussPts_scale * peak_intensity_scale)
          #result$gaussian[i] <- gauss_similarity
          
        }
        else {
          gauss_similarity <- 0.5
        }
      }
      else {
        gauss_similarity <- 0.5
      }
      
      result$gaussian[i] <- gauss_similarity
    }
    
    if(nrow(result) != 0) {
      result <- result %>%
        group_by(name, charge, NH3, Na, Mn, FA) %>%
        summarise(
          neutral_mass = round(mean(neutral_mass), 4),
          mono_mw = round(mean(mono_mw), 4),
          exp_mass = round(mean(exp_mass), 4),
          DP = mean(DP),
          mz = round(mean(mz), 4),
          abundance = sum(abundance),
          time = median(time),
          ppm = round(mean(ppm), 2),
          gaussian = max(gaussian),
          scan_count= sum(scan_count),
          scan_range= sum(scan_range)
        )
      
      if(nrow(result) > 30) {
        result <- filter(result, scan_range >= minscan)
        result$log_ms <- log2(result$neutral_mass)
        model <-  lm(log_ms~time, data = result)
        pred.int <- predict(model, interval = "prediction")
        result <- cbind(result, pred.int)
        
        #result <- filter(result, log_ms > lwr & log_ms < upr | DP < 4)
      }
      result$Exep.isotopic.mz <- round((result$mono_mw - 1.0078 * result$charge) / result$charge, 4)
    }
    else {
      result <- filter(result, scan_range >= minscan & scan_count >= 1)
    }
    
    #calculate score
    result$score <- round(result$gaussian * log10(result$abundance), 2)
    result$time <- round(result$time, 2)
    
    outp_db <- result
  }
  else {
    print("No Results")
    outp_db <- result
  }
  
  incProgress(0.5, detail = "Complete!")
  return(outp_db)
}
