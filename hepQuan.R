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
                    window = 100, nups = 3, ndowns = 3, threshold = 2, ticfilter = .95) {
  
  db_path = "db/test_backbone_database.tsv"
  
  data = read.table(file = db_path , sep = '\t', header = TRUE)
  
  # Filter db table based on DP range
  data <- filter(data, between(data$DP, dp_lwr, dp_upr))

  # Step 1: TIC grouping
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
  iso<- iso[which(iso$peak.No!=0),]
  
  #Step2 mass grouping
  
  #round monoMW
  iso$round_mw <- round(iso$monoisotopic_mw, 4) 
  
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
  
  result <- data.frame(stringsAsFactors = FALSE)
  raw_outp <- data.frame(stringsAsFactors = FALSE)
  
  #returns shift combo matrix
  shifts <- getShifts(na = na, nh = nh, mn = mn, fa = fa)
  
  print(paste("number of shift rows:", nrow(shifts)))
  print(paste(shifts$combo, shifts$shift))
  
  for(i in c(1:nrow(iso))) {
    for(x in c(1:nrow(shifts))) {
      raw <- data.frame(stringsAsFactors = FALSE)
      print(paste(i, x))
      # + iso$charge[i] * 1.0078
      exp_mass = iso$mono_mw[i] - shifts$shift[x]
      raw <- raw %>%
        mutate(
          NH3 = shifts$NH3[x],
          Na = shifts$Na[x],
          Mn = shifts$Mn[x],
          FA = shifts$FA[x],
          shift = shifts$shift[x],
          mono_mw = iso$mono_mw[i],
          exp_mass = exp_mass
        ) %>%
        add_row(
          NH3 = shifts$NH3[x],
          Na = shifts$Na[x],
          Mn = shifts$Mn[x],
          FA = shifts$FA[x],
          shift = round(shifts$shift[x], 5),
          mono_mw = round(iso$mono_mw[i], 5),
          exp_mass = round(exp_mass, 5)
        )
      res_temp <- quanSearch(exp_mass, ppm, data, shifts$shift[x]) %>%
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
      # delete high adductive
      # Adduct (NH3+Na+Mn+FA) < S + HexA + charge – 2
      adducts = res_temp$Mn + res_temp$NH3 + res_temp$Na + res_temp$FA
      res_temp <- filter(res_temp, adducts < res_temp$S + res_temp$HexA + res_temp$charge - 2)
      result <- bind_rows(result, res_temp)
      raw_outp <- bind_rows(raw_outp, raw)
    }
  }
  
  rm(res_temp)
  
  write_csv(raw_outp, file = "shift_output.csv")

  #delete not matched peaks
  if(nrow(result) > 0) {
    result <- filter(result, result$neutral_mass != 0) 
    
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
          neutral_mass = mean(neutral_mass),
          mono_mw = round(mean(mono_mw), 4),
          exp_mass = round(mean(exp_mass), 4),
          #Adductive = mean(Adductive),
          DP = mean(DP),
          mz = round(mean(mz), 4),
          abundance = sum(abundance),
          time = median(time),
          ppm = round(mean(ppm), 2),
          gaussian = max(gaussian),
          scan_count=sum(scan_count),
          scan_range=sum(scan_range)
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
    print("no results")
    outp_db <- result
  }
  
  return(outp_db)
}
