# Extract coverage function
extract_coverage <- function(filepath) {
  folder_name <- basename(dirname(filepath))
  content <- tryCatch({
    readLines(filepath)
  }, error = function(e) {
    warning("Could not read file: ", filepath)
    return(NULL)
  })
  
  if (is.null(content)) return(NULL)
  header_line_idx <- grep("GENOME_TERRITORY\\s+MEAN_COVERAGE", content)
  
  if (length(header_line_idx) == 0) {
    warning("Could not find metrics header in file:", filepath)
    return(NULL)
  }
  data_line_idx <- header_line_idx + 1
  
  if (data_line_idx > length(content)) {
    warning("No data line after header in file:", filepath)
    return(NULL)
  }
  header_line <- content[header_line_idx]
  data_line <- content[data_line_idx]
  headers <- unlist(strsplit(trimws(header_line), "\\s+"))
  values <- unlist(strsplit(trimws(data_line), "\\s+"))
  mean_cov_idx <- which(headers == "MEAN_COVERAGE")
  median_cov_idx <- which(headers == "MEDIAN_COVERAGE")
  mean_coverage <- NA
  median_coverage <- NA
  
  if (length(mean_cov_idx) > 0 && mean_cov_idx <= length(values)) {
    mean_coverage <- as.numeric(values[mean_cov_idx])
  }
  
  if (length(median_cov_idx) > 0 && median_cov_idx <= length(values)) {
    median_coverage <- as.numeric(values[median_cov_idx])
  }
  
  result <- data.frame(
    Folder = folder_name,
    File = basename(filepath),
    MeanCoverage = ifelse(!is.na(mean_coverage), round(mean_coverage, 2), NA),
    MedianCoverage = ifelse(!is.na(median_coverage), round(median_coverage, 2), NA),
    Tissue = ifelse(startsWith(folder_name, "CSF"), "CSF", NA),
    Batch = NA,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# Extract metrics function
extract_metrics <- function(filepath) {
  folder_name <- basename(dirname(filepath))
  
  tissue <- ifelse(startsWith(folder_name, "CSF"), "CSF", NA)
  data <- read.delim(filepath, comment.char = "#", stringsAsFactors = FALSE)
  
  pair_row <- data[data$CATEGORY == "PAIR", ]
  
  if (nrow(pair_row) > 0) {
    total_reads <- as.numeric(pair_row$TOTAL_READS)
    pf_reads_aligned <- as.numeric(pair_row$PF_READS_ALIGNED)
    usable_percentage <- (pf_reads_aligned / total_reads) * 100
    
    # high-quality aligned bases (≥ Q20) and mean read length
    pf_hq_aligned_q20_bases <- as.numeric(pair_row$PF_HQ_ALIGNED_Q20_BASES)
    mean_read_length <- as.numeric(pair_row$MEAN_READ_LENGTH)
    
    # (depth ≥10)
    reads_with_depth <- if (!is.na(mean_read_length)) {
      floor(pf_hq_aligned_q20_bases / mean_read_length)  # approx depth
    } else {
      NA
    }
    
    if (!is.na(reads_with_depth) && reads_with_depth >= 10) {
      reads_meeting_criteria <- reads_with_depth
      percentage_reads_meeting_criteria <- (reads_meeting_criteria / total_reads) * 100
    } else {
      reads_meeting_criteria <- 0
      percentage_reads_meeting_criteria <- 0
    }
    
    return(data.frame(
      Folder = folder_name,
      File = basename(filepath),
      TotalReads = total_reads,
      PFReadsAligned = pf_reads_aligned,
      UsablePercentage = round(usable_percentage, 2),
      ReadsMeetingCriteria = reads_meeting_criteria,
      PercentageReadsMeetingCriteria = round(percentage_reads_meeting_criteria, 2),
      Tissue = tissue
    ))
  } else {
    return(NULL)
  }
}

