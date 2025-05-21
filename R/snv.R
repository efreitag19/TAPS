# Convert annotated snpeff outputs to csv
convert_all_granges_to_csv <- function(base_dir, dir_pattern = "^CSF-", file_pattern = "annotated_snpeff.rds$") {
  # Define our core conversion function
  convert_single_file <- function(rds_file, output_csv_file) {
    tryCatch({
      # Load the RDS file
      gr_obj <- readRDS(rds_file)
      
      # First extract the basic GRanges information
      basic_info <- data.frame(
        chromosome = as.character(seqnames(gr_obj)),
        start = start(gr_obj),
        end = end(gr_obj),
        strand = as.character(strand(gr_obj))
      )
      
      # Get metadata
      meta <- mcols(gr_obj)
      meta_df <- data.frame(row.names=1:length(gr_obj))
      
      # Process each column individually with safer conversion
      for (col_name in names(meta)) {
        col_data <- meta[[col_name]]
        
        # Handle different types of data
        if (is(col_data, "List") || is.list(col_data)) {
          # Process lists element by element to avoid coercion issues
          safe_values <- character(length(col_data))
          
          for (i in seq_along(col_data)) {
            x <- col_data[[i]]
            if (length(x) == 0 || all(is.na(x))) {
              safe_values[i] <- "NA"
            } else {
              # Try to convert each element to string
              tryCatch({
                if (length(x) == 1) {
                  safe_values[i] <- as.character(x)
                } else {
                  safe_values[i] <- paste(as.character(x), collapse=",")
                }
              }, error = function(e) {
                safe_values[i] <- "ERROR_CONVERTING"
              })
            }
          }
          meta_df[[col_name]] <- safe_values
        } else {
          # For regular columns, add directly
          meta_df[[col_name]] <- col_data
        }
      }
      
      # Combine data frames
      result_df <- cbind(basic_info, meta_df)
      
      # Add the sample name (derived from the directory name)
      sample_name <- basename(dirname(rds_file))
      result_df$sample <- sample_name
      
      # Write to CSV
      write.csv(result_df, file=output_csv_file, row.names=FALSE, quote=TRUE)
      
      cat("Converted:", rds_file, "→", output_csv_file, "\n")
      return(TRUE)
    }, error = function(e) {
      cat("ERROR processing", rds_file, ":", conditionMessage(e), "\n")
      return(FALSE)
    })
  }
  
  # List all subdirectories matching the pattern
  ngs_dirs <- list.dirs(base_dir, recursive=TRUE)
  ngs_dirs <- grep(dir_pattern, basename(ngs_dirs), value=TRUE)
  ngs_dirs <- file.path(base_dir, ngs_dirs)
  
  # Find all matching RDS files
  rds_files <- c()
  for (dir in ngs_dirs) {
    matching_files <- list.files(
      dir, 
      pattern=file_pattern,
      full.names=TRUE,
      recursive=FALSE
    )
    rds_files <- c(rds_files, matching_files)
  }
  
  # Create output directory if it doesn't exist
  output_dir <- file.path(base_dir, "csv_exports")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive=TRUE)
  }
  
  # Process each file
  results <- list()
  cat("Found", length(rds_files), "files to process\n")
  
  for (rds_file in rds_files) {
    sample_name <- basename(dirname(rds_file))
    output_file <- file.path(output_dir, paste0(sample_name, "_snpeff.csv"))
    
    success <- convert_single_file(rds_file, output_file)
    results[[rds_file]] <- success
  }
  
  # If successful conversions exist, combine them into one file
  success_files <- names(results)[unlist(results)]
  if (length(success_files) > 0) {
    # Extract the prefix from dir_pattern for naming the combined file
    pattern_prefix <- gsub("\\^", "", dir_pattern)
    pattern_prefix <- gsub("-.*$", "", pattern_prefix)
    if (pattern_prefix == "") pattern_prefix <- "combined"
    
    # Create combined file path
    combined_file <- file.path(output_dir, paste0(pattern_prefix, "_combined_snpeff.csv"))
    
    # Get all successfully created CSV files
    csv_files <- list()
    for (rds_file in success_files) {
      sample_name <- basename(dirname(rds_file))
      csv_file <- file.path(output_dir, paste0(sample_name, "_snpeff.csv"))
      if (file.exists(csv_file)) {
        csv_files[[length(csv_files) + 1]] <- csv_file
      }
    }
    
    # Read and combine all CSVs if any exist
    if (length(csv_files) > 0) {
      combined_data <- NULL
      
      for (csv_file in csv_files) {
        tryCatch({
          data <- read.csv(csv_file, stringsAsFactors=FALSE)
          
          # If this is the first CSV, set combined_data to this data
          if (is.null(combined_data)) {
            combined_data <- data
          } else {
            # Check if column names match
            missing_cols <- setdiff(names(combined_data), names(data))
            if (length(missing_cols) > 0) {
              # Add missing columns with NA values
              for (col in missing_cols) {
                data[[col]] <- NA
              }
            }
            
            # Check for columns in data that aren't in combined_data
            extra_cols <- setdiff(names(data), names(combined_data))
            if (length(extra_cols) > 0) {
              # Add these columns to combined_data with NA values
              for (col in extra_cols) {
                combined_data[[col]] <- NA
              }
            }
            
            # Now columns should match - combine using rbind
            combined_data <- rbind(combined_data, data[names(combined_data)])
          }
        }, error = function(e) {
          cat("Error reading", csv_file, ":", conditionMessage(e), "\n")
        })
      }
      
      # Write combined file
      if (!is.null(combined_data) && nrow(combined_data) > 0) {
        write.csv(combined_data, file=combined_file, row.names=FALSE)
        cat("Created combined file:", combined_file, "with", nrow(combined_data), "rows\n")
      }
    }
  }
  
  # Summary
  success_count <- sum(unlist(results))
  cat("\nSummary: Successfully converted", success_count, "of", length(rds_files), "files\n")
  cat("CSV files saved to:", output_dir, "\n")
  
  return(invisible(results))
}

# Convert snv vcf outputs to csv
convert_all_vcfs_to_csv <- function(base_dir, dir_pattern = "^CSF-") {
  # Load required libraries
  library(vcfR)
  library(stringr)
  library(dplyr)
  
  # Define our core conversion function
  convert_single_file <- function(vcf_file, output_csv_file) {
    tryCatch({
      # Load the VCF file
      vcf_obj <- read.vcfR(vcf_file)
      
      # Extract sample name from the directory name
      sample_name <- basename(dirname(vcf_file))
      
      # Extract basic variant information
      variant_data <- data.frame(
        SAMPLE = sample_name,
        CHROM = getCHROM(vcf_obj),
        POS = getPOS(vcf_obj),
        REF = getREF(vcf_obj),
        ALT = getALT(vcf_obj),
        QUAL = getQUAL(vcf_obj),
        stringsAsFactors = FALSE
      )
      
      # Process INFO field to extract gene and variant type
      info_fields <- vcf_obj@fix[, "INFO"]
      
      # Initialize columns for gene information
      variant_data$GENE <- NA
      variant_data$AberrationType <- NA
      variant_data$HGVS_c <- NA
      variant_data$HGVS_p <- NA
      variant_data$GT <- NA
      variant_data$AD <- NA
      variant_data$DP <- NA
      
      # Extract information from INFO field
      for (i in 1:nrow(variant_data)) {
        info <- info_fields[i]
        
        # Extract Gene name
        gene_match <- str_match(info, "\\|([^\\|]+)\\|transcript\\|")
        if (!is.na(gene_match[1,2])) {
          variant_data$GENE[i] <- gene_match[1,2]
        }
        
        # Extract Annotation type (variant type)
        annot_match <- str_match(info, "ANN=[^\\|]+\\|([^\\|]+)\\|")
        if (!is.na(annot_match[1,2])) {
          variant_data$AberrationType[i] <- annot_match[1,2]
        }
        
        # Extract HGVS.c notation
        hgvs_c_match <- str_match(info, "c\\.([^\\|]+)\\|")
        if (!is.na(hgvs_c_match[1,1])) {
          variant_data$HGVS_c[i] <- paste0("c.", hgvs_c_match[1,1])
        }
        
        # Extract HGVS.p notation
        hgvs_p_match <- str_match(info, "p\\.([^\\|]+)\\|")
        if (!is.na(hgvs_p_match[1,1])) {
          variant_data$HGVS_p[i] <- paste0("p.", hgvs_p_match[1,1])
        }
      }
      
      # Extract FORMAT fields (GT, AD, DP)
      # First, get the FORMAT string for each variant
      format_strings <- vcf_obj@gt[, 1]
      
      # Process each variant
      for (i in 1:nrow(variant_data)) {
        format_fields <- unlist(strsplit(format_strings[i], ":"))
        sample_data <- unlist(strsplit(vcf_obj@gt[i, 2], ":"))
        
        # Match fields with their values
        for (j in 1:length(format_fields)) {
          if (j <= length(sample_data)) {
            field <- format_fields[j]
            value <- sample_data[j]
            
            if (field == "GT") {
              variant_data$GT[i] <- value
            } else if (field == "AD") {
              variant_data$AD[i] <- value
            } else if (field == "DP") {
              variant_data$DP[i] <- value
            }
          }
        }
      }
      
      # Create a copy of the data frame to work with
      processed_data <- variant_data
      
      # Safely extract AD_REF and AD_ALT
      processed_data$AD_REF <- NA
      processed_data$AD_ALT <- NA
      
      # Process AD field safely
      for (i in 1:nrow(processed_data)) {
        if (!is.na(processed_data$AD[i])) {
          ad_values <- unlist(strsplit(processed_data$AD[i], ","))
          if (length(ad_values) >= 1) {
            processed_data$AD_REF[i] <- as.numeric(ad_values[1])
          }
          if (length(ad_values) >= 2) {
            processed_data$AD_ALT[i] <- as.numeric(ad_values[2])
          }
        }
      }
      
      # Safely convert DP to numeric and calculate VAF
      processed_data$DP <- as.numeric(processed_data$DP)
      
      # Calculate VAF safely
      processed_data$VAF <- NA
      for (i in 1:nrow(processed_data)) {
        if (!is.na(processed_data$AD_REF[i]) && !is.na(processed_data$AD_ALT[i])) {
          total <- processed_data$AD_REF[i] + processed_data$AD_ALT[i]
          if (total > 0) {
            processed_data$VAF[i] <- processed_data$AD_ALT[i] / total
          }
        }
      }
      
      # Apply filtering separately without relying on dplyr pipes
      filtered_data <- processed_data[
        !is.na(processed_data$GT) & 
        processed_data$GT %in% c("0/1", "1/1", "1|0", "0|1", "1|1") &
        !is.na(processed_data$AD_ALT) & 
        processed_data$AD_ALT >= 5 &
        !is.na(processed_data$VAF) & 
        processed_data$VAF >= 0.2 &
        !is.na(processed_data$DP) & 
        processed_data$DP >= 10,
      ]
      
      # Only write CSV if there are variants left after filtering
      if (nrow(filtered_data) > 0) {
        write.csv(filtered_data, file=output_csv_file, row.names=FALSE, quote=TRUE)
        cat("Converted:", vcf_file, "→", output_csv_file, "\n")
        return(TRUE)
      } else {
        cat("Skipped (no significant variants):", vcf_file, "\n")
        return(FALSE)
      }
    }, error = function(e) {
      cat("ERROR processing", vcf_file, ":", conditionMessage(e), "\n")
      return(FALSE)
    })
  }
  
  # Find all matching directories based on the provided pattern
  ngs_dirs <- list.dirs(base_dir, recursive=TRUE)
  ngs_dirs <- grep(dir_pattern, basename(ngs_dirs), value=TRUE)
  ngs_dirs <- file.path(base_dir, ngs_dirs)
  
  # Find all matching VCF files
  file_pattern <- "\\.vcf$"
  vcf_files <- c()
  for (dir in ngs_dirs) {
    matching_files <- list.files(
      dir, 
      pattern=file_pattern,
      full.names=TRUE,
      recursive=FALSE
    )
    vcf_files <- c(vcf_files, matching_files)
  }
  
  # Create output directory if it doesn't exist
  output_dir <- file.path(base_dir, "csv_exports")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive=TRUE)
  }
  
  # Process each file
  results <- list()
  cat("Found", length(vcf_files), "files to process\n")
  
  for (vcf_file in vcf_files) {
    sample_name <- basename(dirname(vcf_file))
    output_file <- file.path(output_dir, paste0(sample_name, "_variants.csv"))
    
    success <- convert_single_file(vcf_file, output_file)
    results[[vcf_file]] <- success
  }
  
  # Combine all CSVs into one master file
  if (sum(unlist(results)) > 0) {
    combined_file <- file.path(output_dir, "combined_variants.csv")
    
    # Get list of successfully created CSV files
    csv_files <- list.files(
      output_dir,
      pattern="_variants\\.csv$",
      full.names=TRUE
    )
    
    # Read and combine all CSVs
    combined_data <- NULL
    for (csv_file in csv_files) {
      # Skip the combined file if it already exists
      if (csv_file == combined_file) next
      
      data <- read.csv(csv_file)
      
      # If this is the first CSV, set combined_data to this data
      if (is.null(combined_data)) {
        combined_data <- data
      } else {
        # Check if column names match
        missing_cols <- setdiff(names(combined_data), names(data))
        if (length(missing_cols) > 0) {
          # Add missing columns with NA values
          for (col in missing_cols) {
            data[[col]] <- NA
          }
        }
        
        # Check for columns in data that aren't in combined_data
        extra_cols <- setdiff(names(data), names(combined_data))
        if (length(extra_cols) > 0) {
          # Add these columns to combined_data with NA values
          for (col in extra_cols) {
            combined_data[[col]] <- NA
          }
        }
        
        # Now columns should match - combine using rbind
        combined_data <- rbind(combined_data, data[names(combined_data)])
      }
    }
    
    # Write combined file
    if (!is.null(combined_data) && nrow(combined_data) > 0) {
      write.csv(combined_data, file=combined_file, row.names=FALSE)
      cat("Created combined variant file:", combined_file, "\n")
    }
  }
  
  # Return the path to the combined file or NULL if none was created
  if (sum(unlist(results)) > 0) {
    return(file.path(output_dir, "combined_variants.csv"))
  } else {
    return(NULL)
  }
}

# Combine all csf snv csvs
combine_all_csvs <- function(csv_dir, pattern = "^CSF", output_file = NULL) {
  # Generate default output filename based on pattern if not provided
  if (is.null(output_file)) {
    output_file <- paste0(gsub("\\^", "", pattern), "_combined_data.csv")
  }
  
  # Construct the file pattern to match specified prefix and _variants.csv suffix
  file_pattern <- paste0(pattern, ".*_variants\\.csv$")
  
  # List all matching CSV files in the directory
  csv_files <- list.files(csv_dir, pattern = file_pattern, full.names = TRUE)
  
  if (length(csv_files) == 0) {
    stop(paste0("No CSV files found matching pattern '", file_pattern, "' in the specified directory"))
  }
  
  # Initialize an empty list to store all data frames
  all_data <- list()
  
  # Read each CSV and add sample name
  for (file in csv_files) {
    tryCatch({
      # Read the CSV file
      data <- read.csv(file, stringsAsFactors = FALSE)
      
      # Extract sample name from filename
      sample_name <- sub("_variants\\.csv$", "", basename(file))
      
      # Check if a sample_name column already exists
      if (!"sample_name" %in% names(data)) {
        # Add sample name column
        data$sample_name <- sample_name
      }
      
      # Add to the list
      all_data[[length(all_data) + 1]] <- data
      
      cat("Processed:", file, "\n")
    }, error = function(e) {
      cat("ERROR processing", file, ":", conditionMessage(e), "\n")
    })
  }
  
  # Combine all data frames
  if (length(all_data) == 0) {
    stop("No data could be processed from the CSV files")
  }
  
  # Use rbindlist from data.table to handle different columns
  library(data.table)
  combined_data <- rbindlist(all_data, fill = TRUE)
  
  # Write combined data to a new CSV file
  output_path <- file.path(csv_dir, output_file)
  write.csv(combined_data, file = output_path, row.names = FALSE)
  
  cat("\nCombined data written to:", output_path, "\n")
  cat("Total rows:", nrow(combined_data), "\n")
  cat("Total files combined:", length(all_data), "\n")
  
  return(invisible(combined_data))
}

# Function that creates CSF BED files when possible, NGS files only when needed
ngs_to_bed <- function(input_file, output_dir = "bed", pairs_taps) {
  # Load packages
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Simple log function
  log <- function(msg) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  
  # 1. Check required columns in pairs_taps
  if (!"NGS" %in% names(pairs_taps)) stop("pairs_taps must have an 'NGS' column")
  if (!"pair" %in% names(pairs_taps)) stop("pairs_taps must have a 'pair' column")
  
  # 2. Create mapping from NGS to CSF IDs
  log("Creating NGS to CSF mapping...")
  ngs_to_csf <- list()
  for (i in 1:nrow(pairs_taps)) {
    ngs_id <- pairs_taps$NGS[i]
    csf_id <- pairs_taps$pair[i]
    if (!is.na(ngs_id) && !is.na(csf_id)) {
      ngs_to_csf[[as.character(ngs_id)]] <- as.character(csf_id)
    }
  }
  log(paste("Created mapping for", length(ngs_to_csf), "NGS IDs to CSF IDs"))
  
  # 3. Read variants with proper settings
  log(paste("Reading variants:", input_file))
  variants <- tryCatch({
    read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) {
    # Try reading it manually
    lines <- readLines(input_file)
    header <- strsplit(lines[1], ",")[[1]]
    
    # Create empty data frame with header names
    result <- data.frame(matrix(ncol = length(header), nrow = 0))
    colnames(result) <- header
    
    # Parse each line
    for (i in 2:length(lines)) {
      fields <- unlist(strsplit(lines[i], ","))
      if (length(fields) >= length(header)) {
        # Take only the fields we need
        fields <- fields[1:length(header)]
        result[nrow(result) + 1, ] <- fields
      }
    }
    return(result)
  })
  
  # 4. Extract required columns
  if (!"CHROM" %in% names(variants)) stop("CHROM column not found")
  if (!"POS" %in% names(variants)) stop("POS column not found")
  
  # Find NGS ID column
  ngs_col <- NULL
  for (col in c("Test.Number", "Test.Number_datadump", "Test.Number_selected")) {
    if (col %in% names(variants)) {
      ngs_col <- col
      break
    }
  }
  if (is.null(ngs_col)) stop("No Test.Number column found")
  
  # 5. Clean variant data
  log("Processing variants...")
  variants$CHROM <- as.character(variants$CHROM)
  variants$CHROM <- ifelse(grepl("^chr", variants$CHROM), variants$CHROM, paste0("chr", variants$CHROM))
  variants$POS <- as.numeric(as.character(variants$POS))
  variants$NGS_ID <- as.character(variants[[ngs_col]])
  
  # Extract NGS IDs if they are in a complex format
  for (i in 1:nrow(variants)) {
    if (grepl("NGS-[0-9]+-[0-9]+", variants$NGS_ID[i])) {
      match <- regexpr("NGS-[0-9]+-[0-9]+", variants$NGS_ID[i])
      if (match > 0) {
        variants$NGS_ID[i] <- regmatches(variants$NGS_ID[i], match)
      }
    }
  }
  
  # Remove rows with missing data
  variants <- variants[!is.na(variants$CHROM) & !is.na(variants$POS) & !is.na(variants$NGS_ID), ]
  
  # 6. Create GRanges and perform liftOver
  log("Creating genomic ranges and performing liftOver...")
  gr <- GRanges(
    seqnames = variants$CHROM,
    ranges = IRanges(start = variants$POS, end = variants$POS),
    NGS_ID = variants$NGS_ID
  )
  
  chain_file <- "/gpfs/home/freite01/DB/references/chain_files/hg19ToHg38.over.chain"
  chain <- import.chain(chain_file)
  gr_hg38 <- unlist(liftOver(gr, chain))
  log(paste("LiftOver complete with", length(gr_hg38), "variants"))
  
  # 7. Create BED files - CSF when available, NGS only when no CSF mapping exists
  unique_ngs_ids <- unique(gr_hg38$NGS_ID)
  log(paste("Found", length(unique_ngs_ids), "unique NGS IDs"))
  
  csf_files_created <- character()
  ngs_files_created <- character()
  skipped_files <- character() # Track skipped files due to already existing
  
  for (ngs_id in unique_ngs_ids) {
    # Get variants for this NGS ID
    variants_subset <- gr_hg38[gr_hg38$NGS_ID == ngs_id]
    
    # Create BED data
    bed_data <- data.frame(
      chrom = seqnames(variants_subset),
      start = start(variants_subset) - 1,  # BED is 0-based
      end = end(variants_subset)
    )
    
    # Check if we have a CSF mapping
    if (ngs_id %in% names(ngs_to_csf)) {
      # Create CSF version
      csf_id <- ngs_to_csf[[ngs_id]]
      csf_file <- file.path(output_dir, paste0(csf_id, "_hg38.bed"))
      
      # Check if file already exists - don't overwrite
      if (file.exists(csf_file)) {
        log(paste("Skipping", csf_file, "- file already exists"))
        skipped_files <- c(skipped_files, csf_file)
      } else {
        write.table(
          bed_data, 
          file = csf_file, 
          quote = FALSE, 
          sep = "\t", 
          row.names = FALSE, 
          col.names = FALSE
        )
        log(paste("Created", csf_file, "with", nrow(bed_data), "variants (from NGS ID", ngs_id, ")"))
        csf_files_created <- c(csf_files_created, csf_file)
      }
    } else {
      # No CSF mapping, create NGS version
      ngs_file <- file.path(output_dir, paste0(ngs_id, "_hg38.bed"))
      
      # Check if file already exists - don't overwrite
      if (file.exists(ngs_file)) {
        log(paste("Skipping", ngs_file, "- file already exists"))
        skipped_files <- c(skipped_files, ngs_file)
      } else {
        write.table(
          bed_data, 
          file = ngs_file, 
          quote = FALSE, 
          sep = "\t", 
          row.names = FALSE, 
          col.names = FALSE
        )
        log(paste("Created", ngs_file, "with", nrow(bed_data), "variants (no CSF mapping available)"))
        ngs_files_created <- c(ngs_files_created, ngs_file)
      }
    }
  }
  
  # 8. Update the pairs_taps table (even with existing files)
  if (!"bed" %in% names(pairs_taps)) {
    pairs_taps[, bed := NA_character_]
  }
  
  # First try to match by CSF ID - including both newly created and existing files
  update_count <- 0
  for (i in 1:nrow(pairs_taps)) {
    if (!is.na(pairs_taps$pair[i])) {
      csf_id <- pairs_taps$pair[i]
      bed_file <- file.path(output_dir, paste0(csf_id, "_hg38.bed"))
      if (file.exists(bed_file)) {
        pairs_taps[i, bed := bed_file]
        update_count <- update_count + 1
      }
    }
  }
  
  # Then try by NGS ID for any remaining rows - including both newly created and existing files
  for (i in 1:nrow(pairs_taps)) {
    if (is.na(pairs_taps$bed[i]) && !is.na(pairs_taps$NGS[i])) {
      ngs_id <- pairs_taps$NGS[i]
      bed_file <- file.path(output_dir, paste0(ngs_id, "_hg38.bed"))
      if (file.exists(bed_file)) {
        pairs_taps[i, bed := bed_file]
        update_count <- update_count + 1
      }
    }
  }
  
  log(paste("Updated", update_count, "rows in pairs_taps with bed paths"))
  
  # 9. Clean up any invalid files (not NGS or CSF format) - but only ones we just created
  all_files <- c(csf_files_created, ngs_files_created)  # Only check files we created in this run
  valid_pattern <- paste0("^", output_dir, "/(NGS-[0-9]+-[0-9]+|CSF-[0-9]+-[0-9]+)_hg38.bed$")
  valid_files <- grep(valid_pattern, all_files, value = TRUE)
  invalid_files <- setdiff(all_files, valid_files)
  
  if (length(invalid_files) > 0) {
    log(paste("Removing", length(invalid_files), "invalid BED files"))
    file.remove(invalid_files)
  }
  
  log(paste("Final result:", length(csf_files_created), "CSF BED files and", 
            length(ngs_files_created), "NGS BED files created (no duplication)"))
  log(paste("Skipped", length(skipped_files), "files that already existed"))
  
  return(list(
    csf_files = csf_files_created,
    ngs_files = ngs_files_created,
    skipped_files = skipped_files,
    pairs_taps = pairs_taps
  ))
}

## Filtering snpeff results
filter_snpeff_results <- function(snv,
                                min_dp = 15,
                                max_dp = 150,          
                                min_tlod = 6.0,
                                min_alt_reads = 5, 
                                keep_high = TRUE,
                                keep_high_moderate = TRUE,
                                keep_upstream_modifiers = TRUE,
                                keep_all_impacts = FALSE, # Add this new parameter
                                keep_one_per_variant = TRUE,
                                remove_duplicates = TRUE,
                                return_granges = TRUE) {
  
  snv_df <- as.data.frame(snv)

  cat("Filtering with columns:", paste(head(colnames(snv_df), 10), collapse=", "), "...\n")

  original_count <- nrow(snv_df)
  
  filtered_snv <- snv_df
  
  # Filter by depth
  if ("DP" %in% colnames(snv_df)) {
    filtered_snv <- filtered_snv %>%
      filter(DP >= min_dp & DP <= max_dp)
    cat(sprintf("Filtered by depth: %d/%d variants remain\n", 
                nrow(filtered_snv), original_count))
  } else {
    cat("Warning: 'DP' column not found. Skipping depth filter.\n")
  }
  
  # Filter by TLOD
  if ("TLOD" %in% colnames(filtered_snv)) {
    # Extract numeric values for TLOD
    if (is.list(filtered_snv$TLOD)) {
      # For list column
      filtered_snv$tlod_numeric <- sapply(filtered_snv$TLOD, function(x) {
        if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA_real_)
        if (is.numeric(x)) return(x[1])
        tryCatch(as.numeric(as.character(x)[1]), error = function(e) NA_real_)
      })
    } else {
      # For regular column
      filtered_snv$tlod_numeric <- filtered_snv$TLOD
    }
    
    filtered_snv <- filtered_snv %>%
      filter(tlod_numeric >= min_tlod)
    
    cat(sprintf("Filtered by TLOD: %d/%d variants remain\n", 
                nrow(filtered_snv), original_count))
  } else {
    cat("Warning: 'TLOD' column not found. Skipping TLOD filter.\n")
  }
  
  # Filter by alt read count
  if ("AS_SB_TABLE" %in% colnames(filtered_snv)) {
    # Extract the total alt read count from AS_SB_TABLE
    filtered_snv$alt_reads <- as.numeric(sapply(filtered_snv$AS_SB_TABLE, function(x) {
      as.numeric(x)
    }))
    
    filtered_snv <- filtered_snv %>%
      filter(alt_reads >= min_alt_reads)
    
    cat(sprintf("Filtered by alt read count: %d/%d variants remain\n", 
                nrow(filtered_snv), original_count))
  } else {
    cat("Warning: 'AS_SB_TABLE' column not found. Skipping alt read count filter.\n")
  }
  
  # Apply impact and annotation filter (general rule approach)
  if ("impact" %in% colnames(filtered_snv) && !keep_all_impacts) {
    # Only apply impact filtering if keep_all_impacts is FALSE
    # Start with an empty result
    impact_filtered <- data.frame()
    
    # Keep HIGH and MODERATE impact variants if requested
    if (keep_high_moderate) {
      high_moderate <- filtered_snv %>%
        filter(impact %in% c("HIGH", "MODERATE"))
      
      impact_filtered <- bind_rows(impact_filtered, high_moderate)
    }

    # Keep HIGH
    if (keep_high) {
      high <- filtered_snv %>%
        filter(impact %in% c("HIGH"))
      
      impact_filtered <- bind_rows(impact_filtered, high)
    }
    
    # Keep upstream gene variants with MODIFIER impact if requested
    if (keep_upstream_modifiers && "annotation" %in% colnames(filtered_snv)) {
      upstream_variants <- filtered_snv %>%
        filter(
          impact == "MODIFIER" & 
          grepl("upstream_gene_variant", annotation)
        )
      
      impact_filtered <- bind_rows(impact_filtered, upstream_variants)
    }
    
    # Update filtered_snv with the results
    filtered_snv <- impact_filtered
    
    cat(sprintf("Filtered by impact and annotation: %d/%d variants remain\n", 
                nrow(filtered_snv), original_count))
  } else if ("impact" %in% colnames(filtered_snv) && keep_all_impacts) {
    cat("Keeping all impacts as requested, no impact filtering applied\n")
  } else {
    cat("Warning: 'impact' column not found. Skipping impact filter.\n")
  }
  
  # Remove duplicates if requested
  if (remove_duplicates) {
    before_dedup <- nrow(filtered_snv)
    filtered_snv <- filtered_snv %>%
      distinct(seqnames, start, end, REF, ALT, .keep_all = TRUE)
    
    cat(sprintf("After removing duplicates: %d/%d variants remain\n", 
                nrow(filtered_snv), before_dedup))
  }
  
  # Keep one annotation per variant if requested
  if (keep_one_per_variant && "impact" %in% colnames(filtered_snv)) {
    # Set impact priority
    impact_order <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
    filtered_snv$impact_rank <- match(filtered_snv$impact, impact_order)
    
    # Count before grouping
    before_count <- nrow(filtered_snv)
    
    # Group and select highest impact annotation
    filtered_snv <- filtered_snv %>%
      group_by(seqnames, start, end, REF, ALT) %>%
      arrange(impact_rank) %>%
      slice(1) %>%
      ungroup() %>%
      select(-impact_rank)
    
    cat(sprintf("After keeping one annotation per variant: %d/%d variants remain\n", 
                nrow(filtered_snv), before_count))
  }
  
  # Clean up temporary columns
  filtered_snv <- filtered_snv %>%
    select(-any_of(c("tlod_numeric", "alt_reads")))
  
  # Convert back to GRanges if requested
  if (return_granges) {
    # Convert to GRanges object
    gr <- makeGRangesFromDataFrame(filtered_snv, 
                                  keep.extra.columns = TRUE, 
                                  ignore.strand = FALSE,
                                  seqnames.field = "seqnames",
                                  start.field = "start",
                                  end.field = "end")
    
    # Make sure to set the genome information for the GRanges object
    genome(gr) <- genome(snv)
    
    return(gr)
  } else {
    return(filtered_snv)
  }
}
