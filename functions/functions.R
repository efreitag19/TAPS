# Comprehensive function to construct, find and update paths in pairs_taps
construct_paths <- function(
  pairs_taps,                  # The data.table to update
  base_dir,                    # Base directory path
  pattern,                     # File pattern to search for or append (NO DEFAULT)
  id_column = "pair",          # Column containing identifiers
  target_column = NULL,        # Column to update (if NULL, derived from pattern)
  clean_invalid = TRUE,        # Whether to set invalid paths to NA
  recursive = FALSE,           # Whether to search directories recursively
  match_directories = TRUE,    # Whether to match directory names to pair IDs
  create_if_missing = FALSE,   # Whether to create directories/files if missing
  preserve_existing = TRUE     # Whether to preserve existing entries in target column
) {
  # Required packages
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required")
  }
  
  # Ensure pairs_taps is a data.table
  if (!data.table::is.data.table(pairs_taps)) {
    pairs_taps <- data.table::as.data.table(pairs_taps)
  }
  
  # If target_column is NULL, derive it from pattern
  if (is.null(target_column)) {
    # Extract name before extension
    if (grepl("\\.", pattern)) {
      target_column <- sub("\\.[^.]+$", "", pattern)
    } else {
      target_column <- pattern
    }
    message(paste("Using derived target column:", target_column))
  }
  
  # Create target column if it doesn't exist
  if (!target_column %in% names(pairs_taps)) {
    pairs_taps[, (target_column) := NA_character_]
  }
  
  # Save existing valid entries if preserve_existing is TRUE
  existing_entries <- NULL
  if (preserve_existing) {
    # Get indices of rows with non-NA values in target column
    valid_indices <- which(!is.na(pairs_taps[[target_column]]))
    if (length(valid_indices) > 0) {
      existing_entries <- data.table(
        row_index = valid_indices,
        value = pairs_taps[valid_indices, get(target_column)]
      )
      message(sprintf("Preserving %d existing entries in column '%s'", 
                    length(valid_indices), target_column))
    }
  }
  
  # Construct paths for each pair
  paths <- sapply(pairs_taps[[id_column]], function(pair) {
    if (is.na(pair)) return(NA)
    file_path <- file.path(base_dir, pair, pattern)
    return(file_path)
  })
  
  # Update the target column with constructed paths
  pairs_taps[, (target_column) := paths]
  
  # Clean invalid paths if requested
  if (clean_invalid) {
    # Check which paths exist
    valid <- sapply(pairs_taps[[target_column]], file.exists)
    # Set invalid paths to NA
    pairs_taps[!valid, (target_column) := NA]
  }
  
  # If match_directories is TRUE, check for existing folders and update paths
  if (match_directories) {
    # Get list of existing directories
    if (!dir.exists(base_dir)) {
      if (create_if_missing) {
        dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
      } else {
        warning(paste("Base directory does not exist:", base_dir))
        return(pairs_taps)
      }
    }
    
    existing_folders <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
    existing_folders <- existing_folders[existing_folders != ""]
    
    # Loop through folders and find matching rows
    for (folder in existing_folders) {
      matching_rows <- which(pairs_taps[[id_column]] == folder)
      
      if (length(matching_rows) > 0) {
        # Construct path with the folder
        path <- file.path(base_dir, folder, pattern)
        
        # Check if file exists
        if (file.exists(path) || create_if_missing) {
          pairs_taps[matching_rows, (target_column) := path]
        }
      }
    }
  }
  
  # Report results
  valid_count <- sum(!is.na(pairs_taps[[target_column]]))
  total_count <- nrow(pairs_taps)
  message(sprintf("Updated %d/%d rows with valid paths in column '%s'", 
                  valid_count, total_count, target_column))
  
  # Restore existing entries if preserve_existing is TRUE
  if (preserve_existing && !is.null(existing_entries) && nrow(existing_entries) > 0) {
    for (i in 1:nrow(existing_entries)) {
      row_idx <- existing_entries[i, row_index]
      value <- existing_entries[i, value]
      
      # Only restore if the current value is NA
      if (is.na(pairs_taps[row_idx, get(target_column)])) {
        pairs_taps[row_idx, (target_column) := value]
      }
    }
    
    # Report how many entries were preserved
    preserved_count <- sum(!is.na(pairs_taps[existing_entries$row_index, get(target_column)]))
    message(sprintf("Preserved %d existing entries in column '%s'", 
                    preserved_count, target_column))
  }
  
  return(pairs_taps)
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
  for (col in c("Test.Number_datadump", "Test.Number_selected")) {
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
    } else {
      # No CSF mapping, create NGS version
      ngs_file <- file.path(output_dir, paste0(ngs_id, "_hg38.bed"))
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
  
  # 8. Update the pairs_taps table
  if (!"bed" %in% names(pairs_taps)) {
    pairs_taps[, bed := NA_character_]
  }
  
  # First try to match by CSF ID
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
  
  # Then try by NGS ID for any remaining rows
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
  
  # 9. Clean up any invalid files (not NGS or CSF format)
  all_files <- list.files(output_dir, pattern = "_hg38.bed$", full.names = TRUE)
  valid_pattern <- paste0("^", output_dir, "/(NGS-[0-9]+-[0-9]+|CSF-[0-9]+-[0-9]+)_hg38.bed$")
  valid_files <- grep(valid_pattern, all_files, value = TRUE)
  invalid_files <- setdiff(all_files, valid_files)
  
  if (length(invalid_files) > 0) {
    log(paste("Removing", length(invalid_files), "invalid BED files"))
    file.remove(invalid_files)
  }
  
  log(paste("Final result:", length(csf_files_created), "CSF BED files and", 
            length(ngs_files_created), "NGS BED files created (no duplication)"))
  
  return(list(
    csf_files = csf_files_created,
    ngs_files = ngs_files_created,
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
  if ("impact" %in% colnames(filtered_snv)) {
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
