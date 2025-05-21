# Comprehensive function to construct, find and update paths with enhanced flexibility
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
  preserve_existing = TRUE,    # Whether to preserve existing entries in target column
  flat_structure = FALSE,      # If TRUE, creates paths like base_dir/pairid_pattern instead of base_dir/pairid/pattern
  id_in_filename = FALSE       # If TRUE, creates paths like base_dir/id/id.pattern
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
  
  # Construct paths for each pair based on parameters
  paths <- sapply(pairs_taps[[id_column]], function(pair) {
    if (is.na(pair)) return(NA)
    
    if (flat_structure) {
      # Flat structure: base_dir/pairid_pattern
      file_path <- file.path(base_dir, paste0(pair, pattern))
    } else if (id_in_filename) {
      # ID in both directory and filename: base_dir/id/id.pattern
      file_path <- file.path(base_dir, pair, paste0(pair, pattern))
    } else {
      # Original behavior: base_dir/pairid/pattern
      file_path <- file.path(base_dir, pair, pattern)
    }
    
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
  
  # If match_directories is TRUE and not using flat_structure, check for existing folders
  if (match_directories && (!flat_structure || id_in_filename)) {
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
        # Construct path with the folder based on selected format
        if (id_in_filename) {
          path <- file.path(base_dir, folder, paste0(folder, pattern))
        } else {
          path <- file.path(base_dir, folder, pattern)
        }
        
        # Check if file exists
        if (file.exists(path) || create_if_missing) {
          pairs_taps[matching_rows, (target_column) := path]
        }
      }
    }
  }
  # For flat_structure, we can optionally check for existing files matching our pattern
  else if (match_directories && flat_structure) {
    # Get list of files in the base directory
    if (!dir.exists(base_dir)) {
      if (create_if_missing) {
        dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
      } else {
        warning(paste("Base directory does not exist:", base_dir))
        return(pairs_taps)
      }
    }
    
    # List all files in the directory
    existing_files <- list.files(base_dir, full.names = FALSE)
    
    # Check for each pair if a matching file exists
    for (i in 1:nrow(pairs_taps)) {
      pair_id <- pairs_taps[[id_column]][i]
      if (!is.na(pair_id)) {
        expected_file <- paste0(pair_id, pattern)
        if (expected_file %in% existing_files) {
          pairs_taps[i, (target_column) := file.path(base_dir, expected_file)]
        }
      }
    }
  }
  
  # Create files/directories if requested and they don't exist
  if (create_if_missing) {
    for (i in 1:nrow(pairs_taps)) {
      path <- pairs_taps[[target_column]][i]
      if (!is.na(path) && !file.exists(path)) {
        # Create directory if it doesn't exist
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        
        # Create empty file (touch)
        file.create(path)
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
