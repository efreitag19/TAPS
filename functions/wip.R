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
