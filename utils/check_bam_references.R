# More robust BAM reference checker
library(Rsamtools)

# Function to get detailed reference information
examine_bam_reference <- function(bam_file) {
  tryCatch({
    # Get BAM header
    header <- scanBamHeader(bam_file)
    
    # Extract reference sequence data
    ref_seqs <- names(header[[1]]$targets)
    ref_lengths <- unname(header[[1]]$targets)
    
    # Get sequence count
    seq_count <- length(ref_seqs)
    
    # Check first few chromosomes and their length
    first_chroms <- head(ref_seqs, 5)
    first_lengths <- head(ref_lengths, 5)
    
    # Check if likely standard reference (chr1, chr2, etc.)
    has_standard_chroms <- any(grepl("^chr[0-9XY]+$", ref_seqs))
    
    # Check if primary chromosomes are present
    has_chr1 <- any(grepl("^chr1$", ref_seqs))
    
    # Try to identify reference type
    ref_type <- "Unknown"
    if (seq_count > 3000) {
      ref_type <- "WMG-like (many contigs)"
    } else if (seq_count < 500 && has_standard_chroms) {
      ref_type <- "Standard reference"
    }
    
    # Specifically check for length of chr1 which can be distinctive
    chr1_len <- NA
    if (has_chr1) {
      chr1_idx <- which(ref_seqs == "chr1")[1]
      chr1_len <- ref_lengths[chr1_idx]
    }
    
    # Return comprehensive info
    return(list(
      file = bam_file,
      seq_count = seq_count,
      ref_type = ref_type,
      first_chroms = first_chroms,
      first_lengths = first_lengths,
      has_standard_chroms = has_standard_chroms,
      chr1_length = chr1_len
    ))
  }, error = function(e) {
    return(list(
      file = bam_file,
      error = as.character(e)
    ))
  })
}

# Function to run on a vector of BAM files
check_bam_files <- function(bam_files) {
  results <- list()
  
  for (i in seq_along(bam_files)) {
    bam_file <- bam_files[i]
    cat(sprintf("[%d/%d] Examining: %s\n", i, length(bam_files), basename(bam_file)))
    
    result <- examine_bam_reference(bam_file)
    results[[i]] <- result
    
    if (!is.null(result$error)) {
      cat("  ERROR:", result$error, "\n")
    } else {
      cat("  Sequences:", result$seq_count, "\n")
      cat("  Reference type:", result$ref_type, "\n")
      cat("  First chromosomes:", paste(result$first_chroms, collapse=", "), "\n")
      cat("  Chr1 length:", result$chr1_length, "\n")
    }
    cat("\n")
  }
  
  return(results)
}

# Main execution function
run_bam_check <- function() {
  # Check if samples$bam_input exists
  if (!exists("samples") || !("bam_input" %in% names(samples))) {
    stop("samples$bam_input not found")
  }
  
  # Get BAM files
  bam_files <- samples$bam_input
  
  # Print total files to check
  cat("Found", length(bam_files), "BAM files to check\n\n")
  
  # Run check
  results <- check_bam_files(bam_files)
  
  # Group by reference type and count
  ref_types <- sapply(results, function(r) if(!is.null(r$ref_type)) r$ref_type else "Error")
  seq_counts <- sapply(results, function(r) if(!is.null(r$seq_count)) r$seq_count else NA)
  
  cat("\n=== SUMMARY ===\n")
  cat("Reference types:\n")
  print(table(ref_types))
  
  cat("\nSequence counts:\n")
  print(table(seq_counts))
  
  # Return results
  return(results)
}

# Call the function
# results <- run_bam_check()
