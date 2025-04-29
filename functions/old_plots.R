#!/usr/bin/env Rscript

#' Genomic Data Plot Generation
#' 
#' This script analyzes and visualizes genomic data, particularly focusing on 
#' segmentation statistics from RDS files. It creates gene plots with customizable
#' parameters.
#' 
#' @param seg_file Path to the segmentation RDS file
#' @param cov_file Path to the coverage RDS file
#' @param gene_of_interest Gene name to visualize (e.g., "SMARCB1")
#' @return An SVG plot visualizing the gene data
#' 
#' @author Converted from Jupyter notebook code
#' @date April 24, 2025

# Load required libraries
suppressPackageStartupMessages({
  library(JaBbA)
  library(Flow)
  library(gUtils)
  library(gTrack)
  library(RKernel)
  library(dplyr)
  library(data.table)
  library(rtracklayer)
  library(GenomeInfoDb)
  library(khtools)
  library(signal)
  library(parallel)
})

#' Color mapping function for numeric values
#' 
#' @param x Numeric vector to color
#' @param colors Color palette to use
#' @param alpha Alpha transparency (0-1)
#' @param capvals Boolean whether to cap values
#' @param maxVal Maximum value for capping
#' @param minVal Minimum value for capping
#' @return Vector of color codes
numeric2color <- function(x, colors, alpha, capvals = TRUE, maxVal = 8, minVal = 0) {
  maxx = max(x, na.rm = TRUE)
  minx = min(x, na.rm = TRUE)
  
  if (identical(capvals, TRUE)) {
    minx = max(minVal, min(x, na.rm = TRUE), na.rm = TRUE)
    maxx = min(maxVal, max(x, na.rm = TRUE), na.rm = TRUE)
  }
  
  x_0to1 = (x - minx) / (maxx - minx)
  
  if (identical(capvals, TRUE)) {
    x_0to1 = pmax(x_0to1, 0)
    x_0to1 = pmin(x_0to1, 1)
  }
  
  cols = colorRamp(colors)(x_0to1)
  outCols = rep(NA_character_, length(x_0to1))
  notNaIdx = which(!is.na(x_0to1))
  goodCols = cols[notNaIdx,]
  
  if (length(notNaIdx) > 0) {
    goodCols = rgb(red = goodCols[,1], green = goodCols[,2], blue = goodCols[,3], alpha = alpha, maxColorValue = 255)
    outCols[notNaIdx] = goodCols
  }
  
  return(outCols)
}

#' Main function to generate genomic plot
#' 
#' @param seg_file Path to segmentation RDS file
#' @param cov_file Path to coverage RDS file
#' @param gene_of_interest Gene to visualize (character string)
#' @param output_width Plot width (default: 18)
#' @param output_file Output file path (default: NULL for display)
#' @param cytoband_file Path to cytoband file (default path provided)
#' @param gencode_file Path to gencode annotation file (default path provided)
#' @param gencode_cache_dir Directory for gencode cache (default path provided)
#' @param pdf_filename Optional PDF filename to save plot (default: NULL)
#' @return Plot object (if output_file is NULL) or saves to file
plot_gene_coverage <- function(
  seg_file,
  cov_file,
  gene_of_interest,
  output_width = 18,
  output_file = NULL,
  cytoband_file = "/gpfs/data/imielinskilab/DB/UCSC/hg38.cytoband.txt",
  gencode_file = "~/DB/GENCODE/hg38/v29/gencode.v29.annotation.gff3",
  gencode_cache_dir = "~/DB/GENCODE/hg38/v29/",
  pdf_filename = NULL
) {
  # Suppress all warnings
  options(warn = -1)
  # Read input files
  seg = readRDS(seg_file)
  cov = readRDS(cov_file)
  
  # Set sequence style
  GenomeInfoDb::seqlevelsStyle(seg) = "NCBI"
  GenomeInfoDb::seqlevelsStyle(cov) = "NCBI"
  
  # Calculate segmentation statistics
  gr_segstats = JaBbA:::segstats(target = seg, signal = cov, field = "foreground")
  
  # Process cytoband information
  cyto = fread(cytoband_file)
  names(cyto) = c("seqnames", "start", "end", "band", "stain")
  isZeroStart = any(cyto[, length(intersect(start, end)), by = seqnames]$V1 > 0) || any(cyto$start == 0)
  if (isZeroStart) cyto$start = cyto$start + 1
  cyto = gUtils::dt2gr(cyto)
  GenomeInfoDb::seqlevelsStyle(cyto) = "NCBI"
  cyto$chrom_name = as.character(seqnames(cyto))
  cyto$chromband = paste(cyto$chrom_name, cyto$band, sep = "")
  cyto$arm = gsub("[0-9.]+", "", cyto$band)
  
  # Find overlaps between coverage and cytoband
  # Just use the original code exactly as provided
  armval = gr.findoverlaps(cov %Q% (!is.na(foreground)), cyto, qcol = names(mcols(cov)), scol = "arm")
  
  # Recalculate segmentation statistics
  gr_segstats = JaBbA:::segstats(target = seg, signal = cov, field = "foreground")
  
  # Process arm values
  armvalSplit = armval %>% split(paste(seqnames(.), .$arm, sep = "__"))
  
  # Suppress warnings for smoothing operations
  suppressWarnings(
    smoothSplit <- mclapply(armvalSplit, function(x) {
      if (length(x$foreground) < 3) return(x$foreground)
      signal::sgolayfilt(x$foreground)
    }, mc.cores = 8)
  )
  
  armsmooth = unlist(armvalSplit)
  armsmooth$foreground_smooth = unlist(smoothSplit)
  armsmooth = unname(sort(armsmooth))
  
  # Calculate smooth segmentation statistics
  gr_segstats_smooth = JaBbA:::segstats(target = seg, signal = armsmooth, field = "foreground_smooth")
  
  # Calculate PP fit
  pp_fit = JaBbA:::ppgrid(gr_segstats)
  bestpp_fit = pp_fit[1,]
  purity = bestpp_fit$purity
  ploidy = bestpp_fit$ploidy
  gamma = bestpp_fit$gamma
  beta = bestpp_fit$beta
  
  # Calculate PP fit for smooth data
  pp_fit = JaBbA:::ppgrid(gr_segstats_smooth)
  bestpp_fit = pp_fit[1,]
  purity = bestpp_fit$purity
  ploidy = bestpp_fit$ploidy
  gamma = bestpp_fit$gamma
  beta = bestpp_fit$beta
  
  # Calculate CN values
  gr_segstats$cn = skitools::rel2abs(gr_segstats, purity = purity, ploidy = ploidy, field = "mean")
  gr_segstats_smooth$cn = skitools::rel2abs(gr_segstats_smooth, purity = purity, ploidy = ploidy, field = "mean")
  
  # Import gencode data
  gencode = rtracklayer::import.gff3(gencode_file)
  GenomeInfoDb::seqlevelsStyle(gencode) = "NCBI"
  
  # Create gencode track
  gt_gc = gTrack::track.gencode(
    cached = TRUE,
    gene.collapse = TRUE,
    cached.dir = gencode_cache_dir,
    cached.path = paste0(gencode_cache_dir, "gencode.composite.rds"),
    cached.path.collapsed = paste0(gencode_cache_dir, "gencode.composite.collapsed.rds")
  )
  
  # Filter genes
  genes = gencode %Q% (type == "gene") %Q% (gene_type == "protein_coding")
  genes = gr.val(genes, gr_segstats_smooth, "cn", na.rm = TRUE)
  
  # Set amplification threshold
  genes$amp = FALSE
  genes$amp = (genes$cn / ploidy) > 2
  
  genes$homdel = genes$cn <= 0.5
  if (ploidy <= 2.7) {
    genes$amp = genes$cn >= 5
  } else {
    genes$amp = genes$cn >= 8
  }
  
  # Calculate CN for armsmooth
  armsmooth$cn = skitools::rel2abs(armsmooth, purity = purity, ploidy = ploidy, field = "foreground_smooth")
  
  # Set color and width
  gr_segstats_smooth$ywid = 2
  gr_segstats_smooth$col = numeric2color(
    gr_segstats_smooth$cn,
    c("darkblue", "grey", "red"),
    capvals = TRUE
  )
  
  armsmooth$col = numeric2color(
    armsmooth$cn,
    c("darkblue", "grey", "red"),
    capvals = TRUE
  )
  
  # Create tracks
  gt_signal = gTrack::gTrack(armsmooth, "cn", circles = TRUE, lwd.border = 0.25, y0 = -1)
  gt_segs = gTrack::gTrack(gr_segstats_smooth, "cn", y0 = -1)
  
  # Filter to manageable regions
  gr_segstats_smooth = gr_segstats_smooth[width(gr_segstats_smooth) < 10000]
  
  # Create the plot
  GENES = gene_of_interest # e.g., "SMARCB1"
  
  if (!is.null(pdf_filename)) {
    # Use skitools::ppdf for PDF output
    skitools::ppdf(
      filename = pdf_filename,
      width = output_width
    )
    
    par(oma = c(0, 0, 0, 0), mai = c(0, 2, 1, 3))
    plot(
      c(
        gt_gc %>% {x = .; x@data[[1]] = x@data[[1]][GENES]; x},
        gt_signal,
        gt_segs
      ),
      win = (genes %Q% (gene_name %in% GENES) %>% khtools::gr.noval()) + 1e6,
      legend.params = list(plot = FALSE),
      y.quantile = 0.95,
      xaxis.suffix = "Mb",
      xaxis.unit = 1e6,
      xaxis.round = 1,
      xaxis.width = FALSE,
      ylab.las = 1,
      ylab.adj = 0
    )
    dev.off()
    message(paste("PDF saved to:", pdf_filename))
  } else if (!is.null(output_file)) {
    svg(output_file, width = output_width, height = output_width * 0.6)
    
    par(oma = c(0, 0, 0, 0), mai = c(0, 2, 1, 3))
    plot(
      c(
        gt_gc %>% {x = .; x@data[[1]] = x@data[[1]][GENES]; x},
        gt_signal,
        gt_segs
      ),
      win = (genes %Q% (gene_name %in% GENES) %>% khtools::gr.noval()) + 1e6,
      legend.params = list(plot = FALSE),
      y.quantile = 0.95,
      xaxis.suffix = "Mb",
      xaxis.unit = 1e6,
      xaxis.round = 1,
      xaxis.width = FALSE,
      ylab.las = 1,
      ylab.adj = 0
    )
    
    dev.off()
    message(paste("SVG saved to:", output_file))
  } else {
    show_svg_plot <- function() {
      par(oma = c(0, 0, 0, 0), mai = c(0, 2, 1, 3))
      plot(
        c(
          gt_gc %>% {x = .; x@data[[1]] = x@data[[1]][GENES]; x},
          gt_signal,
          gt_segs
        ),
        win = (genes %Q% (gene_name %in% GENES) %>% khtools::gr.noval()) + 1e6,
        legend.params = list(plot = FALSE),
        y.quantile = 0.95,
        xaxis.suffix = "Mb",
        xaxis.unit = 1e6,
        xaxis.round = 1,
        xaxis.width = FALSE,
        ylab.las = 1,
        ylab.adj = 0
      )
    }
    
    show_svg_plot()
  }
  
  # Restore warning level
  options(warn = 0)
}

#' Function to display SVG plot in different contexts
#' 
#' @param ... Arguments passed to plot_gene_coverage function
#' @return SVG plot object
show_svg <- function(...) {
  plot_gene_coverage(...)
}

#' Cache mechanism for plots
#' 
#' @param file Path to cache file
#' @param FUN Function to execute if cache doesn't exist
#' @param ... Arguments passed to FUN
#' @return Result of FUN or cached result
cache_plot <- function(file, FUN, ...) {
  if (file.exists(file)) {
    message(paste("Loading cached plot from:", file))
    return(readRDS(file))
  } else {
    result <- FUN(...)
    saveRDS(result, file)
    return(result)
  }
}

# Add a function to check package availability and install if needed
check_and_install_packages <- function() {
  required_packages <- c(
    "JaBbA", "Flow", "gUtils", "gTrack", "RKernel", "dplyr", 
    "data.table", "rtracklayer", "GenomeInfoDb", "khtools", 
    "signal", "parallel"
  )
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    
    # First try Bioconductor for bioinformatics packages
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    
    for (pkg in missing_packages) {
      tryCatch({
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      }, error = function(e) {
        message("Failed to install via BiocManager, trying CRAN: ", pkg)
        install.packages(pkg)
      })
    }
  }
}

# If script is run directly (not sourced), execute example
if (sys.nframe() == 0) {
  # Example usage when script is run directly
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) >= 3) {
    seg_file <- args[1]
    cov_file <- args[2]
    gene <- args[3]
    
    output_file <- NULL
    pdf_filename <- NULL
    
    if (length(args) >= 4) {
      if (grepl("\\.pdf$", args[4], ignore.case = TRUE)) {
        pdf_filename <- args[4]
      } else {
        output_file <- args[4]
      }
    }
    
    suppressWarnings(
      tryCatch({
        plot_gene_coverage(seg_file, cov_file, gene, output_file = output_file, pdf_filename = pdf_filename)
      }, error = function(e) {
        cat("Error in plot_gene_coverage: ", e$message, "\n")
      })
    )
  } else {
    cat("Usage: Rscript genomic_plot.R <seg_file.rds> <cov_file.rds> <gene_name> [output_file.svg|pdf]\n")
    cat("Example: Rscript genomic_plot.R seg.rds cov.rds SMARCB1 output.svg\n")
    cat("Example: Rscript genomic_plot.R seg.rds cov.rds SMARCB1 output.pdf\n")
  }
}
