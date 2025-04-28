#!/usr/bin/env Rscript

#' Genomic Data Plot Generation for Copy Number Variation
#' 
#' This script analyzes and visualizes genomic copy number variation data, 
#' particularly focusing on segmentation statistics from RDS files. 
#' It creates gene plots with customizable parameters.
#' 
#' @param seg_file Path to the segmentation RDS file
#' @param cov_file Path to the coverage RDS file
#' @param gene_of_interest Gene name to visualize (e.g., "SMARCB1")
#' @return An SVG or PDF plot visualizing the gene CNV data

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
  library(skitools)
})

#' Color mapping function for numeric values
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

#' Main function to generate genomic copy number variation plot
plot_gene_cnv <- function(
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
  # Read input files
  seg = readRDS(seg_file)
  cov = readRDS(cov_file)
  
  # Disable all warnings globally during execution
  oldw <- getOption("warn")
  options(warn = -1)
  
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
  armval = gr.findoverlaps(cov %Q% (!is.na(foreground)), cyto, qcol = names(mcols(cov)), scol = "arm")
  
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
  GENES = gene_of_interest
  
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
  
  # Restore warning level
  options(warn = oldw)
}

#' Function to display SVG plot in different contexts
show_svg <- function(...) {
  plot_gene_cnv(...)
}

#' Cache mechanism for plots
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

# If script is run directly (not sourced), execute example
if (sys.nframe() == 0) {
  # Suppress all warnings for the script execution
  options(warn = -1)
  
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
    
    # Run without showing warnings
    plot_gene_cnv(seg_file, cov_file, gene, output_file = output_file, pdf_filename = pdf_filename)
    
  } else {
    cat("Usage: Rscript genomic_plot.R <seg_file.rds> <cov_file.rds> <gene_name> [output_file.svg|pdf]\n")
    cat("Example: Rscript genomic_plot.R seg.rds cov.rds SMARCB1 output.svg\n")
    cat("Example: Rscript genomic_plot.R seg.rds cov.rds SMARCB1 output.pdf\n")
  }
  
  # Restore warnings
  options(warn = 0)
}
