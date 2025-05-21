#!/usr/bin/env Rscript

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

# Import external functions if needed
import::from("~/scripts/r_sources/jupyter_plot.R", show_svg, cache_plot)

# Set Jupyter plot width globally
options("jupyter.plot.width" = 20)

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

# Main function to generate genomic copy number variation plot
plot_gene_cnv <- function(
  seg_file,
  cov_file,
  gene_of_interest,
  output_width = 18,
  output_file = NULL,
  cytoband_file = "/gpfs/data/imielinskilab/DB/UCSC/hg38.cytoband.txt",
  gencode_file = "~/DB/GENCODE/hg38/v29/gencode.v29.annotation.gff3",
  gencode_cache_dir = "~/DB/GENCODE/hg38/v29/",
  filename = NULL
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
  
  if (!is.null(filename)) {
    # Use skitools::ppdf for PDF output
    skitools::ppdf(
      {
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
      },
      filename = filename,
      width = output_width
    )
    message(paste("PDF saved to:", filename))
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

# We're importing show_svg and cache_plot from external file
# import::from("~/scripts/r_sources/jupyter_plot.R", show_svg, cache_plot)
# So we don't redefine them here

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
    png_filename <- NULL
    
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

# To analyze gene fusions
analyze_fusion_junction <- function(bam_path, 
                                    transcript1 = "ENST00000302326",  # CXXC5
                                    transcript2 = "ENST00000302517",  # MN1
                                    transcript1_exon = 2, # exon 2
                                    transcript1_intron = 1,
                                    transcript2_exon = 1, # exon 1
                                    transcript2_intron = 1,
                                    filename = NULL,
                                    output_format = "svg",
                                    flip_transcripts = FALSE) {
  # Load required libraries exactly as in your code
  library(gGnome)
  library(gUtils)
  library(khtools)
  library(BiocManager)
  library(rtracklayer)
  
  # Import GENCODE annotations exactly as in your code
  gr_v29 = rtracklayer::import.gff3("~/DB/GENCODE/hg38/v29/gencode.v29.annotation.gff3")
  GenomeInfoDb::seqlevelsStyle(gr_v29) = "NCBI"
  gencode = readRDS("~/DB/GENCODE/hg38/v29/gencode.v29.annotation.nochr.rds")
  
  # Load RefSeq lookup data
  refseq_lookup = readRDS("~/DB/GENCODE/refseq_ensembl_lookup.rds")
  data.table::setDT(refseq_lookup)
  
  # Extract exons for each transcript - match your exact approach
  transcript1_exons = gencode %Q% which(grepl(transcript1, transcript_id)) %Q% (type == "exon")
  transcript2_exons = gencode %Q% which(grepl(transcript2, transcript_id)) %Q% (type == "exon")
  
  # Process all exons
  exons = gencode %Q% (type == "exon")
  exons_by_tx <- split(exons, exons$transcript_id)
  
  # Get exon and intron ranges
  exon_ranges = range(exons_by_tx)
  intron_ranges = GenomicRanges::setdiff(exon_ranges, exons_by_tx)
  
  # Process intron ranges exactly as in your code
  indices <- rep(seq_along(intron_ranges), lengths(intron_ranges))
  tx_ids <- names(intron_ranges)[indices]
  
  intron_ranges_clean <- unlist(intron_ranges)
  mcols(intron_ranges_clean)$transcript_id <- tx_ids
  
  # Order introns exactly as in your code
  ordered_introns <- intron_ranges_clean[order(ifelse(strand(intron_ranges_clean) == "+", 1, -1) * start(intron_ranges_clean))]
  
  ordered_introns <- (
    gUtils::gr2dt(ordered_introns) %>%
      .[, intron_number := seq_len(.N), by = transcript_id] %>%
      gUtils::dt2gr()
  )
  
  # Process each transcript's introns - match your exact approach
  transcript2_introns = ordered_introns %Q% (grepl(transcript2, transcript_id))
  transcript2_introns$intron_number = 1:NROW(transcript2_introns)
  transcript2_introns$type = "intron"
  
  transcript1_introns = ordered_introns %Q% (grepl(transcript1, transcript_id))
  transcript1_introns$intron_number = 1:NROW(transcript1_introns)
  transcript1_introns$type = "intron"
  
  # Combine exons and introns
  transcript2_ei = gUtils::grbind(transcript2_exons, transcript2_introns) %>% sort()
  transcript1_ei = gUtils::grbind(transcript1_exons, transcript1_introns) %>% sort()
  
  # Create windows around specified exon junctions - match your exact approach
  transcript2_windows = GenomicRanges::reduce(transcript2_ei %Q% (exon_number %in% transcript2_exon | intron_number %in% transcript2_intron))
  transcript1_windows = GenomicRanges::reduce(transcript1_ei %Q% (exon_number %in% transcript1_exon | intron_number %in% transcript1_intron))
  
  # Combine windows - maintain the exact order as in your code (MN1 first, then CXXC5)
  mc_windows = c(transcript2_windows, transcript1_windows)
  
  # Set genome style
  GenomeInfoDb::seqlevelsStyle(mc_windows) = "UCSC"
  
  # Read BAM file exactly as in your code
  reads = bamUtils::read.bam(bam_path, intervals = mc_windows + 1000, )
  reads_unlisted = gUtils::grl.unlist(reads)
  
  # Fix chromosome naming
  grgenome = khtools::gr.genome()
  reads_unlisted = gUtils::gr.fix(reads_unlisted, grgenome)
  
  # Create junction window - option to flip
  if(flip_transcripts) {
    junction_window = GRangesList(
      c(gUtils::gr.flipstrand(mc_windows[2]), mc_windows[1])
    )
  } else {
    junction_window = GRangesList(
      c(gUtils::gr.flipstrand(mc_windows[1]), mc_windows[2])
    )
  }
  
  # Create junction object
  jj_window = gGnome::jJ(junction_window)
  
  # Modify junction support function for custom analysis
  my_junction_support <- skitools::junction.support
  
  body_lines <- deparse(body(my_junction_support))
  modified_lines <- gsub("gGnome::merge", "gGnome:::merge", body_lines)
  new_body <- parse(text = modified_lines)
  body(my_junction_support) <- new_body
  
  # Run junction support analysis
  jun = my_junction_support(reads_unlisted, junctions = jj_window, realign = FALSE)
  
  # Track visualization setup
  junc_reads_grl = GenomicRanges::reduce(jun %>% split(jun$grl.ix))
  names(junc_reads_grl) = NULL
  
  # Create gTracks
  gt_junc = gTrack::gTrack(junc_reads_grl, draw.paths = TRUE)
  gt_mc = gTrack::gTrack(mc_windows, draw.paths = TRUE)
  
  # Get gene annotations exactly as in your code
  td.ge.collapsed = gTrack::track.gencode(
    cached = TRUE,
    gene.collapse = TRUE,
    cached.dir="~/DB/GENCODE/hg38/v29/",
    cached.path.collapsed="~/DB/GENCODE/hg38/v29/gencode.composite.collapsed.rds", 
    cached.path="~/DB/GENCODE/hg38/v29/gencode.composite.rds"
  )
  
  td.ge.notcollapsed = gTrack::track.gencode(
    cached = TRUE,
    gene.collapse = FALSE,
    cached.dir="~/DB/GENCODE/hg38/v29/",
    cached.path.collapsed="~/DB/GENCODE/hg38/v29/gencode.composite.collapsed.rds", 
    cached.path="~/DB/GENCODE/hg38/v29/gencode.composite.rds"
  )
  
  # Process data for visualization
  td.dat = td.ge.notcollapsed@data[[1]]
  td.ge.copy = khtools::copy3(td.ge.notcollapsed)
  
  # Filter for transcripts of interest exactly as in your code
  terms <- c(transcript1, transcript2)
  pattern <- paste(terms, collapse = "|")
  td.ge.copy@data[[1]] <- td.dat[grepl(pattern, names(td.dat))]
  
  # Formatting for visualization
  df = gTrack::formatting(gt_junc)
  df$height = 30
  gTrack::formatting(gt_junc) = df
  
  # Set default output filename
  if(is.null(filename) && output_format == "png") {
    base_filename = basename(bam_path)
    filename = paste0(
      "fusion_", 
      gsub("\\.bam$", "", base_filename), 
      "_", 
      paste(c(transcript2, transcript1), collapse="_"), 
      ".png"
    )
  }
  
  # Generate visualization based on requested format
  if(output_format == "svg") {
    show_svg(
      plot(
        c(td.ge.copy, gt_mc, gt_junc),
        window = GenomicRanges::reduce(jun + 100, ignore.strand = TRUE) + 5000
      ),
      height = 20,
      width = 15
    )
    message("SVG visualization created")
  } else if(output_format == "pdf") {
    # If filename is provided, use it with pdf extension
    pdf_file = if(!is.null(filename)) {
      gsub("\\.png$", ".pdf", filename)
    } else {
      base_filename = basename(bam_path)
      paste0(
        "fusion_", 
        gsub("\\.bam$", "", base_filename), 
        "_", 
        paste(c(transcript2, transcript1), collapse="_"), 
        ".pdf"
      )
    }
    
    pdf(pdf_file, height = 20, width = 15)
    plot(
      c(td.ge.copy, gt_mc, gt_junc),
      window = GenomicRanges::reduce(jun + 100, ignore.strand = TRUE) + 5000
    )
    dev.off()
    message(paste("PDF saved to:", pdf_file))
  } else if(output_format == "png") {
    # Use ppng as in your example when filename is provided
    ppng(
      plot(
        c(td.ge.copy, gt_mc, gt_junc),
        window = GenomicRanges::reduce(jun + 100, ignore.strand = TRUE) + 5000
      ),
      height = 20,
      width = 15,
      filename = filename
    )
    message(paste("PNG saved to:", filename))
  }
  
  # Return useful objects
  return(list(
    junction = jun,
    reads = reads_unlisted,
    mc_windows = mc_windows,
    transcripts = c(transcript1, transcript2),
    exon_numbers = c(transcript1_exon, transcript2_exon),
    intron_numbers = c(transcript1_intron, transcript2_intron)
  ))
}
