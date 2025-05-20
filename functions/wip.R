analyze_fusion_junction <- function(bam_path, 
                                    transcript1 = "ENST00000302326",  # CXXC5
                                    transcript2 = "ENST00000302517",  # MN1
                                    transcript1_exon = 2,
                                    transcript1_intron = 1,
                                    transcript2_exon = 1, 
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
