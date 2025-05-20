analyze_fusion_junction <- function(bam_path, 
                                   enst_ids = c("ENST00000302326", "ENST00000302517"),
                                   exon_numbers = c(1, 2),
                                   intron_numbers = c(1, 1),
                                   output_file = NULL,
                                   output_format = "svg") {
  # Load required libraries exactly as in the original code
  library(gGnome)
  library(gUtils)
  library(khtools)
  library(BiocManager)
  library(rtracklayer)
  
  # Import GENCODE annotations exactly as in the original code
  gr_v29 = rtracklayer::import.gff3("~/DB/GENCODE/hg38/v29/gencode.v29.annotation.gff3")
  GenomeInfoDb::seqlevelsStyle(gr_v29) = "NCBI"
  gencode = readRDS("~/DB/GENCODE/hg38/v29/gencode.v29.annotation.nochr.rds")
  
  # Load RefSeq lookup data
  refseq_lookup = readRDS("~/DB/GENCODE/refseq_ensembl_lookup.rds")
  data.table::setDT(refseq_lookup)
  
  # Get transcript IDs from the parameters
  refseq_lookup[refseq_mrna %in% c("NM_002430", "NM_016463")]$ensembl_transcript_id %>% base::dput()
  c(enst_ids)
  
  # Extract exons for each transcript
  cxxc5_exons = gencode %Q% which(grepl(enst_ids[1], transcript_id)) %Q% (type == "exon")
  mn1_exons = gencode %Q% which(grepl(enst_ids[2], transcript_id)) %Q% (type == "exon")
  
  # Re-load khtools to ensure operator availability
  library(khtools)
  
  # Process all exons
  exons = gencode %Q% (type == "exon")
  # exons by tx = exons %>% gGnome::gr_construct_by("transcript_id")
  exons_by_tx <- split(exons, exons$transcript_id)
  
  # Get exon and intron ranges
  exon_ranges = range(exons_by_tx)
  intron_ranges = GenomicRanges::setdiff(exon_ranges, exons_by_tx)
  
  # Process intron ranges
  indices <- rep(seq_along(intron_ranges), lengths(intron_ranges))
  tx_ids <- names(intron_ranges)[indices]
  
  intron_ranges_clean <- unlist(intron_ranges)
  mcols(intron_ranges_clean)$transcript_id <- tx_ids
  
  # Order introns exactly as in the original code
  ordered_introns <- intron_ranges_clean[order(ifelse(strand(intron_ranges_clean) == "+", 1, -1) * start(intron_ranges_clean))]
  
  ordered_introns <- (
    gUtils::gr2dt(ordered_introns) %>%
      .[, intron_number := seq_len(.N), by = transcript_id] %>%
      gUtils::dt2gr()
  )
  
  # Process each transcript's introns
  mn1_introns = ordered_introns %Q% (grepl(enst_ids[2], transcript_id))
  mn1_introns$intron_number = 1:NROW(mn1_introns)
  mn1_introns$type = "intron"
  
  cxxc5_introns = ordered_introns %Q% (grepl(enst_ids[1], transcript_id))
  cxxc5_introns$intron_number = 1:NROW(cxxc5_introns)
  cxxc5_introns$type = "intron"
  
  # Combine exons and introns
  mn1_ei = gUtils::grbind(mn1_exons, mn1_introns) %>% sort()
  cxxc5_ei = gUtils::grbind(cxxc5_exons, cxxc5_introns) %>% sort()
  
  # Create windows around specified exon junctions
  mn1_windows = GenomicRanges::reduce(mn1_ei %Q% (exon_number %in% exon_numbers[2] | intron_number %in% intron_numbers[2]))
  cxxc5_windows = GenomicRanges::reduce(cxxc5_ei %Q% (exon_number %in% exon_numbers[1] | intron_number %in% intron_numbers[1]))
  
  # For debugging, display what we found
  print(mn1_windows)
  print(cxxc5_windows)
  
  # Combine windows
  mc_windows = c(mn1_windows, cxxc5_windows)
  
  # Set genome style
  GenomeInfoDb::seqlevelsStyle(mc_windows) = "UCSC"
  
  # Read BAM file
  reads = bamUtils::read.bam(bam_path, intervals = mc_windows + 1000, )
  reads_unlisted = gUtils::grl.unlist(reads)
  
  # Fix chromosome naming
  grgenome = khtools::gr.genome()
  reads_unlisted = gUtils::gr.fix(reads_unlisted, grgenome)
  
  # Create junction window
  junction_window = GRangesList(
    c(gUtils::gr.flipstrand(mc_windows[1]), mc_windows[2])
  )
  junction_window = gUtils::gr.fix(junction_window, grgenome)
  
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
  
  # Get gene annotations
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
  
  # Filter for transcripts of interest
  terms <- c(enst_ids)
  pattern <- paste(terms, collapse = "|")
  td.ge.copy@data[[1]] <- td.dat[grepl(pattern, names(td.dat))]
  
  # Formatting for visualization
  df = gTrack::formatting(gt_junc)
  df$height = 30
  gTrack::formatting(gt_junc) = df
  
  # Set default output filename if not provided
  if(is.null(output_file)) {
    base_filename = basename(bam_path)
    output_file = paste0(
      "fusion_", 
      gsub("\\.bam$", "", base_filename), 
      "_", 
      paste(enst_ids, collapse="_"), 
      ".", 
      output_format
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
    pdf(output_file, height = 20, width = 15)
    plot(
      c(td.ge.copy, gt_mc, gt_junc),
      window = GenomicRanges::reduce(jun + 100, ignore.strand = TRUE) + 5000
    )
    dev.off()
    message(paste("PDF saved to:", output_file))
  } else if(output_format == "png") {
    png(
      output_file,
      height = 20 * 72,
      width = 15 * 72,
      res = 72
    )
    plot(
      c(td.ge.copy, gt_mc, gt_junc),
      window = GenomicRanges::reduce(jun + 100, ignore.strand = TRUE) + 5000
    )
    dev.off()
    message(paste("PNG saved to:", output_file))
  }
  
  # Return useful objects
  return(list(
    junction = jun,
    reads = reads_unlisted,
    mc_windows = mc_windows,
    transcripts = enst_ids,
    exon_numbers = exon_numbers,
    intron_numbers = intron_numbers
  ))
}
