# TAPS Workflow

> **TAPS** (TET-assisted pyridine borane sequencing) workflow by Eloise

## Processing using Nextflow

Follow instructions in `/gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/TAPS_new.ipynb` after receiving link to fastq's, which takes you through the following steps:

1. **Nextflow** - Results give you:
   - `/parabricks`: bams
   - `/coverage`: fragcounter, dryclean
   - `/sv_calling`: gridss

Update paths in pairs table

2. **Dedup bams** (including merging if necessary)

3. **Picard** - get wgs_metrics for number of reads, coverage ->
      update smb://research-cifs.nyumc.org/Research/snudem01lab/snudem01labspace/$PROJECTS/CSF cfDNA/TAPS LIB PREPS_RESEARCH/ND Summary Data.xlsx
   ```
   ~/tasks/tumor_only/hg38/PicardBamQC.task
   ```
   
4. **CBS** - get seg.rds and cov.rds for plotting CNVs
   ```
   ~/tasks/CBS_ZC.task
   ```

5. **Rastair** - TAPS-specific methylation caller to get rastair_output.mods which serves as methylation beta matrix
   ```
   ~/tasks/Rastair_GRCh38.task
   ```

## Analyses to do per each new batch

6. **ichorCNA** - generate genomeWide.pdfs from bam, and put plots into smb://research-cifs.nyumc.org/Research/snudem01lab/snudem01labspace/$PROJECTS/CSF cfDNA/TAPS_cnv
   ```
   /gpfs/data/imielinskilab/git/mskilab/flows/tasks/ichorCNA.task
   ```

7. **Classifier** - create bed file from rastair that then inputs into methylation classifier, with output being a .csv with predictions and accompanying plots 
   ```
   /gpfs/data/imielinskilab/Git/flows/tasks/Classifier_EF.task
   ```

8. **SNV** - following along in `/gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/SNV.ipynb` 
   
   ngs_to_bed() to get bed files from PACT data
   
   ```
   /gpfs/home/freite01/lab/Git/flows/tasks/SNV_taps.task
   ```
   followed by
   ```
   /gpfs/data/imielinskilab/Git/flows/tasks/SnpEff_taps.task
   ```
   followed by parsesnpeff & other custom filtering
   
## Other downstream analyses

These, and a few others, are in `/gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/EF_downstream_analyses.ipynb`

9. **CNVs**  
   
   plot_gene_cnv()
   
10. **Fusions**

   analyze_fusion_junction()

11. **Complex events**  
    ```
    /gpfs/home/freite01/lab/Git/flows/tasks/JaBbA_taps.task
    ```
    followed by
    ```
    /gpfs/home/freite01/lab/Git/flows/tasks/Events_taps.task
    ```

## Cohort-wide summaries

12. **PCA / UMAP** - create PCA and UMAP plots from rastair outputs, intersected with epic array defined methylation sites 
    ```
    /gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/pca.ipynb
    ```

13. **Methylation unsupervised clustering**  - create methylation heatmaps with top 5,000 variable sites via unsupervised clustering
    ```
    /gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/methylation_heatmap.ipynb
    ```

14. **Oncoprint**  - plug output from SNV and metadata and (for now) manual visualization of CNVs from genomeWide ichor plots to create oncoprint summarizing current cohort
    ```
    /gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/oncoprint.ipynb
    ```
