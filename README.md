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

7. **Classifier**  
   ```
   /gpfs/data/imielinskilab/Git/flows/tasks/Classifier_EF.task
   ```

8. **SNV**  
   ```
   /gpfs/home/freite01/lab/Git/flows/tasks/SNV_taps.task
   ```
   followed by
   ```
   /gpfs/data/imielinskilab/Git/flows/tasks/SnpEff_taps.task
   ```
   followed by
   parsesnpeff & other filtering
   
## Other downstream analyses

These, and a few others, are in `/gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/EF_downstream_analyses.ipynb`

9. **CNVs**  
   ```R
   source("~/Projects/TAPS/functions/plots.R") -> plot_gene_cnv
   ```
10. **Fusions**

11. **Complex events**  
    ```
    /gpfs/home/freite01/lab/Git/flows/tasks/JaBbA_taps.task
    ```
    followed by
    ```
    /gpfs/home/freite01/lab/Git/flows/tasks/Events_taps.task
    ```

## Cohort-wide summaries

12. **PCA / UMAP**  
    ```
    /gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/pca.ipynb
    ```

13. **Methylation unsupervised clustering**  
    ```
    /gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/methylation_heatmap.ipynb
    ```

14. **Oncoprint**  
    ```
    /gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/oncoprint.ipynb
    ```
