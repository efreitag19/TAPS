# TAPS

TAPS workflow xx Eloise

## processing using next flow
follow instructions in "/gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/TAPS_new.ipynb" after receiving link to fastq's, takes you through the following steps
1. Nextflow
    - parabricks: bams
    - coverage: fragcounter, dryclean
    - sv_calling: gridss
2. Deduping bams (incl. merging if necessary)
3. Picard : "~/tasks/tumor_only/hg38/PicardBamQC.task"
4. CBS : "~/tasks/CBS_ZC.task"
5. Rastair : "~/tasks/Rastair_GRCh38.task"

## analyses to do per each new batch
6. ichorCNA : "/gpfs/data/imielinskilab/git/mskilab/flows/tasks/ichorCNA.task"
7. Classifier : "/gpfs/data/imielinskilab/Git/flows/tasks/Classifier_EF.task"
8. SNV : "/gpfs/home/freite01/lab/Git/flows/tasks/SNV_taps.task"
9. cool amplifications : source("~/Projects/TAPS/functions/plots.R") -> plot_gene_cnv
10. complex events : "/gpfs/home/freite01/lab/Git/flows/tasks/JaBbA_taps.task" -> "/gpfs/home/freite01/lab/Git/flows/tasks/Events_taps.task"

## optional summarizing things
11. PCA / umap : "gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/pca.ipynb"
12. methylation unsupervised clustering : "gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/methylation_heatmap.ipynb"
13. oncoprint : "gpfs/data/imielinskilab/projects/TAPS/wmg-nyu-matija/oncoprint.ipynb"
