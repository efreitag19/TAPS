# eftools.R - Loads all functions from TAPS repository

# Source individual function files
source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/R/snv.R")
source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/R/qual.R")
source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/R/plots.R")
source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/R/functions.R")

load_eftools <- function() {
  source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/R/snv.R")
  source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/R/qual.R")
  source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/R/plots.R")
  source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/R/functions.R")
  
  cat("All eftools functions loaded successfully!\n")
}
