# eftools.R - Loads all functions from TAPS repository

# Source individual function files
source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/functions/snv.R")
source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/functions/qual.R")
source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/functions/plots.R")
source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/functions/functions.R")

load_eftools <- function() {
  source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/functions/snv.R")
  source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/functions/qual.R")
  source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/functions/plots.R")
  source("/gpfs/data/imielinskilab/home/freite01/git/TAPS/functions/functions.R")
  
  cat("All eftools functions loaded successfully!\n")
}
