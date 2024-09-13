library("rtracklayer")
library('GenomicRanges')
library(dplyr)
library(purrr)

f_genes <- "data/hg38_GencodeBasic_V43.gz"

source("src/plot_coverage_functions.R")


# --------------------------------------
# ATAC
# --------------------------------------
scale_up = 1.15
alpha = FALSE
ATAC = TRUE;sup = FALSE; bcl11a = FALSE; myc = FALSE
source("src/plot_gamma_locus.R")

# ATAC MYC (random gene should not change)
scale_up = 1.15
myc = TRUE
alpha = FALSE
ATAC = TRUE;sup = FALSE; bcl11a = FALSE
source("src/plot_gamma_locus.R")

# ATAC sup
myc = FALSE
ATAC = TRUE;sup = TRUE; bcl11a = FALSE
source("src/plot_gamma_locus.R")

# ATAC alpha
alpha = TRUE
ATAC=TRUE; sup = FALSE; bcl11a = FALSE
source("src/plot_gamma_locus.R")

# ATAC BCL11a
scale_up = 1.10
alpha=FALSE;ATAC = TRUE;sup=FALSE;bcl11a = TRUE
source("src/plot_gamma_locus.R")

scale_up = 1.15

# --------------------------------------
# CR
# ---------------------------------------
ATAC = FALSE;sup = FALSE; bcl11a = FALSE; myc=FALSE
source("src/plot_gamma_locus.R")

#CR sup
ATAC = FALSE;sup = TRUE; bcl11a = FALSE
source("src/plot_gamma_locus.R")

# CR alpha
alpha = TRUE
ATAC=FALSE; sup = FALSE; bcl11a = FALSE
source("src/plot_gamma_locus.R")


