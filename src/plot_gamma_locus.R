library("rtracklayer")
library('GenomicRanges')
library(dplyr)
library(purrr)
library(zoo)

f_genes <- "data/hg38_GencodeBasic_V43.gz"

source("src/plot_coverage_functions.R")

# ATAC = TRUE
# sup = FALSE
# bcl11a = FALSE
# alpha = FALSE
# myc = FALSE
# scale_up = 1.15 # scale factor for coverage pot (to scale op the max y-axis)

# ---------------------------------------------
# ATAC

transcripts_of_interest <- c()
if (ATAC) {
  p_bigwigs <- "/Users/richardgremmen/Documents/projects/HbF_react/bigwigs atac/output_HbF_react_no_dedup"
  seq = "ATAC"
  if (bcl11a) {
    f_fig_specs <- "data/ATAC_coverage_figure_bcl11a.tsv"
    post_fix <- "_ATAC_coverage.pdf"
    # https://www.ncbi.nlm.nih.gov/gene/?term=AJ404611
    # XL = "ENST00000335712.11"
    # https://www.ncbi.nlm.nih.gov/gene/?term=AJ404612
    # L = "ENST00000356842.9"
    
    # ENST00000642384.2 (top of UCSC browser)
    
    # chosen: ENST00000356842.9
    transcripts_of_interest <- c("ENST00000356842.9")
  }
  else if (sup) {
    f_fig_specs <- "data/ATAC_coverage_figure_sup.tsv"
    post_fix <- "_ATAC_coverage_sup.pdf"
  } else {
    f_fig_specs <- "data/ATAC_coverage_figure.tsv"
    post_fix <- "_ATAC_coverage.pdf"
  }
} else { # CR
  p_bigwigs <- "/Users/richardgremmen/Documents/projects/VER9574/CR/low_q_mapping/"
  seq = "CutnRun"
  if (alpha) {
    f_fig_specs <- "data/CR_coverage_figure.tsv"
    post_fix <- "_CR_coverage.pdf"    
  } else if (sup) {
    f_fig_specs <- "data/CR_coverage_figure_sup.tsv"
    post_fix <- "_CR_coverage_sup.pdf"
  } else {
    f_fig_specs <- "data/CR_coverage_figure.tsv"
    post_fix <- "_CR_coverage.pdf"    
  }
}

genes <- read.delim(f_genes, header = TRUE, sep = "\t", stringsAsFactors = F)

paths_bigwig = list(
  list(bw_folder = p_bigwigs)
)

df_all_bws <- get_df_of_all_bws(paths_bigwig = paths_bigwig,
                                prefix = "")

admin_df <- read.table("data/HbF_sample_admin.tsv",
                       header=TRUE, comment.char = "#")
nrow(admin_df %>% filter(for_publication == "Yes"))
# 15

pub_df <- admin_df %>% 
  filter(library == seq)

nrow(pub_df)

df_all_bws_annot <- df_all_bws %>% 
  left_join(admin_df, by = c("name" = "expname"))

if (ATAC) {

  # cache TSS signal for all ATAC samples (independent of locus)
  # that will be published
  f_atac_tss_signal <- "data/temp/ATAC_tss_signal.tsv"
  
  if (!file.exists(f_atac_tss_signal)) {
  
    non_overlapping_tss_regions <- get_non_overlapping_tss_regions(genes, width = 1000)
    
    df_all_bws_annot <- df_all_bws_annot %>%
      rowwise() %>% 
      mutate(
        tss_coverage_1000 = calc_tss_coverage(
          bw = bw,
          tss_regions = non_overlapping_tss_regions )    
      ) %>% ungroup()
    
    df_all_bws_annot %>% write.table(
      file = f_atac_tss_signal,
      sep = '\t',
      row.names = F,
      quote = F
    )

  } else {
    df_all_bws_annot <- read.table(f_atac_tss_signal, header = T, sep = "\t")
  }  
}

# read figure specifications
df_fig_specs <- read.delim(f_fig_specs, header = TRUE, sep = "\t", 
                           stringsAsFactors = F, comment.char = "#")

df_fig_specs_2 <- df_fig_specs %>% 
  left_join(df_all_bws_annot, by=c("expname" = "name"))


if (ATAC) {
  df_sorted_bws <- df_fig_specs_2 %>% 
      mutate(name = ifelse(is.na(name_publication), expname, name_publication)) %>% 
      dplyr::select(name, bw, tss_coverage_1000)
} else {
  df_sorted_bws <- df_fig_specs_2 %>% 
    mutate(name = ifelse(is.na(name_publication), expname, name_publication)) %>% 
    dplyr::select(name, bw)
}

# ---------------------------------------------------------------
locus_gamma <- list(locus_name = "gammaglobin", chr = "chr11", 
                    # start = 5210200, end = 5313900, 
                    start = 5200000, end = 5320000, 
                    genes_of_interest = c("HBB", "HBD", "HBG2", "HBG1", "HBE1"))
locus_alpha <- list(locus_name = "alphaglobin", chr = "chr16", start = 150000, end = 190000,
                    genes_of_interest = c("HBZ", "HBM", "HBA1", "HBA2"))

#chr2:60,363,406-60,645,549
locus_bcl11a <- list(locus_name = "bcl11a", chr = "chr2", start = 60350000, end = 60645549,
                     genes_of_interest = c())

# scale on either HBA1/HBA2 or HS3 of the LCR
# chr16:172,452-179,127
locus_hba1_hba2 <- list(locus_name = "hba1_hba2", chr = "chr16", start = 172452, end = 179127,
                         genes_of_interest = c("HBA1", "HBA2"))

locus_myc <- list(locus_name = "myc", chr = "chr8", start =127730413, end=127747382,
                  genes_of_interest = c("MYC"))

zoom_gamma <- get_zoom(locus_gamma)
zoom_alpha <- get_zoom(locus_alpha)
zoom_hba1_hba2 <- get_zoom(locus_hba1_hba2)
zoom_bcl11a <- get_zoom(locus_bcl11a)

zoom_myc <- get_zoom(locus_myc)

if (ATAC) {
  df_sorted_bws <- df_sorted_bws %>% dplyr::select(1:3)
} else {
  df_sorted_bws <- df_sorted_bws %>% dplyr::select(1:2)
}


df_sorted_bws$bw_base <- basename(df_sorted_bws$bw)

if (bcl11a) {
  zoom <- zoom_bcl11a
  locus <- locus_bcl11a
} else if (alpha ) {
  zoom <- zoom_alpha
  locus <- locus_alpha
} else if (myc) {
  zoom <- zoom_myc
  locus <- locus_myc
} else {
  zoom <- zoom_gamma
  locus <- locus_gamma
}

# 
df_sorted_bws <- df_sorted_bws %>%
  rowwise() %>% 
  mutate(
    max_signal = get_max_locus(bw = bw, zoom = zoom)
  ) %>% ungroup()

# scale on either TSS coverage or HS3 of the LCR
if (ATAC) {

  df_sorted_bws <- df_sorted_bws %>%
    rowwise() %>% 
    mutate(
      max_gamma = get_max_locus(bw = bw, zoom = zoom_gamma)
    ) %>% ungroup()
  
  # scaling for ATAC based on the coverage around all TSS sites (-500bp - +500bp)
  # previously we used the hs3 of the LCR
  non_overlapping_tss_regions <- get_non_overlapping_tss_regions(genes, width = 1000)
  
  sum_length <- sum(width(non_overlapping_tss_regions))
  
  df_sorted_bws$tss_average_1000 <- df_sorted_bws$tss_coverage_1000 / sum_length
  
  max_scaled <- round( max(df_sorted_bws$max_signal / df_sorted_bws$tss_average_1000) / 10) * 10
  scale_tick_round <- 1
  
} else {

  df_sorted_bws <- df_sorted_bws %>%
    rowwise() %>% 
    mutate(
      sum_hba = get_sum_locus(bw = bw, zoom = zoom_hba1_hba2)
    ) %>% ungroup()
  
  # scaling for CR based on the coverage at HBA1/HBA2
  df_sorted_bws$hba_average <- df_sorted_bws$sum_hba / width(zoom_hba1_hba2)

  max_scaled <- round(max(df_sorted_bws$max_signal / df_sorted_bws$hba_average))
  scale_tick_round <- 1
}


y_max <- max_scaled
hline <- max_scaled / 2
# hline <- 0 (if you want to remove the red line)

locus_name <- locus$locus_name
genes_of_interest <- locus$genes_of_interest


p_output <- paste0("output/","figures/")

f_name <- paste0(p_output, locus_name, post_fix)
dir.create(p_output, showWarnings = FALSE)


# --------------------------------
# start plot
# --------------------------------
pdf(f_name, height=11.6, width=10)

mat.row <- nrow(df_sorted_bws) + 1

layout.matrix <- matrix(1:mat.row, nrow = mat.row, ncol = 1, byrow = TRUE)
layout(mat = layout.matrix
       , height = c(rep(lcm(2), mat.row - 1), lcm(4))
       , width = c(rep(lcm(20), mat.row))
)

width_label <- 12
  
par(oma = c(0, 0, 0, 0))
par(mar=c(0, width_label, 0.5, 2) +0.1)


gRNA_114<-GRanges(seqnames='chr11', ranges=IRanges(start=5280665, end=5280684))
gRNA_503<-GRanges(seqnames='chr11', ranges=IRanges(start=5255265, end=5255284))

# For ATAC bcl11a we want the max scaled over all samples 
if (ATAC & bcl11a) {
  # max scaled over all samples
  max_max_scaled <- 0
  
  for (i in 1:nrow(df_sorted_bws)) {
    
    scale_factor <- df_sorted_bws$tss_average_1000[i]
    
    # for bcl11a: 
    max_signal <- df_sorted_bws$max_signal[i]
    
    max_scaled <- max_signal / scale_factor
    max_max_scaled <- max(max_scaled, max_max_scaled)    
  }
}

for (i in 1:nrow(df_sorted_bws)) {
  
  name <- df_sorted_bws$name[i]
  name <- gsub("_", " ", name)

  bw = df_sorted_bws$bw[i]

  cov<-import(bw, which=zoom)

  if (ATAC) {
    scale_factor <- df_sorted_bws$tss_average_1000[i]
    
    # we scale on max_gamma, also for other loci (expect for bcl11a, see below)
    max_signal <- df_sorted_bws$max_gamma[i]
    
    max_scaled <- max_signal / scale_factor

    y_max <- max_scaled
    hline <- 0

    if (bcl11a) {

      # Overrule the above: use max_max_scaled over all experiments in figure
      y_max <- max_max_scaled
      hline <- 0
    }

  } else {
    scale_factor <- df_sorted_bws$hba_average[i]
  }

  # for BCL11a use binned_bars / running_mean with 50bp window
  # for all other loci use "raw"
  plot_bw(cov = cov, zoom = zoom, 
          method = "raw", # raw, binned , binned_bars
          field = "score", # score, running_mean
          bin_size = 50,
          yMax=y_max, name = name, 
          scale_up = scale_up,
          scale_factor = scale_factor,
          scale_tick_round = scale_tick_round,
          hline = hline)

  show_single_guide <- (locus_name == "gammaglobin" && 
                         grepl("115", name))

  if (show_single_guide) {

    #115_HBG2	chr11	-	5249959	5249978	 CCTTGACCAATAGCCTTGAC
    #115_HBG1	chr11	-	5254883	5254902	 CCTTGACCAATAGCCTTGAC

    x_coords <- c(5249959/1e6, 5254883/1e6)

    line_height <- 0.3  # Height of the lines as a fraction of the maximum y-axis value

    for (x in x_coords) {
      y_line <- max(cov$score) * line_height
      segments(x0 = x, y0 = 0, x1 = x, y1 = y_line, lwd = 2, col = 'blue')
    }

  }
  
  show_bcl11a_guide <- (locus_name == "bcl11a" && 
                          grepl("Bcl11a", name))
  
  if (show_bcl11a_guide) {
    
    #chr2	-	60495264	60495283	 GTGATAAAAGCAACTGTTAG
    x_coords <- c(60495264/1e6, 60495283/1e6)
    
    line_height <- 0.3  # Height of the lines as a fraction of the maximum y-axis value
    
    for (x in x_coords) {
      y_line <- max(cov$score) * line_height
      segments(x0 = x, y0 = 0, x1 = x, y1 = y_line, lwd = 2, col = 'blue')
    }
    
  }

  show_del <- F  
  if (show_del) {
    # Calculate the height of the rectangle (20% of max_score)
    rect_height <- 0.2 * max(cov$score)
    
    # Find the coordinates for the rectangle
    rect_x_start <- start(gRNA_503)/1e6
    rect_x_end <- start(gRNA_114)/1e6
    rect_y <- 0  # Start from the bottom of the plot
    
    # Draw the rectangle with diagonal stripes
    rect(rect_x_start, rect_y - (rect_height/2), rect_x_end, rect_y + (rect_height/2), 
         col = "blue", density = 15,
         border = "transparent")
  }
    
}

par(mar=c(6, width_label, 0, 2) +0.1) 

longestOnly = (locus_name %in% c("gammaglobin") )

plotGlobinGenes (genes=genes, plotZoom=zoom, 
                 genes_of_interest = genes_of_interest,
                 transcripts_of_interest = transcripts_of_interest, 
                 xaxisUnit ="Mb", genename=TRUE, 
                 colour='black', longestOnly=longestOnly,
                 xlabel = "")

dev.off()
