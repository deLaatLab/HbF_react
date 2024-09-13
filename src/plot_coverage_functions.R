
# exact plotting of the bigwigs
plot_bars <- function(cov, zoom, yMax, scale_factor,
                      field = c("score", "running_mean"), bin_size) {
  # Set plot parameters
  ylim <- c(0,yMax)
  xlim <- c(start(zoom) / 1e6, end(zoom) / 1e6)
  
  if (field == "running_mean") {

    binned_cov_df <- as.data.frame(cov)
    
    # plot as barplot since we have fixed bins
    binned_cov_df$midpoint <- binned_cov_df$start + (bin_size / 2)

    barplot( height = binned_cov_df$running_mean / scale_factor
            , names.arg = binned_cov_df$midpoint
            , space = 0
            , ylim = c(0, yMax)
            , las = 2, cex.names = 0.5
            , xaxt = "n"
            , axes = FALSE
            , ann = FALSE
            )

  } else {
    
    # old way for unequal bin sizes (takes a long time!)
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, xlab = "", ylab = "",
                cex.axis = 0.9, cex.lab = 1, cex.main = 0.9)

    # Get the number of bins
    num_bins <- length(cov)
    
    # Plot each bin separately
    for (i in 1:num_bins) {
      # Extract bin start and end
      bin_start <- start(cov[i]) / 1e6
      bin_end <- end(cov[i]) / 1e6

      rect(bin_start, 0, bin_end, cov$score[i] / scale_factor, col = "black", border = NA)

    }
  }
}


plot_bw<-function(cov, zoom, yMax=100, name=NULL, scale_up = 1.15, 
                  scale_factor = 1, hline = 40, scale_tick_round = 10,
                  method = c("raw","binned","binned_bars"),
                  field = c("score", "running_mean"),
                  bin_size = 200){
  
  # cov<-import(bw, which=zoom)
  cov$col<-'black'

  # Hack to scale everything up by 15% to display the highest tick (see below)
  # above the maximum score:
  yMax = yMax * scale_up


  if (method == "raw") {
    plot(start(cov)/1e6,  cov$score / scale_factor
         ,ylim=c(0,yMax)
         ,xlim=c(start(zoom)/1e6,end(zoom)/1e6)
         ,xlab="", ylab="",type="h"
         ,frame.plot=FALSE,cex.axis=0.9,cex.lab=1
         ,cex.main=0.9
         ,col=cov$col
         ,axes=FALSE
         ,ann=FALSE
    )
  }
  
  if (method == "binned" | method == "binned_bars") {
  
    # Define bins for smoothing out the signal
    bin_size <- bin_size
    bins <- tile(zoom, n = ceiling(width(zoom) / bin_size))
    
    merged_bins <- unlist(bins)
    
    gr.data.cov <- GenomicRanges::coverage(cov, weight="score")
    chr_data <- gr.data.cov[as.character(seqnames(zoom)[1])]
    
    avg <- binnedAverage(merged_bins, chr_data, "score")    
  
    if (field == "running_mean") {
      
      # Calculate the running mean over 5 bins
      window_size <- 5
      avg$running_mean <- rollapply(avg$score, width = window_size, FUN = mean, align = "center", fill = NA)
      
    }
    avg$col<-'black'

  }

  if (method == "binned") {
    
    if (field == "running_mean") {
      plot(start(avg)/1e6,  avg$running_mean / scale_factor
           ,ylim=c(0,yMax)
           ,xlim=c(start(zoom)/1e6,end(zoom)/1e6)
           ,xlab="", ylab="",type="h"
           ,frame.plot=FALSE,cex.axis=0.9,cex.lab=1
           ,cex.main=0.9
           ,col=avg$col
           ,axes=FALSE
           ,ann=FALSE)
      } else {
        plot(start(avg)/1e6,  avg$score / scale_factor
             ,ylim=c(0,yMax)
             ,xlim=c(start(zoom)/1e6,end(zoom)/1e6)
             ,xlab="", ylab="",type="h"
             ,frame.plot=FALSE,cex.axis=0.9,cex.lab=1
             ,cex.main=0.9
             ,col=avg$col
             ,axes=FALSE
             ,ann=FALSE)
      }
      
    
  }

  if (method == "binned_bars") {
    # optional: exact plotting of the bigwigs (replace plot above)
    plot_bars(cov = avg, zoom = zoom, yMax = yMax, scale_factor = scale_factor,
              field = field, bin_size = bin_size)
  }

  # Calculate even integer values within the range
  even_values <- seq(0, floor(yMax/2)*2, by=scale_tick_round)
  # Determine the highest tick value
  highest_tick <- max(even_values[length(even_values)], 0)
  # Add highest tick label and zero
  axis(2, at=c(0, highest_tick), labels=c(0, highest_tick), las=2)

  mtext(name, side=2, line=3.5, at=yMax/2, adj=1, cex=0.5, las=1)

  if (hline > 0) {
    abline(h = hline, lty = 2, col = "red")
  }

}

something_empty <- function() {
  # xaxt makes it empty
  plot(NULL,
       xlim = c(start(zoom) / 1e6, end(zoom) / 1e6),
       ylim = c(0, 2), ylab = "", xlab = "", yaxt = "none",
       frame.plot = FALSE, cex.axis = 0.9, cex.lab = 1
       , xaxt = "n"
  )
}

plotGlobinGenes <- function(genes, plotZoom, genename = TRUE, 
                            # please provide either a gene OR a transcript
                            # combinations are possible but if they overlap still all transcripts of a gene are chosen
                            genes_of_interest = c("HBD", "HBG2", "HBG1", "HBE1"),
                            transcripts_of_interest = c(),
                            colour = "black", 
                            xaxisUnit = c("Mb", "Kb", "bp"), 
                            longestOnly = TRUE, xlabel = NULL) {
  if (xaxisUnit == "Mb") {
    Div <- 1e6
  } else if (xaxisUnit == "Kb") {
    Div <- 1e3
  } else {
    Div <- 1
  }
  
  chrom <- as.character(seqnames(plotZoom))
  start <- as.numeric(start(plotZoom))
  end <- as.numeric(end(plotZoom))
  
  genes.chr <- genes[genes$chrom == chrom & genes$txStart < end & genes$txEnd > start, ]
  
  # please provide either a gene OR a transcript
  # combinations are possible but if they overlap still all transcripts of a gene are chosen
  genes.chr <- genes.chr[genes.chr$name2 %in% genes_of_interest | 
                           genes.chr$name %in% transcripts_of_interest, ]
  

  if (longestOnly) {
    # keep longest variant
    genes.chr$len <- genes.chr$txEnd - genes.chr$txStart
    genes.chr <- genes.chr[genes.chr$len < 2000, ] # remove long globin genes
    genes.chr <- genes.chr[order(genes.chr$len, decreasing = TRUE), ]
    genes.chr <- genes.chr[!duplicated(genes.chr$name2), ]
    genes.chr <- genes.chr[order(genes.chr$txStart), ]
    
    # longestOnly selects these transcripts
    # View(genes.chr %>% select(name, name2))
    # ENST00000292901.7 HBD
    # ENST00000647020.1 HBB
    # ENST00000396895.3 HBE1
    # ENST00000336906.6 HBG2
    # ENST00000330597.5 HBG1
  }
  
  # plot the exons
  if (is.null(xlabel)) {
    xlabel <- paste0("Position ", chrom, " (", xaxisUnit, ")")
  }

  something_empty()
  
  es <- (as.numeric(unlist(strsplit(as.character(genes.chr$exonStarts), ","))))
  ee <- (as.numeric(unlist(strsplit(as.character(genes.chr$exonEnds), ","))))
  
  strands <- rep(genes.chr$strand, unlist(lapply(strsplit(as.character(genes.chr$exonStarts), ","), length)))
  
  if (genename) {
    yPosexon <- ifelse(strands == "+", 1.05, 0.45)
    yPosGene <- ifelse(genes.chr$strand == "+", 1.3, 0.7)
    yPosGeneName <- ifelse(genes.chr$strand == "+", 1.75, 0.25)
  } else {
    yPosexon <- ifelse(strands == "+", 1.25, 0.25)
    yPosGene <- ifelse(genes.chr$strand == "+", 1.5, 0.5)
    yPosGeneName <- ifelse(genes.chr$strand == "+", 1.65, 0.35)
  }

  rect(xleft = es / Div, ybottom = yPosexon, 
       xright = ee / Div, ytop = yPosexon + 0.5, 
       col = colour, border = colour)
  
  
  # plot the gene
  segments(genes.chr$txStart / Div, yPosGene, genes.chr$txEnd / Div, yPosGene, col = colour)

  if (genename) {
    genesCenter <- (genes.chr$txStart + genes.chr$txEnd) / 2
    if (length(genesCenter) > 0) {
      text(x = genesCenter / Div, y = yPosGeneName, labels = genes.chr$name2, cex = 0.8)
    }
  }
  
  if ("HBG1" %in% genes_of_interest) {
  
    # LCR coordinates: just a stripe just like the gene (without exons)
    lcr_start <-c(5275366)
    lcr_end <- c(5291803)
    
    segments(lcr_start / Div, yPosGene, lcr_end / Div, yPosGene, col = colour)
    
    if (genename) {
      genesCenter <- (lcr_start + lcr_end) / 2
      text(x = genesCenter / Div, y = yPosGeneName, labels = c("LCR"), cex = 0.8)
    }
  }
}


custom_order <- function(x) {
  if (grepl("WT", x)) return(1)
  if (grepl("503", x)) return(2)
  if (grepl("115", x)) return(3)
  if (grepl("Gata1", x)) {
    return(4)
  } else {
    return (9)
  }
  
}

get_df_of_all_bws <- function(paths_bigwig,
    prefix = "VER9574_HbF_react_AF_CutnRun_Hudep-2_") {
  
  df_all_bws <- data.frame()
  
  for (path_bigwig in paths_bigwig) {
    bw <- list.files(path = path_bigwig$bw_folder, pattern = '*.bw', full.names = TRUE)
    sorted_bw <- bw[order(sapply(bw, custom_order))]
    
    df_bws <- data.frame(bw)
    df_bws$mapq <- path_bigwig$mapq
    
    df_bws$name <- df_bws$bw %>% basename() %>% 
      sub(pattern = prefix, replacement = "") %>% 
      sub(pattern = ".bw", replacement = "")
    
    df_all_bws <- rbind.data.frame(df_all_bws, df_bws)
    
    print(path_bigwig$bw_folder)
    
  }
  
  df_all_bws <- df_all_bws %>% arrange(name)
  
  return (df_all_bws)
  
}


get_zoom <- function(locus) {
  zoom <- GRanges(seqnames = locus$chr,
                  ranges = IRanges(start = locus$start, end = locus$end))
}


get_max_locus <- function(bw, zoom) {
  cov <- rtracklayer::import(bw, which=zoom)
  
  max_value <- max(cov$score)
  return(max_value)
} 

get_sum_locus <- function(bw, zoom) {
  cov <- rtracklayer::import(bw, which=zoom)
  
  bin_widths <- width(cov)
  bin_signal <- cov$score
  
  integrated_signal <- sum(bin_signal * bin_widths)
  
  return(integrated_signal)
} 


get_non_overlapping_tss_regions <- function(genes, width = 2000) {
  
  # Convert the data frame to GRanges object
  gene_gr <- GRanges(seqnames = genes$chr,
                     ranges = IRanges(start = genes$txStart, end = genes$txEnd),
                     strand = genes$strand,
                     gene_id = genes$name)
  
  # Extract TSS coordinates based on strand
  pos_strand <- gene_gr[strand(gene_gr) == "+"]
  neg_strand <- gene_gr[strand(gene_gr) == "-"]
  
  tss_gr_pos <- resize(pos_strand, width = 1, fix = "start")
  tss_gr_neg <- resize(neg_strand, width = 1, fix = "end")
  
  tss_gr <- c(tss_gr_pos, tss_gr_neg)
  
  tss_gr_uniq <- tss_gr %>% unique()
  
  # Define regions around TSS (1000bp up- and downstream)
  tss_regions <- resize(tss_gr_uniq, width = width, fix = "center")
  
  # Fuse overlapping or adjacent ranges into non-overlapping ranges
  non_overlapping_tss_regions <- GenomicRanges::reduce(tss_regions)
  
  return(non_overlapping_tss_regions)
}


calc_tss_coverage <- function(bw, tss_regions) {
  
  # Importing as RleList because binnenAverage needs that
  cov <- rtracklayer::import(bw, as = "RleList")
  
  tss_regions <- keepSeqlevels(tss_regions, seqlevels(cov), pruning.mode = "coarse")
  
  avg <- binnedAverage(tss_regions, cov, "score")
  
  bin_widths <- width(avg)
  bin_signal <- avg$score
  
  integrated_signal <- sum(bin_signal * bin_widths)
  
  return (integrated_signal)
  
}
