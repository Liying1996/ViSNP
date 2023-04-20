#' Return barplots of cCREs that SNPs overlapped.
#'
#' @param snps_loc Required.  The locations of SNPs, the format should be like: "chr1:1014863".
#' @param assembly Optional. "hg19" and hg38" can be selected. Default is "hg38".
#' @param show_unclassified Optional. Whether to display SNPs without overlapping cCREs. Default is FALSE.
#' @param enrichment Optional. Whether to do the enrichment analysis of cCREs. Default is FALSE. (Users can do the enrichment analysis by analyze_ccre_enrich() as well)
#' @param show_p Optional.Show p-values or "***" (significant) on enrichment plot. Default is FALSE.
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' plot_batch_cCREs(snps_loc, assembly = "hg38", show_unclassified=FALSE, enrichment=FALSE, show_p=FALSE)


plot_batch_cCREs <- function(snps_loc, assembly = "hg38", show_unclassified=FALSE, enrichment=FALSE, show_p=FALSE){
  # data <- data[data$VARIANT_CLASS=="SNV",]
  if (assembly=="hg19"){
    cCRE_data <- cCRE_data_hg19
  }else{
    cCRE_data <- cCRE_data_hg38
  }
  colnames(cCRE_data) <- c("chrom", "start", "end", "accession", "SCREEN_accession", "type")

  # location <- unique(data$Location)

  # cCREs_overlap <- c()
  # for (loc in location){
  #   chrom <- str_split_fixed(loc, ":", 2)[1]
  #   pos <- as.numeric(str_split_fixed(loc, ":", 2)[2])
  #   pos_bed <- pos - 1
  #
  #   intersect_cCRE_type <- cCRE_data[(cCRE_data$chrom==chrom)&(cCRE_data$start<=pos_bed)&(cCRE_data$end>=pos_bed), 6]
  #   if (length(intersect_cCRE_type)==0){
  #     cCREs_overlap <- c(cCREs_overlap, "Unclassified")
  #   }else{
  #     cCREs_overlap <- c(cCREs_overlap, intersect_cCRE_type)
  #   }
  # }
  test_bed <- data.frame(str_split_fixed(snps_loc, ':', 2))
  colnames(test_bed) <- c('chr', 'start')
  test_bed$start <- as.numeric(test_bed$start) - 1
  test_bed$end <- test_bed$start
  test_total_overlap <- bt.intersect(a = cCRE_data, b = test_bed, wa = TRUE, wb = TRUE)
  test_df <- test_total_overlap[paste(test_total_overlap$V7, test_total_overlap$V8) %in% paste(test_bed$chr, test_bed$start), ]
  test_overlap <- gsub(",CTCF-bound", "", test_df$V6)
  cCRE_counts <- data.frame(table(test_overlap))
  colnames(cCRE_counts) <- c("cCRE", "Freq")

  if (show_unclassified){
    cCRE_counts$cCRE <- factor(cCRE_counts$cCRE, levels = c("dELS", "pELS","PLS", "DNase-H3K4me3", "CTCF-only", "Unclassified"))
    colors <- c("#FFCD00", "#FFA700", "#FF0000", "#ffaaaa", "#00B0F0", "gray")
  }else{
    cCRE_counts <- cCRE_counts[cCRE_counts$cCRE!="Unclassified",]
    cCRE_counts$cCRE <- factor(cCRE_counts$cCRE, levels = c("dELS", "pELS","PLS", "DNase-H3K4me3", "CTCF-only"))
    colors <- c("#FFCD00", "#FFA700", "#FF0000", "#ffaaaa", "#00B0F0")
  }

  g1 <- ggplot(cCRE_counts, aes(x = cCRE, y = Freq, fill = cCRE)) + geom_bar(stat = "identity") +
      guides(fill = "none") +
      labs(y = "Frequency") +
      theme_snp() +
      theme(axis.text.x = element_text(vjust=0.5), axis.text = element_text(size=12)) +
      scale_fill_manual(values=colors)
  print(g1)

  if (enrichment){
    g2 <- analyze_ccre_enrich(snps_loc, assembly = assembly, show_p=show_p)
    print(g2)
  }
}


