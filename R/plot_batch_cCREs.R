#' Return barplots of cCREs that SNPs overlapped.
#'
#' @param data Required. The annotation results from VEP.
#' @param assembly Optional. "hg19" and hg38" can be selected. Default is "hg38".
#' @param show_unclassified Optional. Whether to display SNPs without overlapping cCREs.
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' plot_batch_cCREs(data, assembly = "hg38", show_unclassified=FALSE)


plot_batch_cCREs <- function(data, assembly = "hg38", show_unclassified=FALSE){
  data <- data[data$VARIANT_CLASS=="SNV",]
  if (assembly=="hg19"){
    cCRE_data <- cCRE_data_hg19
  }else{
    cCRE_data <- cCRE_data_hg38
  }
  colnames(cCRE_data) <- c("chrom", "start", "end", "accession", "SCREEN_accession", "type")

  location <- unique(data$Location)

  cCREs_overlap <- c()
  for (loc in location){
    chrom <- str_split_fixed(loc, ":", 2)[1]
    pos <- as.numeric(str_split_fixed(loc, ":", 2)[2])
    pos_bed <- pos - 1

    intersect_cCRE_type <- cCRE_data[(cCRE_data$chrom==chrom)&(cCRE_data$start<=pos_bed)&(cCRE_data$end>=pos_bed), 6]
    if (length(intersect_cCRE_type)==0){
      cCREs_overlap <- c(cCREs_overlap, "Unclassified")
    }else{
      cCREs_overlap <- c(cCREs_overlap, intersect_cCRE_type)
    }
  }

  cCREs_overlap <- gsub(",CTCF-bound", "", cCREs_overlap)
  cCRE_counts <- data.frame(table(cCREs_overlap))
  colnames(cCRE_counts) <- c("cCRE", "Freq")

  if (show_unclassified){
    cCRE_counts$cCRE <- factor(cCRE_counts$cCRE, levels = c("dELS", "pELS","PLS", "DNase-H3K4me3", "CTCF-only", "Unclassified"))
    colors <- c("#FFCD00", "#FFA700", "#FF0000", "#ffaaaa", "#00B0F0", "gray")
  }else{
    cCRE_counts <- cCRE_counts[cCRE_counts$cCRE!="Unclassified",]
    cCRE_counts$cCRE <- factor(cCRE_counts$cCRE, levels = c("dELS", "pELS","PLS", "DNase-H3K4me3", "CTCF-only"))
    colors <- c("#FFCD00", "#FFA700", "#FF0000", "#ffaaaa", "#00B0F0")
  }
  g <- ggplot(cCRE_counts, aes(x = cCRE, y = Freq, fill = cCRE)) + geom_bar(stat = "identity") +
      guides(fill = "none") +
      labs(y = "Frequency") +
      theme_snp() +
      theme(axis.text.x = element_text(vjust=0.5), axis.text = element_text(size=10)) +
      scale_fill_manual(values=colors)
  print(g)
}


