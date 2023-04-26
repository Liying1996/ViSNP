#' Return barplots of associated GWAS phenotypes of SNPs.
#'
#' @param snps_loc Required. snps_loc Required. The locations of SNPs, the format should be "chr1:1014863".
#' @param assembly Optional. The assembly version of the input SNPs. "hg19" and "hg38" can be selected. Default is "hg38". Default is "hg38".
#' @param show_num Optional. The number of associated phenotypes shown on the plot. Default is 5.
#' @param enrichment Optional. Whether to do the enrichment analysis of GWAS studies. Default is FALSE. (Users can do the enrichment analysis by analyze_gwas_enrich() as well)
#'
#' @return A barplot.
#' @export
#'
#' @examples
#' plot_batch_gwas(data)

plot_batch_gwas <- function(snps_loc, assembly="hg38",show_num=5, enrichment=FALSE){
  gwas_data <- gwas_data
  colors <- pal_igv(alpha = 0.8)(show_num)

  if (assembly=="hg38"){
    input_snps <- unique(snps_loc)
  }else{
    hg19_loc <- unique(snps_loc)
    input_snps <- c()
    for (s in hg19_loc){
      input_snps <- c(input_snps, hg19tohg38(s))
    }
  }
  chrom <- str_split_fixed(input_snps, ":", 2)[,1]
  pos <- as.numeric(str_split_fixed(input_snps, ":", 2)[,2])

  gwas_data$location <- paste("chr", gwas_data$CHR_ID, ":", gwas_data$CHR_POS, sep="")

  gwas_df <- gwas_data[gwas_data$location %in% input_snps, c('location', 'DISEASE.TRAIT')]
  gwas_df <- unique(gwas_df)

  gwas_counts <- data.frame(table(gwas_df$DISEASE.TRAIT))
  colnames(gwas_counts) <- c("Phenotype", "Freq")
  gwas_counts <- gwas_counts[order(gwas_counts$Freq, decreasing = T), ]
  gwas_counts <- gwas_counts[1:show_num,]
  gwas_counts <- gwas_counts[order(gwas_counts$Freq, decreasing = F), ]
  gwas_counts$Phenotype <- factor(gwas_counts$Phenotype, levels = gwas_counts$Phenotype)

  g <- ggplot(gwas_counts, aes(x = Phenotype, y = Freq, fill = Phenotype)) + geom_bar(stat = "identity") +
      coord_flip() +
      guides(fill = "none") +
      labs(y = "Frequency") +
      labs(title = "Associated GWAS phenotypes") +
      theme_bw() +
      theme(axis.text = element_text(size=10), axis.text.x = element_text(vjust=0.5)) +
      theme(plot.title = element_text(hjust = 0.5, size=16)) +
      scale_fill_manual(values=colors) +
      scale_x_discrete(labels=function(x) str_wrap(x, width=20))

  print(g)
  return(g)

  if (enrichment){
    g2 <- analyze_gwas_enrich(snps_loc=snps_loc, assembly=assembly)
    print(g2)
  }

}
