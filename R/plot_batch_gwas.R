#' Return barplots of associated GWAS phenotypes of SNPs.
#'
#' @param data Required. The annotation results from VEP.
#' @param assembly Optional. The assembly version of the input SNPs. "hg19" and "hg38" can be selected. Default is "hg38". Default is "hg38".
#' @param show_num Optional. The number of associated phenotypes shown on the plot. Default is 5.
#'
#' @return A barplot.
#' @export
#'
#' @examples
#' plot_batch_gwas(data)


plot_batch_gwas <- function(data, assembly="hg38", show_num=5){
  gwas_data <- gwas_data
  colors <- pal_igv(alpha = 0.8)(show_num)

  if (assembly=="hg38"){
    input_snps <- data$Location
  }else{
    hg19_loc <- data$Location
    input_snps <- c()
    for (s in hg19_loc){
      input_snps <- c(input_snps, hg19tohg38(s))
    }
  }
  chrom <- str_split_fixed(input_snps, ":", 2)[,1]
  pos <- as.numeric(str_split_fixed(input_snps, ":", 2)[,2])

  gwas_data$location <- paste("chr", gwas_data$CHR_ID, ":", gwas_data$CHR_POS, sep="")

  gwas_df <- gwas_data[gwas_data$location %in% input_snps, ]

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
      theme_snp() +
      labs(title = "Associated GWAS phenotypes") +
      theme(axis.text.x = element_text(vjust=0.5)) +
      scale_fill_manual(values=colors) +
      scale_x_discrete(labels=function(x) str_wrap(x, width=10))

  print(g)

  return(gwas_df)
}
