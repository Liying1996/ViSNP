#' Return the enrichment of GWAS studies.
#'
#' @param snps_loc Required. The locations of SNPs, the format should be "chr1:1014863".
#' @param assembly Optional. "hg19" and hg38" can be selected. Default is "hg38".
#' @param
#'
#' @return A plot.
#' @export
#'
#' @examples
#' analyze_gwas_enrich(snps_loc)


analyze_gwas_enrich <- function(snps_loc, assembly='hg38'){
  gwas_data <- gwas_data

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
  gwas_df <- gwas_data[gwas_data$location %in% input_snps, ]
  gwas_df <- unique(gwas_df[c('location', 'DISEASE.TRAIT')])


  # eg. blood protein levels
  total_count <- length(pos)
  # input_enrich_count <- length(unique(gwas_df[gwas_df$DISEASE.TRAIT=="Blood protein levels","location"]))
  input_enrich_count <- length(unique(gwas_df$location))
  # input_percentage <- input_enrich_count/total_count * 100

  control_samples <- control_samples # Load rds control
  colnames(control_samples) <- c('chrom', 'pos', 'rsID')
  gwas_data_control <- gwas_data[paste('rs', gwas_data$SNP_ID_CURRENT, sep='') %in% control_samples$rsID, ]
  gwas_data_control$SNP_ID_CURRENT <- paste('rs', gwas_data_control$SNP_ID_CURRENT, sep='')

  # control_percents <- c()
  control_num <- c()
  for (i in 1:500){
    controls <- control_samples[sample(nrow(control_samples), total_count),]
    control_df <- gwas_data_control[gwas_data_control$SNP_ID_CURRENT %in% controls$rsID, ]
    control_df <- unique(control_df[c('location', 'DISEASE.TRAIT')])
    # control_enrich_count <- length(unique(control_df[control_df$DISEASE.TRAIT=="Blood protein levels","location"]))
    control_enrich_count <- length(unique(control_df$location))
    control_num <- c(control_num, control_enrich_count)
  }

  control_mean <- mean(control_num)
  control_sd <- sd(control_num)

  zscore <- (input_enrich_count - control_mean) / control_sd
  pvalue <- min(2 * pnorm(-zscore), 1)
  pvalue <- signif(pvalue, 3)

  plot_df <- data.frame(disease=rep("ALL", 2), type=c('Input', 'Control') , Enrich=c(input_enrich_count, control_mean), sd=c(0, control_sd))

  plot_df$type <- factor(plot_df$type, levels=c('Input', 'Control'))

  g <- ggplot(plot_df,aes(type, Enrich, fill=type)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=Enrich-sd, ymax=Enrich+sd), width=0.4, color=c('white', 'black')) +
    scale_fill_jco() +
    labs(x = '', y = 'Frequency') +
    annotate('text', x = 1, y = max(plot_df$Enrich) * 1.15, label=substitute(paste(italic('p'), ' = ', x), list(x=pvalue)), size=4.5) +
    guides(fill=guide_legend(NULL)) +
    theme_bw() +
    theme(axis.title = element_text(size=15), axis.text = element_text(size=12))
  return(g)
}
