#' Analyze the enrichment of cCREs.
#'
#' @param snps_loc Required. The genomic locations of SNPs, the format should be like: "chr1:1014863".
#' @param assembly Optional. "hg19" and hg38" can be selected. Default is "hg38".
#' @param show_p Optional. Show p-values or "***" (significant) on output plot. Default is FALSE.
#'
#' @return A plot.
#' @export
#'
#' @examples
#' analyze_ccre_enrich(snps_loc)


# options(bedtools.path = "~/anaconda3/bin/")
# library('bedtoolsr')

analyze_ccre_enrich <- function(snps_loc, assembly='hg38', show_p=FALSE){
  if (assembly=="hg19"){
    cCRE_data <- cCRE_data_hg19
  }else{
    cCRE_data <- cCRE_data_hg38
  }

  snps_loc <- unique(snps_loc)

  test_bed <- data.frame(str_split_fixed(snps_loc, ':', 2))
  colnames(test_bed) <- c('chr', 'start')
  test_bed$start <- as.numeric(test_bed$start) - 1
  test_bed$end <- test_bed$start
  test_total_overlap <- bt.intersect(a = cCRE_data, b = test_bed, wa = TRUE, wb = TRUE)
  test_df <- test_total_overlap[paste(test_total_overlap$V7, test_total_overlap$V8) %in% paste(test_bed$chr, test_bed$start), ]
  test_overlap <- gsub(",CTCF-bound", "", test_df$V6)
  cCRE_counts <- data.frame(table(test_overlap))
  colnames(cCRE_counts) <- c("cCRE", "Freq")

  total_count <- nrow(test_bed)

  control_samples <- control_samples # Load controls
  colnames(control_samples) <- c('chrom', 'pos', 'rsID')

  control_samples <- control_samples[sort(sample(nrow(control_samples), min(total_count*100), nrow(control_samples))),]

  control_bed <- control_samples[, c('chrom', 'pos')]
  control_bed$pos <- control_bed$pos - 1
  control_bed$end <- control_bed$pos
  colnames(control_bed) <- c('chr', 'start', 'end')

  control_total_overlap <- bt.intersect(a = cCRE_data, b = control_bed, wa = TRUE, wb = TRUE)

  cCRE_type <- c("CTCF-only", "DNase-H3K4me3", "PLS", "dELS", "pELS")
  control_all <- c()
  for (i in 1:500){
    controls <- control_samples[sample(nrow(control_samples), total_count),]
    control_df <- control_total_overlap[paste(control_total_overlap$V7, control_total_overlap$V8) %in% paste(controls$chrom, controls$pos-1), ]
    control_overlap <- gsub(",CTCF-bound", "", control_df$V6)
    control_counts <- data.frame(table(control_overlap))
    colnames(control_counts) <- c("cCRE", "Freq")
    for (c in cCRE_type){
      if (c %in% control_counts$cCRE){
        control_all <- c(control_all, control_counts[control_counts$cCRE==c, 2])
      }else{
        control_all <- c(control_all, 0)
      }
    }
    control_all <- c(control_all, total_count - sum(control_counts$Freq))
  }

  control_all_df <- as.data.frame(matrix(control_all, ncol=6, byrow=TRUE))
  colnames(control_all_df) <- c(cCRE_type, 'Unclassified')
  control_all_df$cCRE <- total_count - control_all_df$Unclassified
  control_all_df <- subset(control_all_df, select = -c(Unclassified))

  mean_control <- apply(control_all_df, 2, mean)
  sd_control <- apply(control_all_df, 2, sd)

  new_cCRE_counts <- c()
  for (i in cCRE_type){
    new_cCRE_counts <- c(new_cCRE_counts, cCRE_counts[cCRE_counts$cCRE==i, 2])
  }
  new_cCRE_counts <- c(new_cCRE_counts, sum(cCRE_counts$Freq))
  ccRE_counts_sd <- rep(0, 6)

  cCRE_df_mean <- data.frame(rbind(new_cCRE_counts, mean_control))
  cCRE_df_sd <- data.frame(rbind(ccRE_counts_sd, sd_control))
  cCRE_df_mean$type <- c('Input', 'Control')
  cCRE_df_sd$type <- c('Input', 'Control')

  p_values <- c()
  for (i in c("cCRE", "dELS", "pELS", "PLS", "DNase.H3K4me3", "CTCF.only")){
    z <- (cCRE_df_mean[1, i] - cCRE_df_mean[2, i])/cCRE_df_sd[2, i]
    p <- 2 * pnorm(-z)
    p_values <- c(p_values, signif(min(1, p), 2))
  }


  melt_mean <- melt(cCRE_df_mean, id.vars = 'type')
  melt_sd <- melt(cCRE_df_sd, id.vars = 'type')

  melt_df <- cbind(melt_mean[, c('type', 'variable', 'value')], melt_sd[, 'value'])
  colnames(melt_df) <- c('type1', 'type2', 'mean', 'sd')
  melt_df$type1 <- factor(melt_df$type1, levels = c('Input', 'Control'))
  melt_df$type2 <- factor(melt_df$type2, rev(levels(melt_df$type2)))


  text_loc_y1 = apply(cCRE_df_mean[,c("cCRE", "dELS", "pELS", "PLS", "DNase.H3K4me3", "CTCF.only")], 2, max)
  text_loc_y = text_loc_y1 + max(text_loc_y1) * 0.1

  segments <- tibble(
    x1 = 1:6-0.2,
    x2 = 1:6+0.2,
    y1 = text_loc_y1 + max(text_loc_y1) * 0.05,
    y2 = text_loc_y1 + max(text_loc_y1) * 0.05
  )


  star_p <- c()
  for (x in p_values){
    if (x <= 1e-7){
      star_p <- c(star_p, '***')
    }else if ((x > 1e-7) & (x <= 1e-3)){
      star_p <- c(star_p, '**')
    }else if ((x > 1e-3) & (x <= 0.05)){
      star_p <- c(star_p, '*')
    }else{
      star_p <- c(star_p, 'NS')
    }
  }

  if (show_p){
    show_label <- paste('p = ', p_values, sep='')
    # show_label <- sapply(p_values, function(x) substitute(paste(italic('p'), ' = ', x), list(x=x)))
  }else{
    show_label <- star_p
  }

    melt_df$type2 <- factor(melt_df$type2, levels = c("cCRE", "dELS", "pELS", "PLS", "DNase.H3K4me3", "CTCF.only"))
    g <- ggplot() +
    geom_bar(data=melt_df,aes(x=type2, y=mean, fill=type1), stat="identity", position = position_dodge(0.6)) +
    geom_errorbar(data=melt_df, aes(x=type2, ymin=mean-sd, ymax=mean+sd, color=type1, alpha=type1), stat='identity', position = position_dodge(0.6), width=0.2) +
      scale_color_manual(values = c('white', 'black')) +
      scale_alpha_manual(values = c(0, 0.9)) +
      scale_fill_igv() +
    guides(fill=guide_legend(NULL), color=guide_legend(NULL), alpha="none") +
    geom_segment(data=segments, aes(x=x1, y=y1, xend=x2, yend=y2)) +
    annotate('text', x=1:6,y=text_loc_y,label=show_label) +
    labs(x = '', y = 'Frequency', title="cCRE enrichment") +
    theme_bw() +
    theme(axis.title = element_text(size=14), axis.text = element_text(size=9)) +
    theme(plot.title = element_text(hjust = 0.5, size=16))

    return(g)
}





