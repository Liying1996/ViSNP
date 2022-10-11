#' Return a table and a line chart of the physcores of +/-10bp of SNPs
#'
#' @param snps Required.
#' @param colors Optional. Colors of lines. If not specified, default is to use the colors of `rainbow()`.
#' @param input_type Optional. The assembly version of the input SNP. "rsID", "hg19" and "hg38" can be selected. Default is "rsID".
#'
#' @return A data.frame of physcores and a line chart.
#' @export
#'
#' @examples
#' snps <- c("rs2769466", "rs4661232")
#' physcores_table <- plot_physcores(snps, colors = c("blue", "orange"))


plot_physcores <- function(snps, colors="default", input_type="rsID"){ # 1 snp is okay
  if (colors[1]!="default" & length(colors)!=length(snps)){
    return(message("Please provide the same number of colors as the snps!"))
  }
  cat("Get PhyloP scores from 3D-SNP v2 ...\n")

  if (input_type=="hg19"){
    new_snps <- c()
    for (s in snps){
      s <- get_snp_rsID(s, assembly = "hg19")
      new_snps <- c(new_snps, s)
    }
    snps <- new_snps
  }else if (input_type=="hg38"){
    new_snps <- c()
    for (s in snps){
      s <- get_snp_rsID(s, assembly = "hg38")
      new_snps <- c(new_snps, s)
    }
    snps <- new_snps
  }

  physcores <- c()
  snps_list <- c()
  for (snp in snps){
    url <- paste("https://www.omic.tech/3dsnpv2/api.do?id=", snp, "&format=json&type=basic,physcores", sep="")
    res = GET(url, config = httr::config(ssl_verifypeer = FALSE))
    res_data = fromJSON(rawToChar(res$content))
    for (p in res_data$physcores_update){
      if (length(p)!= 0){
        single_physcores <- p
        break
      }
    }
    physcores <- c(physcores, single_physcores)
    snps_list <- c(snps_list, rep(snp, 21))
  }
  phylop_df <- data.frame(position=rep(c(-10:10), length(snps)), physcore=physcores, snp = snps_list)

  if (colors[1]=="default"){
    colors <- rainbow(length(snps))
  }
    g <- ggplot(phylop_df, aes(x = position, y = physcore, color=snp, shape=snp)) +
     geom_line() + geom_point() +
     geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
     theme_snp() +
     scale_color_manual(values = colors) +
     guides(color = guide_legend(title=NULL), shape=guide_legend(title=NULL)) +
     labs(x = "Relative Position", y = "PhyloP score")
    print(g)

    return(phylop_df)

}

