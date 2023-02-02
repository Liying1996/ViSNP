#' Generates pairwise linkage disequilibrium statistics.
#' It should be noted that LDlinkR is called here, so a Personal Access Token must be required. Please refer to \href{https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html}{LDlinR} for details.
#' @param snps Required. between 1 - 10 variants, using an rsID or chromosome coordinate.
#' @param pop Optional. Default = “CEU”. list_pop() return all the options.
#' @param r2d Optional. Default = “r2”. Either “r2” for LD R2 (R-squared) or “d” for LD D’ can be selected.
#' @param assembly Optional. Output UCSC Assembly versions. "hg19" or "hg38" can be selected. Default is "hg19".
#' @param plot_heatmap Optional. Whether to output the heatmap of the LD block. Default is TRUE.
#'
#' @return A data.frame of pairwise linkage disequilibrium statistics. Desired output can be based on estimates of R2 or D’.
#' @export
#'
#' @examples
#' ld_block <- get_LD(c("chr7:24962419", "rs114", "rs127", "rs7805287", "rs60676332", "rs10239961"), assembly = "hg19")
#'


get_LD <- function(snps, pop="CEU", r2d="r2", assembly = "hg19", plot_heatmap=TRUE){
  if (assembly=="hg19"){
    g <- "grch37"
  }else{
    g <- "grch38"
  }
  ld_mat <- LDmatrix(snps = snps,
           pop = pop,
           r2d = r2d,
           token = Sys.getenv("LDLINK_TOKEN"),
           genome_build = g
          )
  rownames(ld_mat) <- ld_mat$RS_number
  ld_mat2 <- as.matrix(ld_mat[,-1])

  if (plot_heatmap){
    key_range <- seq(0, 1,by=0.2)
    colors = colorRampPalette(colors = c("white", "#FF4500"))(length(key_range))
    p <- pheatmap(ld_mat2, scale = "none", col = colors, display_numbers = T, angle_col = "45", cellwidth = 45, cellheight  = 45, breaks = key_range, legend_breaks = key_range, legend_labels = key_range, border_color ="white")
    print(p)
  }
    return(ld_mat2)
}
