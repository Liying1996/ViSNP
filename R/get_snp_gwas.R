#' Return the associated phenotypes and related research (from GWAS catalog) of the SNP.
#' @param snp Required.
#' @param input_type Optional. The type of the input SNP. "rsID", hg19" or "hg38" can be selected. Default is "rsID".
#'
#' @return A data.frame of associated GWAS phenpotypes.
#' @export
#'
#' @examples
#' gwas_info <- get_snp_gwas("rs1891906")



get_snp_gwas <- function(snp, input_type="rsID"){ # rsID, hg19, hg38
  gwas_data <- gwas_data

  if (input_type=="hg19"){
    snp <- hg19tohg38(snp)
  }else if(input_type=="rsID"){
    snp <- get_snp_loc(snp, assembly = "hg38")
  }
  chrom <- str_split_fixed(snp, ":", 2)[1]
  pos <- str_split_fixed(snp, ":", 2)[2]
  gwas_info <- gwas_data[(gwas_data$CHR_ID==gsub("chr","", chrom)) & (gwas_data$CHR_POS==pos), ]
  return(gwas_info)
}
