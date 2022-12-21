#' Return the associated phenotypes and related research (from GWAS catalog) of the SNP.
#' @param snp Required.
#' @param input_type Optional. The type of the input SNP. "rsID", hg19" or "hg38" can be selected. Default is "rsID".
#' @param output_type Optional. The type of the output matrix. "core" or "full" can be selected. The output of "core" is subset of "full".Default is "core".
#'
#' @return A data.frame of associated GWAS phenpotypes.
#' @export
#'
#' @examples
#' gwas_info <- get_snp_gwas("rs1891906")



get_snp_gwas <- function(snp, input_type="rsID", output_type="core"){

  if  (!input_type %in% c("rsID", "hg19", "hg38")){
    return(message("Please select 1 input type from rsID, hg19 and hg38!"))
  }
  if  (!output_type %in% c("core", "full")){
    return(message("Please select 1 output type from core and full!"))
  }

  gwas_data <- gwas_data

  if (input_type=="hg19"){
    snp <- hg19tohg38(snp)
  }else if(input_type=="rsID"){
    snp <- get_snp_loc(snp, assembly = "hg38")
  }
  chrom <- str_split_fixed(snp, ":", 2)[1]
  pos <- str_split_fixed(snp, ":", 2)[2]
  gwas_info <- gwas_data[(gwas_data$CHR_ID==gsub("chr","", chrom)) & (gwas_data$CHR_POS==pos), ]

  if (output_type == "full"){
    return(gwas_info)
  }else{
    gwas_info_subset <- gwas_info[c("CHR_ID", "CHR_POS", "SNPS", "STRONGEST.SNP.RISK.ALLELE", "MAPPED_GENE", "PUBMEDID", "DISEASE.TRAIT", "STUDY")]
    return(gwas_info_subset)
  }
}
