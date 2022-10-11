#' Return GTEx significant single tissue eQTLs.
#'
#' Using the \href{https://gtexportal.org/home/api-docs/index.html#!/association/singleTissueEqtl}{GTEx api} to generate the information.
#'
#' @param snp Required.
#' @param input_type Optional. The type of the input SNP. "rsID", "hg19" or "hg38" can be selected. Default is "rsID".
#' @param eqtl_tissue Optional. Tissue ID of the tissue of interest. Default is "Whole_Blood".
#'
#' @return A data.frame of associated phenotypes and reserach of the SNP.
#' @export
#'
#' @examples
#' gtex_data <- get_snp_eqtl(snp, eqtl_tissue = "Whole_Blood")


get_snp_eqtl <- function(snp, input_type="rsID", eqtl_tissue="Whole_Blood"){
  if  (!input_type %in% c("rsID", "hg19", "hg38")){
    return(message("Please select from rsID, hg19 or hg38!"))
  }
  if (input_type=="hg38"){
    snp <- get_snp_rsID(snp, assembly = "hg38")
  }else if(input_type=="hg19"){
    snp <- get_snp_rsID(snp, assembly = "hg19")
  }

  gtex_url <- paste("https://gtexportal.org/rest/v1/association/singleTissueEqtl?format=json&snpId=", snp,"&tissueSiteDetailId=", eqtl_tissue, "&datasetId=gtex_v8", sep="")
  gtex_page = GET(gtex_url) # 3D-SNP (Location : hg19)
  gtex_data = fromJSON(rawToChar(gtex_page$content))
  return(gtex_data$singleTissueEqtl)
}


