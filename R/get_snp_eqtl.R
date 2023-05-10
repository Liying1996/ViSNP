#' Return GTEx significant single tissue eQTLs.
#'
#' Using the \href{https://gtexportal.org/api/v2/association/singleTissueEqtl}{GTEx api} to generate the information.
#'
#' @param snp Required.
#' @param input_type Optional. The type of the input SNP. "rsID", "hg19" or "hg38" can be selected. Default is "rsID".
#' @param eqtl_tissue Optional. Tissue ID of the tissue of interest. Default is "Whole_Blood".
#'
#' @return A data.frame of associated phenotypes and reserach of the SNP.
#' @export
#'
#' @examples
#' gtex_data <- get_snp_eqtl(snp="rs10040658", eqtl_tissue = "Whole_Blood")


get_snp_eqtl <- function(snp, input_type="rsID", eqtl_tissue="Whole_Blood"){
  if  (!input_type %in% c("rsID", "hg19", "hg38")){
    return(message("Please select from rsID, hg19 or hg38!"))
  }
  if (input_type=="hg38"){
    snp <- get_snp_rsID(snp, assembly = "hg38")
  }else if(input_type=="hg19"){
    snp <- get_snp_rsID(snp, assembly = "hg19")
  }
  # GTEx V1
  # gtex_url <- paste("https://gtexportal.org/rest/v1/association/singleTissueEqtl?format=json&snpId=", snp,"&tissueSiteDetailId=", eqtl_tissue, "&datasetId=gtex_v8", sep="")
  # gtex_page = GET(gtex_url) # 3D-SNP (Location : hg19)
  # gtex_data = fromJSON(rawToChar(gtex_page$content))
  # return(gtex_data$singleTissueEqtl)

  # GTEx V2
  url <- paste("https://www.omic.tech/3dsnpv2/api.do?id=", snp, "&format=json&type=basic", sep="")
  res = GET(url, config = httr::config(ssl_verifypeer = FALSE))

  res_data = fromJSON(rawToChar(res$content))
  location <- res_data$disp_chr_pos[!is.na(res_data$disp_chr_pos)]
  location_hg38 <- get_snp_loc(snp, assembly = "hg38")
  ref_allele <- res_data$alleles_Ref[!is.na(res_data$alleles_Ref)]
  alt_allele <- res_data$alleles_Alt[!is.na(res_data$alleles_Alt)]

  variant_id <- paste(gsub(':', '_', location_hg38), ref_allele, alt_allele, 'b38', sep='_')
  url_gtex <- "https://gtexportal.org/api/v2/association/singleTissueEqtl"
  params <- list(
    variantId = variant_id,
    datasetId = "gtex_v8",
    tissueSiteDetailId = eqtl_tissue
  )

  response <- GET(url_gtex, query = params)
  json <-  fromJSON(rawToChar(response$content))
  return(json$data)

}


