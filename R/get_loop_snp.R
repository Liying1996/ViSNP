#' Return 3D-interacting SNPs of a SNP
#'
#' @param snp Required.
#' @param input_type Optional. The type of the input SNP. "rsID", hg19" or "hg38" can be selected. Default is "rsID".

#'
#' @return A data.frame of 3D-interacting SNPs
#' @export
#'
#' @examples
#' loop_snps <- get_loop_snp("rs1059196")

get_loop_snp <- function(snp, input_type="rsID"){
  if (input_type=="hg38"){
    snp <- get_snp_rsID(snp, assembly = "hg38")
  }else if(input_type=="hg19"){
    snp <- get_snp_rsID(snp, assembly = "hg19")
  }

  url <- paste("https://www.omic.tech/3dsnpv2/api.do?id=", snp, "&format=json&type=3dsnp", sep="")
  res = GET(url, config = httr::config(ssl_verifypeer = FALSE))
  res_data = fromJSON(rawToChar(res$content))

  loop_snp_col <- res_data$data_loop_snp
  for (i in 1:length(loop_snp_col)){
    if (!is.null(loop_snp_col[[i]])){
      loop_snp_info <- loop_snp_col[[i]]
      break
    }
  }

  if (nrow(loop_snp_info)==0){
    return(data.frame())
  }
  return(loop_snp_info)

}
