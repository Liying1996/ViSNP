#' Convert the chromosome coordinate of a SNP to rsID
#'
#' @param snp_loc Required. The chromosome coordinate of a SNP, eg. "chr22:19712094"
#' @param assembly Optional. UCSC Assembly versions. "hg19" or "hg38" can be selected. Default is "hg19".
#'
#' @return The rsID of the SNP.
#' @export
#'
#' @examples
#' get_snp_rsID("chr22:19712094") # hg19
#'


get_snp_rsID <- function(snp_loc, assembly="hg19"){ # input: location output: rsID
  if (assembly=="hg38"){
    snp_loc <- hg38tohg19(snp_loc)
  }
  chrom <- str_split_fixed(snp_loc, ":", 2)[1]
  pos <- as.numeric(str_split_fixed(snp_loc, ":", 2)[2])
  url <- paste("https://www.omic.tech/3dsnpv2/api.do?position=", pos-2, "-", pos, "&chrom=",chrom ,"&format=json&type=basic", sep="")
  res = GET(url, config = httr::config(ssl_verifypeer = FALSE)) # 3D-SNP (Location : hg19)
  res_data = fromJSON(rawToChar(res$content))
  rsID <- res_data$id[1]
  return(rsID)
}


