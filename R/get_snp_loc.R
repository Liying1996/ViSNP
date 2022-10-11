#' Convert SNP rsID to genomic locations.
#'
#' @param snp Required.
#' @param assembly Optional. The assembly version of the input SNP. hg19" and "hg38" can be selected. Default is "hg19".
#'
#' @return A character.
#' @export
#'
#' @examples
#' get_snp_loc("rs1059196")


get_snp_loc <- function(snp, assembly="hg19"){ # 1-based same as snpDB
  url <- paste("https://www.omic.tech/3dsnpv2/api.do?id=", snp, "&format=json&type=basic", sep="")
  res = GET(url, config = httr::config(ssl_verifypeer = FALSE))
  res_data = fromJSON(rawToChar(res$content))
  location <- res_data$disp_chr_pos[!is.na(res_data$disp_chr_pos)]
  chrom <- str_split_fixed(location, ":", 2)[1]
  pos <- as.numeric(str_split_fixed(location, ":", 2)[2]) + 1
  location <- paste(chrom, pos, sep=":")

  if (assembly=="hg38"){
    location <- hg19tohg38(location)
  }

  return(location)
}

