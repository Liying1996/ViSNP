#' Return SNP-disrupted motifs
#'
#' @param snp Required.
#' @param input_type Optional. The type of the input SNP. "rsID", hg19" or "hg38" can be selected. Default is "rsID".

#'
#' @return A data.frame of disrupted motifs.
#' @export
#'
#' @examples
#' disrupted_motifs <- get_snp_motif("chr1:109676139", "hg38")


get_snp_motif <- function(snp, input_type="rsID"){ # 1-based same as SNPdb

  if (input_type=="hg38"){
    snp <- get_snp_rsID(snp, assembly = "hg38")
  }else if(input_type=="hg19"){
    snp <- get_snp_rsID(snp, assembly = "hg19")
  }

  url <- paste("https://www.omic.tech/3dsnpv2/api.do?id=", snp, "&format=json&type=motif", sep="")
  res = GET(url, config = httr::config(ssl_verifypeer = FALSE))
  res_data = fromJSON(rawToChar(res$content))

  n_motif <- str_count(res_data$motif, ";")
  motifs <- data.frame(t(data.frame(str_split_fixed(res_data$motif, ";", n_motif+1))))
  motifs <- data.frame(str_split_fixed(motifs[,1], ",", 5))
  colnames(motifs) <- c("Source", "Motif", "Strand", "Sequence", "Altered_Pos")

  return(motifs)
}
