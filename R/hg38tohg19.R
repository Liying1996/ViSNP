#' Convert hg38 to hg19 coordinate.
#'
#' @param snp_loc The hg38 chromosome coordinate.
#'
#' @return A character. eg. "chr22:19724571".
#' @export
#'
#' @examples
#' hg38tohg19('chr22:19724571')


hg38tohg19 <- function(snp_loc){
  chrom <- str_split_fixed(snp_loc, ":", 2)[1]
  pos <- as.numeric(str_split_fixed(snp_loc, ":", 2)[2])
  chain <- import.chain("data/hg38ToHg19.over.chain")
  hg38.ver <- GRanges(seqnames=chrom, ranges=IRanges(start=pos, width = 1))
  hg19.ver <- liftOver(hg38.ver, chain)
  snp_hg19_loc <- as.character(hg19.ver@unlistData)
  return(snp_hg19_loc)
}
