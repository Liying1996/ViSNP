#' Convert hg38 to hg19 coordinate.
#'
#' @param snp_loc The hg38 chromosome coordinate.
#'
#' @return A character. eg. "chr22:19724571".
#' @export
#'
#' @examples
#' hg38tohg19('chr22:19724571')

hg19tohg38 <- function(snp_loc){
  chrom <- str_split_fixed(snp_loc, ":", 2)[1]
  pos <- as.numeric(str_split_fixed(snp_loc, ":", 2)[2])
  # chain <- import.chain("data/hg19ToHg38.over.chain")
  chain_path <- system.file("data/hg19ToHg38.over.chain", package = "ViSNP")
  chain <- import.chain(chain_path)
  hg19.ver <- GRanges(seqnames=chrom, ranges=IRanges(start=pos, width = 1))
  hg38.ver <- liftOver(hg19.ver, chain)
  snp_hg38_loc <- as.character(hg38.ver@unlistData)
  return(snp_hg38_loc)
}
