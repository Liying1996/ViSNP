#' Return information of a batch of SNPs (Annotated by query_snps())
#'
#' @param snps Required.
#' @param input_type Optional. The assembly version of the input SNP. "rsID", "hg19" and "hg38" can be selected. Default is "rsID".
#' @param eqtl_tissue Optional. Tissue ID of the tissue of interest. Default is "Whole_Blood".
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' batch_info_table <- get_batch_SNP_info(snps = c("rs1891906" "rs10"), input_type="rsID", eqtl_tissue="Whole_Blood")


get_batch_SNP_info <- function(snps, input_type="rsID", eqtl_tissue="Whole_Blood"){
  batch_df <- as.data.frame(matrix(nrow=15, ncol=length(snps)))
  index = 1
  for (snp in snps){
    info_table <- query_snp(snp, input_type=input_type, eqtl_tissue=eqtl_tissue)
    batch_df[, index] <- info_table$Value
    index <- index + 1
  }
  new_batch_df <- as.data.frame(t(batch_df))
  colnames(new_batch_df) <- info_table$Info
  return(new_batch_df)
}
