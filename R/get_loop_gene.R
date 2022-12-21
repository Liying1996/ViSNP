#' Return 3D-interacting genes of a SNP
#'
#' @param snp Required.
#' @param input_type Optional. The type of the input SNP. "rsID", hg19" or "hg38" can be selected. Default is "rsID".
#' @param output_type Optional. The type of the output matrix. "core" or "full" can be selected. The output of "core" is subset of "full".Default is "core".
#'
#' @return A data.frame of 3D-interacting genes.
#' @export
#'
#' @examples
#' loop_genes <- get_loop_gene("rs1059196")


get_loop_gene <- function(snp, input_type="rsID", output_type="core"){
  if  (!input_type %in% c("rsID", "hg19", "hg38")){
    return(message("Please select 1 input type from rsID, hg19 and hg38!"))
  }
  if  (!output_type %in% c("core", "full")){
    return(message("Please select 1 output type from core and full!"))
  }

  if (input_type=="hg38"){
    snp <- get_snp_rsID(snp, assembly = "hg38")
  }else if(input_type=="hg19"){
    snp <- get_snp_rsID(snp, assembly = "hg19")
  }

  url <- paste("https://www.omic.tech/3dsnpv2/api.do?id=", snp, "&format=json&type=3dgene", sep="")
  res = GET(url, config = httr::config(ssl_verifypeer = FALSE))
  res_data = fromJSON(rawToChar(res$content))

  loop_gene_col <- res_data$data_loop_gene[!is.null(res_data$data_loop_gene)]
  for (i in 1:length(loop_gene_col)){
    if (!is.null(loop_gene_col[[i]])){
      loop_gene_info <- loop_gene_col[[i]]
      break
    }
  }
  if (nrow(loop_gene_info)==0){
    return(data.frame())
  }

  if (output_type == "full"){
    return(loop_gene_info)
  }else{
    loop_gene_info_subset <- loop_gene_info[c("start", "end", "cell", "gene", "tissue",  "type_loop")]
    return(loop_gene_info_subset)
  }
}
