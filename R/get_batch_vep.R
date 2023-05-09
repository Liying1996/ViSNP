#' Return annotation by VEP
#'
#' @param snps Required.
#' @param input_type Optional. The type of the input SNP. "rsID", hg19" or "hg38" can be selected. Default is "rsID".
#' @return A data.frame of VEP annotation
#' @export
#'
#' @examples
#' vep_anno <- get_batch_vep(c("rs56116432", "rs10040658"), input_type = "rsID") # List of annotations

get_batch_vep <- function(snps, input_type='rsID'){
  snps <- unique(str_split_fixed(snps, ',', 2)[,1]) # eg. "rs11074135,COSV61234901"

  if (input_type=="hg38"){
    new_snps <- c()
    for (snp in snps){
      new_snps <- c(new_snps, get_snp_rsID(snp, assembly = "hg38"))
    }
  }else if(input_type=="hg19"){
    new_snps <- c()
    for (snp in snps){
      new_snps <-  c(new_snps, get_snp_rsID(snp, assembly = "hg19"))
    }
  }

  snps <- new_snps
  n <- length(snps)
  mod <- n %% 200
  if (mod == 0){
    fold <- max(floor(n/200)-1, 0)
  }else{
    fold <- floor(n/200)
  }
  res_list <- list()
  for (i in 0:fold){
    range_min <- 200*i + 1
    if (i != fold){
      range_max <- 200*(i+1)
    }else{
      range_max <- n
    }
    snps_curr <- snps[range_min:range_max]
    res <- POST("https://rest.ensembl.org/vep/human/id", content_type("application/json"), accept("application/json"), body = paste('{ "ids" : [ "', paste(snps_curr, collapse = '","'), '" ] }', sep=''))
    res_data = fromJSON(rawToChar(res$content))
    res_list[[i+1]] <- res_data
}
  anno_data <- do.call("rbind", res_list)
  return(anno_data)

}


