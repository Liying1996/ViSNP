#' Return the genomic locations of gene symbols
#'
#' @param genes Required. Gene Symbol(s), eg. "GP1BB"
#' @param assembly Optional. Output UCSC Assembly versions. "hg19" or "hg38" can be selected. Default is "hg19".
#'
#' @return A data.frame including the genomic location.
#' @export
#'
#' @examples
#' get_genes_loc(c('GP1BB', 'TBX1'))
#'

get_genes_loc <- function(genes, assembly="hg19"){
  if (!assembly %in% c("hg19", "hg38")){
    return(message("Please select hg19 or hg38 version!"))
  }
  if (assembly == "hg19"){
    genome_version = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  } else {
    genome_version = useMart(biomart="ensembl", host="https://www.ensembl.org", dataset="hsapiens_gene_ensembl")
  }
  ensembl = useDataset("hsapiens_gene_ensembl",mart=genome_version)
  info <- getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version", "external_gene_name",'chromosome_name','start_position','end_position'), mart = ensembl)
  genes <- as.data.frame(genes)
  colnames(genes)[1]<-"external_gene_name"
  genes_loc_df <- merge(genes, info, by= 'external_gene_name')
  colnames(genes_loc_df)[1] <- "gene_symbol"
  return(genes_loc_df)
}
