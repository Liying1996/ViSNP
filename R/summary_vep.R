#' Return barplots of affected genes of SNPs.
#'
#' @param data Required. The annotation results from VEP.
#' @param plot_type Optional."gene", "snp", "all" and "merged" can be selected. Default is "rsID".
#' @param show_num Optional. The number of affected genes shown on the plot. Default is 7.
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#'

summary_vep <- function(vep_file, output="summmary_vep.txt"){
  data <- read.table(vep_file, header = F)
  colnames(data) <- c("Uploaded_variation", "Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence",	"cDNA_position",	"CDS_position",	"Protein_position",	"Amino_acids",	"Codons",	"Existing_variation",	"Extra")


}


