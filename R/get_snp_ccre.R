#' Return the cCREs (candidate cis-regulatory elements) that the SNPs located in.
#'
#' @param snp Required. The chromosome coordinate or rsID of a SNP, eg. "chr22:19712094"
#' @param input_type Optional. The type of the input SNP. "rsID", hg19" or "hg38" can be selected. Default is "rsID".
#' @param output_assembly Optional. The UCSC Assembly versions of chromosome coordinates of cCREs. "hg19" or "hg38" can be selected. Default is "hg19".
#'
#' @return A data.frame of overlapped cCREs.
#' @export
#'
#' @examples
#' cCRE_intersect <- get_snp_ccre(snp="rs1059196")


get_snp_ccre <- function(snp, input_type="rsID", output_assembly="hg19"){
  if (output_assembly == "hg19"){
    cCRE_data <- cCRE_data_hg19
  }else{
    cCRE_data <- cCRE_data_hg38
  }
  colnames(cCRE_data) <- c("chrom", "start", "end", "accession", "cCRE_accession", "type")

  if (output_assembly=="hg19"){
    if (input_type=="hg38"){
      snp <- hg38tohg19(snp)
    }else if(input_type=="rsID"){
      snp <- get_snp_loc(snp, assembly = "hg19")
    }
  }else{
    if (input_type=="hg19"){
      snp <- hg19tohg38(snp)
    }else if(input_type=="rsID"){
      snp <- get_snp_loc(snp, assembly = "hg38")
    }
  }

  chrom <- str_split_fixed(snp, ":", 2)[1]
  pos <- as.numeric(str_split_fixed(snp, ":", 2)[2])
  pos_bed <- pos - 1

  intersect_cCREs <- cCRE_data[(cCRE_data$chrom==chrom)&(cCRE_data$start<=pos_bed)&(cCRE_data$end>=pos_bed),]
  return(intersect_cCREs)
}

