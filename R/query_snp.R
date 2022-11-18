#' Return information of a SNP
#'
#' @param snp Required.
#' @param input_type Optional. The assembly version of the input SNP. "rsID", "hg19" and "hg38" can be selected. Default is "rsID".
#' @param eqtl_tissue Optional. Tissue ID of the tissue of interest. Default is "Whole_Blood".
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' info_table <- query_snp("chr1:1014863 ", input_type="hg38")



# Get message from 3D SNP2.0 database
query_snp <- function(snp, input_type="rsID", eqtl_tissue="Whole_Blood"){ # or can select output directly | version can set hg38 as well
  if  (!input_type %in% c("rsID", "hg19", "hg38")){
    return(message("Please select 1 input type from rsID, hg19 or hg38!"))
  }

  print("Get information from 3D-SNP v2 ...")

  if (input_type=="hg38"){
    snp <- get_snp_rsID(snp, assembly = "hg38")
  }else if(input_type=="hg19"){
    snp <- get_snp_rsID(snp, assembly = "hg19")
  }


  url <- paste("https://www.omic.tech/3dsnpv2/api.do?id=", snp, "&format=json&type=basic,motif,ccre,eqtl,physcores,3dsnp,3dgene,tfbs", sep="")
  res = GET(url, config = httr::config(ssl_verifypeer = FALSE)) # 3D-SNP (Location : hg19)

  res_data = fromJSON(rawToChar(res$content))
  rsID <- snp
  location <- res_data$disp_chr_pos[!is.na(res_data$disp_chr_pos)]
  chrom <- str_split_fixed(location, ":", 2)[1]
  position_hg19 <- as.numeric(str_split_fixed(location, ":", 2)[2]) + 1
  location_hg19 <- paste(chrom, position_hg19, sep = ":")
  location_hg38 <- get_snp_loc(snp, assembly = "hg38")


  ref_allele <- res_data$alleles_Ref[!is.na(res_data$alleles_Ref)]
  alt_allele <- res_data$alleles_Alt[!is.na(res_data$alleles_Alt)]

  # AF
  EAS_AF <- as.numeric(res_data$alleles_EAS[1])
  AMR_AF <- as.numeric(res_data$alleles_AMR[1])
  AFR_AF <- as.numeric(res_data$alleles_AFR[1])
  EUR_AF <- as.numeric(res_data$alleles_EUR[1])
  SAS_AF <- as.numeric(res_data$alleles_SAS[1])
  AF_info <- paste(EAS_AF,"(EAS),",AMR_AF, "(AMR),", AFR_AF, "(AFR),", EUR_AF, "(EUR),", SAS_AF, "(SAS)",sep="")


  linear_closest_gene <- res_data$data_gene[[1]]$name
  linear_closest_gene_id <- res_data$data_gene[[1]]$id
  linear_closest_gene_type <- res_data$data_gene[[1]]$loc

  linear_closest_gene_merged <- c()
  for (i in 1:length(linear_closest_gene)){
    new <- paste(linear_closest_gene[i], "(", linear_closest_gene_id[i], ")", sep="")
    linear_closest_gene_merged <- c(linear_closest_gene_merged, new)
  }

  linear_closest_gene_type <- paste(linear_closest_gene_type, collapse = ",")
  # linear_closest_gene_description <- res_data$data_gene[[1]]$description

  if (is.null(linear_closest_gene)){
    linear_closest_gene <- "-"
    linear_closest_gene_id <- "-"
    linear_closest_gene_type <- "-"
    linear_closest_gene_description <- "-"
  }

  loop_gene_col <- res_data$data_loop_gene[!is.null(res_data$data_loop_gene)]
  for (i in 1:length(loop_gene_col)){
    if (!is.null(loop_gene_col[[i]])){
      loop_gene_info <- loop_gene_col[[i]]
      break
    }
  }
  loop_gene <- paste(unique(loop_gene_info$gene)[unique(loop_gene_info$gene)!=""], collapse =",")
  if (loop_gene==""){
    loop_gene <- "-"
  }

  loop_snp_col <- res_data$data_loop_snp
  for (i in 1:length(loop_snp_col)){
    if (!is.null(loop_snp_col[[i]])){
      loop_snp_info <- loop_snp_col[[i]]
      break
    }
  }
  loop_snps <- unique(loop_snp_info$SNP_B)[unique(loop_snp_info$SNP_B)!=""]
  if (length(loop_snps)!=0){
    loop_snps <- loop_snps[order(loop_snps)]
    if (length(loop_snps)>=5){
      loop_snps <- paste(loop_snps[1:5], collapse =",") # Only shows 5 SNPs-B
    }else{
      loop_snps <- paste(loop_snps, collapse =",")
    }
  }else{
    loop_snps <- "-"
  }

  if (loop_snps==""){
    loop_snps <- "-"
  }

  tfbs_col <- res_data$data_tfbs
  tfbs <- "-"
  if (!is.null(tfbs_col)){
    for (i in 1:length(tfbs_col)){
      if (!is.null(tfbs_col[[i]])){
        tfbs_info <- tfbs_col[[i]]
        break
      }
    }
    tfbs <- unique(tfbs_info$gene)
    tfbs <- paste(tfbs, collapse = ",")
  }

  print("Get information from SCREEN and GWAS catalog...")
  cCRE_table <- get_snp_ccre(snp)
  if (nrow(cCRE_table)==0){
    cCRE_info <- "-"
  }else{
    cCRE_info <- paste(cCRE_table$type, " (", cCRE_table$SCREEN_accession, ")", sep="")
  }


  gwas_table <- get_snp_gwas(snp)
  if (nrow(gwas_table)==0){
    gwas_info <- "-"
  }else{
    gwas_info <- unique(gwas_table$DISEASE.TRAIT)
    gwas_info <- paste(gwas_info, collapse = "; ")
  }

  # eQTL GTEx api
  print("Get information from GTEx...")
  gtex_data <- get_snp_eqtl(snp, eqtl_tissue = eqtl_tissue)
  if (length(gtex_data)==0){
    gtex_info <- "-"
  }else{
    gtex_info <- paste(gtex_data$geneSymbol, collapse = ",")
  }



  cat(paste("rsID:", rsID, "\n"))
  # cat(paste("Assembly:", assembly, "\n"))
  cat(paste("chrom:", chrom, "\n"))
  cat(paste("Location(hg19):", location_hg19, "\n"))
  cat(paste("Location(hg38):", location_hg38, "\n"))
  cat(paste("Ref allele:", ref_allele, "\n"))
  cat(paste("Alt allele:", alt_allele, "\n"))
  cat(paste("AF:", AF_info, "\n"))
  if (length(linear_closest_gene)<=1){
    if (linear_closest_gene!="-"){
      cat(paste("Linear closest gene: ", linear_closest_gene_merged, "\n", sep=""))
    }else{
      cat(paste("Linear closest gene:", linear_closest_gene,"\n", sep = ""))
    }
  }else{
    cat(paste("Linear closest gene:", paste(linear_closest_gene_merged, collapse = ","), "\n", sep=""))
  }
  cat(paste("Linear closest gene type:", linear_closest_gene_type, "\n"))
  # cat(paste("Linear closest gene description:", linear_closest_gene_description, "\n"))
  cat(paste("3D Loop Genes:", loop_gene, "\n"))
  cat(paste("3D-Interacting SNPs:", loop_snps, "\n"))
  cat(paste("TFBS:", tfbs, "\n"))
  cat(paste("GTEx eQTLs:"), gtex_info, "\n")
  cat(paste("cCRE:", cCRE_info, "\n"))
  cat(paste("GWAS:", gwas_info, "\n"))


  info_table <- data.frame(Info=c("rsID", "Chromosome", "Location(hg19)", "Location(hg38)", "Ref allele", "Alt allele","AF",  "Linear closest gene","Linear closest gene type", "3D-Loop Genes", "3D-Interacting SNPs", "TFBS","eQTLs", "cCRE", "GWAS"), Value=c(rsID, chrom, location_hg19, location_hg38, ref_allele, alt_allele, AF_info, paste(linear_closest_gene_merged, collapse = ","), linear_closest_gene_type, loop_gene, loop_snps, tfbs, gtex_info, cCRE_info, gwas_info))

  return(info_table)

}
