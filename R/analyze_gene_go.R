#' Return the go enrichment results of genes.
#'
#' @param snp_genes Required. The annotated genes of input SNPs.
#' @param output_type Optional. Show dotplot/barplot of GO enrichment results. Default is dotplot.
#'
#' @return A plot.
#' @export
#'
#' @examples
#' analyze_gene_go(snps_genes)


analyze_gene_go <- function(snp_genes, output_type="dotplot"){
  genes = bitr(snp_genes, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

  go <- enrichGO(gene = unique(genes$ENTREZID), OrgDb = org.Hs.eg.db,ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = FALSE, keyType = 'ENTREZID')

  if (nrow(go@result)==0){
    print('0 enriched terms found...')
    break
  }
  if (output_type == "dotplot"){
    g <- dotplot(go,showCategory=10, title="Go enrichment Dotplot") +
      scale_y_discrete(labels=function(x) str_wrap(x, width=20)) +
      theme(plot.title = element_text(hjust = 0.5, size=20))
  }else{
    g <- barplot(go,showCategory=7, drop=T,title  = "Go enrichment Barplot",font.size = 15) +
    scale_y_discrete(labels=function(x) str_wrap(x, width=20)) +
      theme(plot.title = element_text(hjust = 0.5, size=20))
  }
  return(g)
}
