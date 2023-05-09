#' Return the network.
#'
#' @param snp Required.
#' @param input_type Optional. The assembly version of the input SNP. "rsID", "hg19" and "hg38" can be selected. Default is "rsID".
#' @param cells Optional. The output cells. Default is "all".
#' @param pop Required. Population.
#' @param tissue_type Optional. The tissue of interest (for eQTLs). Default is "Whole_Blood".. Default is "Whole_Blood".
#' @param colors Optional. Vectors with 6 colors. "default" is `brewer.pal(6, "Set1")`.
#'
#' @return A Graph.
#' @export
#'
#' @examples
#' plot_snp_network(snp="rs1059196")


plot_snp_network <- function(snp, input_type = "rsID", cells = "all", pop, tissue_type="Whole_Blood", colors="default"){

  if  (!input_type %in% c("rsID", "hg19", "hg38")){
    return(message("Please select 1 input type from rsID, hg19 or hg38!"))
  }

  if (input_type=="hg38"){
    snp <- get_snp_rsID(snp, assembly = "hg38")
  }else if(input_type=="hg19"){
    snp <- get_snp_rsID(snp, assembly = "hg19")
  }

  if (colors[1] == "default"){
    colors <- brewer.pal(7, "Set1")
  }

  ccre_info <- get_snp_ccre(snp)
  gwas_info <- get_snp_gwas(snp)
  loop_genes <- get_loop_gene(snp)
  loop_genes <- loop_genes[complete.cases(loop_genes),]
  loop_snps <- get_loop_snp(snp)
  loop_snps <- loop_snps[loop_snps$type==pop,]
  loop_snps <- unique(loop_snps[complete.cases(loop_snps), 3:ncol(loop_snps)])
  if (cells[1] !="all"){
    loop_snps <- loop_snps[loop_snps$cell %in% cells,]
  }


  eqtl_info <- get_snp_eqtl(snp, eqtl_tissue = tissue_type)

  nodes_gwasid <- gwas_info$PUBMEDID
  nodes_gwas_desease <- gwas_info$DISEASE.TRAIT
  nodes_ccre <- ccre_info$SCREEN_accession
  nodes_ccre_type <- ccre_info$type
  nodes_3dsnp_cell <- loop_snps$cell
  nodes_3dsnp_snp <- loop_snps$SNP_B
  nodes_3dgenes_cell <- loop_genes$cell
  nodes_3dgenes_gene <- loop_genes$gene
  nodes_eqtl_tissue <- eqtl_info$tissueSiteDetailId
  nodes_eqtl_gene <- eqtl_info$geneSymbol
  nodes_eqtl_gene2 <- c()
  for (x in eqtl_info$geneSymbol){
    if (!x %in% nodes_3dgenes_gene){
      nodes_eqtl_gene2 <- c(nodes_eqtl_gene2, x)
    }
  }

  gwas_cate <- rep("GWAS", length(unique(c(nodes_gwasid, nodes_gwas_desease))))
  ccre_cate <- rep("cCREs", length(unique(c(nodes_ccre, nodes_ccre_type))))
  loop_cate <- rep("3D-loop", length(unique(c(nodes_3dsnp_cell, nodes_3dgenes_cell))))
  loop_snp_cate <- rep("3D-loop snp", length(unique(nodes_3dsnp_snp)))
  loop_gene_cate <- rep("3D-loop gene", length(unique(nodes_3dgenes_gene)))
  eqtl_cate <- rep("eQTL", length(unique(c(nodes_eqtl_tissue, nodes_eqtl_gene))))
  eqtl_cate2 <- rep("eQTL", length(unique(c(nodes_eqtl_tissue, nodes_eqtl_gene2))))

  # Label and Group
  nodes_snp <- data.frame(Label = c(snp, unique(c(nodes_gwasid, nodes_gwas_desease)), unique(c(nodes_ccre, nodes_ccre_type)), unique(c(nodes_3dsnp_cell, nodes_3dgenes_cell)), unique(nodes_3dsnp_snp), unique(nodes_3dgenes_gene), unique(c(nodes_eqtl_tissue, nodes_eqtl_gene2))), Category = c("SNP", gwas_cate, ccre_cate, loop_cate, loop_snp_cate, loop_gene_cate, eqtl_cate2))

  # from-to-weight

  directions_to <- c(unique(nodes_gwasid), unique(nodes_ccre_type), unique(c(nodes_3dsnp_cell, nodes_3dgenes_cell)), unique(nodes_eqtl_tissue))
  directions_from <- rep(snp, length(directions_to))

  edges_snp <- data.frame(from=c(directions_from, nodes_gwasid, nodes_ccre_type, nodes_3dsnp_cell, nodes_3dgenes_cell, nodes_eqtl_tissue), to=c(directions_to, nodes_gwas_desease, nodes_ccre, nodes_3dsnp_snp, nodes_3dgenes_gene, nodes_eqtl_gene), weight=rep(1, length(directions_from) + length(c(nodes_gwasid, nodes_ccre, nodes_3dsnp_cell, nodes_3dgenes_cell, nodes_eqtl_tissue))))
edges_snp <- edges_snp[!duplicated(edges_snp),]

  # Draw

  net <- graph_from_data_frame(
    d=edges_snp,vertices=nodes_snp,
    directed=TRUE)

  deg <- igraph::degree(net)
  # vcolor <- c("orange","red","lightblue","purple","yellow", "green")
  vcolor <- colors

  # Nodes colors
  V(net)$color <- vcolor[factor(V(net)$Category)]

  # Edges width
  E(net)$width <- E(net)$weight

  lay <- do.call("layout_with_kk", list(net))
  plot(net, vertex.size=10+deg*2.5, vertex.color = adjustcolor(V(net)$color, alpha.f = .7), edge.arrow.size=.25, edge.curved=.1, vertex.label.cex=.75, vertex.label.dist=0.1, layout=lay)

}
