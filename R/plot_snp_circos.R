#' Return the circos plot.
#'
#' @param snp_info_table Required. The output data.frame of `query_snp()`.
#' @param output_assembly Optional. The version of chromosome coordinates of output circos plots. Default is "hg19".
#' @param window_size Optional. Window size around the SNP displayed. Default is 2e5.
#' @param circos_pos Optional. Circos plot coordinate position. You can choose 'Relative' or 'Absolute'. Relative position represents the window range (0-window_size) of the input SNP, while absolute position is the actual genomic coordinates. The default is 'Relative'.
#' @param savefile Optional. The filename of the circos plot.
#'
#' @return A circos plot.
#' @export
#'
#' @examples
#' info_table <- query_snp("rs1891906")
#' plot_snp_circos(snp_info_table=info_table, window_size = 1e5, savefile="~/circos.pdf")


plot_snp_circos <- function(snp_info_table, output_assembly = "hg19", window_size=2e5, circos_pos="Relative", savefile="circos.pdf"){

  snp <- snp_info_table[snp_info_table$Info=="rsID", 2]
  chrom <- snp_info_table[snp_info_table$Info=="Chromosome", 2]
  window_size <- as.numeric(window_size)

  if (output_assembly=="hg19"){
    location <- snp_info_table[snp_info_table$Info=="Location(hg19)", 2]
    chrom_info <- chrom_info_hg19
  }else{
    location <- snp_info_table[snp_info_table$Info=="Location(hg38)", 2]
    chrom_info <- chrom_info_hg38
  }

  chrom_max_len <- chrom_info[chrom_info$V1==chrom, 2]
  pos <-  as.numeric(str_split_fixed(location, ":", 2)[2])

  # cat("Get the genomic locations of 3D-loop genes...")
  loop_genes <- snp_info_table[snp_info_table$Info=="3D-Loop Genes",2]
  n <- str_count(loop_genes, ",")
  loop_genes <- as.vector(str_split_fixed(loop_genes, ",", n+1))

  if (output_assembly=="hg38"){
   input_type <- "hg38"
  }else{
    input_type <- "hg19"
  }
    gene_loc_df <- get_genes_loc(loop_genes, assembly = input_type)

    gene_loc_input <- gene_loc_df[, c("gene_symbol", "start_position", "end_position")]
    gene_loc_input <- data.frame(chrom=rep(chrom, nrow(gene_loc_input)), start_pos=gene_loc_input$start_position, end_pos=gene_loc_input$end_position, gene=gene_loc_input$gene_symbol)

  loop_snps <- snp_info_table[snp_info_table$Info=="3D-Interacting SNPs",2]
  n <- str_count(loop_snps, ",")
  loop_snps <- as.vector(str_split_fixed(loop_snps, ",", n+1))

  loop_snps_filtered <- c()
  snps_loc <- c()
  for (i in 1:length(loop_snps)){
    s <-loop_snps[i]
    s_loc <- get_snp_loc(s, assembly = input_type)
    if (s_loc != "NA:NA"){
      snps_loc <- c(snps_loc, s_loc)
      loop_snps_filtered <- c(loop_snps_filtered, loop_snps[i])
    }
  }
  snps_loc_input <- data.frame(chrom=str_split_fixed(snps_loc, ":", 2)[,1], start_pos=as.numeric(str_split_fixed(snps_loc, ":", 2)[,2]), end_pos=as.numeric(str_split_fixed(snps_loc, ":", 2)[,2]), rsID=loop_snps_filtered)

  curr_loc_input <- data.frame(chrom=chrom, start_pos=pos, end_pos=pos, rsID=snp)
  snps_loc_input <- rbind(snps_loc_input, curr_loc_input)


  start_pos <- max(0, as.numeric(pos)-window_size)
  end_pos <- min(as.numeric(pos)+window_size, chrom_max_len)
  chrom_input <- data.frame(chrom=chrom, start_pos=start_pos, end_pos=end_pos)

  # cCREs
  if (output_assembly=="hg19"){
    cCRE_data <- cCRE_data_hg19
  }else{
    cCRE_data <-cCRE_data_hg38
  }

  colnames(cCRE_data) <- c("chrom", "start", "end", "accession", "cCRE_accession", "type")
  cCRE_data_intersect <- cCRE_data[(cCRE_data$chrom==chrom) & (cCRE_data$start >= start_pos) & (cCRE_data$end <= end_pos), ]

  cCRE_data_group <- cCRE_data_intersect
  cCRE_data_group$type <- gsub(pattern = ",CTCF-bound", "", cCRE_data_intersect$type)
  cCRE_data_input <- cCRE_data_group[c("chrom", "start", "end", "type")]
  cCRE_all_colors <- data.frame(type=c("dELS", "pELS","PLS", "DNase-H3K4me3", "CTCF-only"), cols=c("#FFCD00", "#FFA700", "#FF0000", "#ffaaaa", "#00B0F0"))

  cCRE_colors <- c()
  for (i in 1:nrow(cCRE_data_input)){
    cCRE_colors <- c(cCRE_colors, cCRE_all_colors[cCRE_all_colors$type==cCRE_data_input$type[i],2])
  }

  # Draw

  pdf(file = savefile, width = 10, height = 7)

  circle_size = unit(1, "snpc")
  circos.par("start.degree" = 90)
  if (circos_pos=="Absolute"){
    circos.genomicInitialize(chrom_input, plotType = 'axis', axis.labels.cex = 0.45, tickLabelsStartFromZero = FALSE)}else{
    circos.genomicInitialize(chrom_input, plotType = 'axis', axis.labels.cex = 0.7, tickLabelsStartFromZero = TRUE)
    }

  circos.genomicTrackPlotRegion(
    cCRE_data_input, track.height = 0.1, stack = TRUE, bg.border = "gray",
        panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = cCRE_colors,border = NA, ...)
          }
    )

  # Gene
  gene_cols <- pal_igv(alpha = 0.6)(nrow(gene_loc_input))

  circos.genomicTrackPlotRegion(
    gene_loc_input, track.height = 0.1, stack = TRUE, bg.border = "gray",
        panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value,col = gene_cols,border = NA, ...)
          }
    )

  # 3D-loop SNPs
  circos.genomicTrackPlotRegion(
    snps_loc_input, track.height = 0.1, stack = TRUE, bg.border = NA,
        panel.fun = function(region, value, ...) {
        circos.genomicPoints(region, value, pch = 18, cex = 1.2, col = c(rep("#00C322", nrow(snps_loc_input)-1), "#FF1300"),border = NA, ...)
        },
    )

  if (nrow(gene_loc_input!=0)){
    curr_snp <- data.frame(chrom=chrom, start_pos=as.numeric(pos), end_pos=as.numeric(pos), rsID=snp, freq=1:nrow(gene_loc_input))
    rcols <- scales::alpha(gene_cols,alpha=0.4)
    circos.genomicLink(curr_snp, gene_loc_input, directional=1,col=rcols)
  }

  circos.genomicLabels(snps_loc_input, labels.column = 4, side = "inside", connection_height = 0.2, padding = 3,
  col = c(rep("#009393", nrow(snps_loc_input)-1), "#FF1300"), cex = 0.8, labels_height = 0.25, niceFacing = TRUE)

  legend_labels <- unique(cCRE_data_input$type)
  at_colors <- c()
  for (i in legend_labels){
    at_colors <- c(at_colors, cCRE_all_colors[cCRE_all_colors$type==i,2])
  }

  cCRE_legend <- Legend(labels = legend_labels, legend_gp = gpar(fill=at_colors, fontsize = 6),
  grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'), title="cCREs")
  pushViewport(viewport(x = 0.9, y = 0.7))
  grid.draw(cCRE_legend)
  upViewport()

  gene_legend <- Legend(labels = gene_loc_input$gene, legend_gp = gpar(fill=gene_cols, fontsize = 6),
  grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'), title="3D-interacting genes")
  pushViewport(viewport(x = 0.9, y = 0.4))
  grid.draw(gene_legend)
  upViewport()

  title(chrom)
  dev.off()
  circos.clear()

}
