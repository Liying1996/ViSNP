#' Plot 3D-loop genes.
#'
#' @param snp Required.
#' @param input_type Optional. The assembly version of the input SNP. "rsID", "hg19" and "hg38" can be selected. Default is "rsID".
#' @param output_assembly Optional. The UCSC Assembly versions of chromosome coordinates of outputs. "hg19" or "hg38" can be selected. Default is "hg19".
#' @param show_cells Optional. The number of cell types shown in the figure. Default is 3 (Top 3 cell types). The names of cell types are available as well.
#'
#' @return A plot.
#' @export
#'
#' @examples
#' plot_loop_gene(snp="rs10")
#' plot_loop_gene(snp="rs10", output_assembly='hg38', show_cells=c("VentricleLeft", "Caki2", "HepG2"))
#' plot_loop_gene(snp="rs10040658", show_cells=c("GM12878", "A549"))

plot_loop_gene <- function(snp, input_type="rsID", output_assembly="hg19", show_cells=3){

  if  (!input_type %in% c("rsID", "hg19", "hg38")){
    return(message("Please select 1 input type from rsID, hg19 or hg38!"))
  }

  if (input_type=="hg38"){
    snp_id <- get_snp_rsID(snp, assembly = "hg38")
  }else if(input_type=="hg19"){
    snp_id <- get_snp_rsID(snp, assembly = "hg19")
  }else{
    snp_id <- snp
  }

  loop_genes <- get_loop_gene(snp_id)
  loop_genes <- loop_genes[loop_genes$gene!="",]
  loop_genes <- loop_genes[!duplicated(loop_genes),]

  cell_types <- unique(loop_genes$cell)

  if (class(show_cells)=="numeric"){
    if (length(cell_types) > show_cells){
      cell_types <- cell_types[1:show_cells]
      loop_genes <-  loop_genes[loop_genes$cell %in% cell_types, ]
    }
  }

  if (class(show_cells)=="character"){
    if (sum(show_cells %in% cell_types)==length(show_cells)){
      loop_genes <- loop_genes[loop_genes$cell %in% show_cells, ]
    }else{
      return(message("Invalid cell type!"))
    }
  }

  cell_types <- unique(loop_genes$cell)

  chrom <- str_split_fixed(loop_genes$start[1], ":", 2)[1]
  if (output_assembly=="hg19"){
    snp_loc <- get_snp_loc(snp_id, "hg19")
  }else{
    snp_loc <- get_snp_loc(snp_id, "hg38")
  }
  snp_pos <- as.numeric(str_split_fixed(snp_loc, ":", 2)[2])

  colors_default <- pal_igv(alpha = 0.8)(9)

  gene_types <- unique(loop_genes$gene)
  colors_genes_df <- data.frame(gene=gene_types, color=pal_jco(alpha = 0.8)(length(gene_types)))
  c_index <- 2

  scales_total_label <- c()
  scales_total_x <- c()
  scales_total_y <- c()
  cell_y <- c()
  cell_labels <- c()
  gene_text_x <- c()
  gene_text_y <- c()
  gene_text_label <- c()

  for (i in 1:length(cell_types)){
    seg_color <- colors_default[1]
    loop_genes_curr <- loop_genes[loop_genes$cell == cell_types[i], ]
    gene_id <- loop_genes_curr$gene

    if (output_assembly=="hg19"){
      gene_loc <- get_genes_loc(gene_id, 'hg19')
    }else{
      gene_loc <- get_genes_loc(gene_id, 'hg38')
    }

    gene_pos <- gene_loc[c("start_position", "end_position")]
    loop_genes_curr <- loop_genes_curr[loop_genes_curr$gene %in% gene_loc$gene_symbol,]

    start_cha <- str_split_fixed(loop_genes_curr$start, ":" , 2)[,2]

    if (length(start_cha)==1){
      start1 <- as.numeric(str_split_fixed(start_cha, "-", 2)[1])
      start2 <- as.numeric(str_split_fixed(start_cha, "-", 2)[2])
      end_cha <- str_split_fixed(loop_genes_curr$end, ":" , 2)[,2]
      end1 <- as.numeric(str_split_fixed(end_cha, "-", 2)[1])
      end2 <- as.numeric(str_split_fixed(end_cha, "-", 2)[2])
    }else{
      start1 <- as.numeric(str_split_fixed(start_cha, "-", 2)[,1])
      start2 <- as.numeric(str_split_fixed(start_cha, "-", 2)[,2])
      end_cha <- str_split_fixed(loop_genes_curr$end, ":" , 2)[,2]
      end1 <- as.numeric(str_split_fixed(end_cha, "-", 2)[,1])
      end2 <- as.numeric(str_split_fixed(end_cha, "-", 2)[,2])
      for (x in 2:length(start_cha)){
        if ((start_cha[x] == start_cha[x-1]) & (end_cha[x] == end_cha[x-1])){
          seg_color <- c(seg_color, seg_color[x-1])
        }else{
          seg_color <- c(seg_color, colors_default[c_index])
          c_index <- c_index + 1
        }
      }
    }

    if ((output_assembly=="hg38")){
      for (p in 1:length(start1)){
        start1[p] <- as.numeric(str_split_fixed(hg19tohg38(paste(chrom, start1[p], sep=':')), ":", 2)[2])
        start2[p] <- as.numeric(str_split_fixed(hg19tohg38(paste(chrom, start2[p], sep=':')), ":", 2)[2])
        end1[p] <- as.numeric(str_split_fixed(hg19tohg38(paste(chrom, end1[p], sep=':')), ":", 2)[2])
        end2[p] <-as.numeric(str_split_fixed(hg19tohg38(paste(chrom, end2[p], sep=':')), ":", 2)[2])
      }
    }

    max_end <- max(gene_loc$end_position, end2)
    min_start <- min(gene_loc$start_position, start1)

    anchor1_start <- (start1 - min_start)/(max_end - min_start) * 10
    anchor1_end <- (start2 - min_start)/(max_end - min_start) * 10
    anchor2_start <- (end1 - min_start)/(max_end - min_start) * 10
    anchor2_end <- (end2 - min_start)/(max_end - min_start) * 10
    snp_input <- (snp_pos - min_start)/(max_end - min_start) * 10


    gene_start <- (gene_pos$start_position - min_start)/(max_end - min_start) * 10
    gene_end <- (gene_pos$end_position - min_start)/(max_end - min_start) * 10
    gene_mid <- (gene_start + gene_end)/2
    gene_mid3 <- gene_mid

    # adjust positions
    gene_mid2 <- c()
    for (x in gene_mid){
      if (abs(x - snp_input) < 0.2){
        if (x > snp_input){
          x = x + 0.35
          snp_input <- snp_input - 0.2
        }else{
          x = x - 0.35
          snp_input <- snp_input + 0.2
        }
      }
      gene_mid2 <- c(gene_mid2, x)

    }
      gene_mid <- gene_mid2

      gene_colors <- c()
      for (gene in gene_loc$gene_symbol){
        gene_colors <- c(gene_colors, colors_genes_df[colors_genes_df$gene == gene, 2])
      }

      if (i == 1){
        snp_input_total <- tibble(
          x = snp_input,
          y = 2 * i
        )

      arrows <- tibble(
        x1 = rep(snp_input, length(gene_mid3)),
        x2 = gene_mid3,
        y1 = rep(2 * i + 0.01, length(gene_mid3)),
        y2 = rep(2 * i + 0.01, length(gene_mid3))
      )

      segments <- tibble(
        x1 = c(anchor1_start, anchor2_start),
        x2 = c(anchor1_end, anchor2_end),
        y1 = rep(2 * i, length(c(anchor1_start, anchor2_start))),
        y2 = rep(2 * i, length(c(anchor1_start, anchor2_start))),
        color = rep(seg_color, 2)
      )

      segments_gene <- tibble(
        x1 = gene_start,
        x2 = gene_end,
        y1 = rep(2*i, length(gene_start)),
        y2 = rep(2*i, length(gene_end)),
        color = gene_colors
      )

    }else{
      snp_input_total <- add_row(snp_input_total, tibble(
        x = snp_input,
        y = 2 * i
      ))

      arrows <- add_row(arrows, tibble(
        x1 = rep(snp_input, length(gene_mid3)),
        x2 = gene_mid3,
        y1 = rep(2 * i + 0.01, length(gene_mid3)),
        y2 = rep(2 * i + 0.01, length(gene_mid3))
      ))

      segments <- add_row(segments, tibble(
        x1 = c(anchor1_start, anchor2_start),
        x2 = c(anchor1_end, anchor2_end),
        y1 = rep(2 * i, length(c(anchor1_start, anchor2_start))),
        y2 = rep(2 * i, length(c(anchor1_start, anchor2_start))),
        color = rep(seg_color, 2)
      ))



      segments_gene <- add_row(segments_gene, tibble(
        x1 = gene_start,
        x2 = gene_end,
        y1 = rep(2*i, length(gene_start)),
        y2 = rep(2*i, length(gene_end)),
        color = gene_colors
      ))
    }
    scales_x1 <- round(min(start1)/1000000, 2)
    scales_gap <- round((max(end2) - min(start1))/1000000/4, 2)
    scales_x <- c(scales_x1, scales_x1+scales_gap, scales_x1+2*scales_gap, scales_x1+3*scales_gap, scales_x1+4*scales_gap)
    scales_total_label <- c(scales_total_label, scales_x)
    scales_total_y <- c(scales_total_y, rep(2 * i - 0.5, 5))
    scales_total_x <- c(scales_total_x, seq(0, 10, 2.5))

    cell_y <- c(cell_y, 1 + 2 * i)
    cell_labels <- c(cell_labels, cell_types[i])

    gene_text_x = c(gene_text_x, snp_input, gene_mid)
    gene_text_y = c(gene_text_y, rep(i * 2 - 0.3, length(gene_mid)+1))
    gene_text_label <- c(gene_text_label, snp_id, gene_loc$gene_symbol)

  }


  sequence <- tibble(
    x1 = rep(0, length(cell_types)),
    x2 = rep(10, length(cell_types)),
    y1 = 2 * c(1:length(cell_types)),
    y2 = 2 * c(1:length(cell_types))
  )


  seq_label <- tibble(
    x = rep(seq(0,10,2.5), length(cell_types)),
    y = rep(2 * c(1:length(cell_types)) - 0.6, rep(5, length(cell_types)))
  )
  tmp_df <- data.frame(gene_text_x=gene_text_x, gene_text_y=gene_text_y, gene_text_label=gene_text_label)
  tmp_df <- tmp_df[!duplicated(tmp_df),]

  segments <- add_row(segments, segments_gene)
  segments <- segments[!duplicated(segments),]
  # pdf("~/Documents/Vi-SNP/Frontier/links_gene2.pdf", width=5, height = 6)
  ggplot() + geom_segment(data = sequence, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.7) +
    geom_segment(data=segments, aes(x=x1, y=y1, xend=x2, yend=y2, color=color), alpha=0.7, size=2) +
    geom_point(data = snp_input_total, aes(x=x, y=y), shape=18, color="#F75000", size=3) +
    geom_point(data = seq_label, aes(x=x, y=y), shape="I", size=1.5) +
    geom_curve(data = arrows, aes(x=x1, y=y1, xend=x2, yend=y2), arrow = arrow(length = unit(0.01, "npc")), size = 0.5, color="purple", curvature=0.7) +
    annotate("text", x=tmp_df$gene_text_x, y=tmp_df$gene_text_y , label=tmp_df$gene_text_label, size=2, color="black", angle=45, alpha=0.6) +
    annotate("text", x=rep(0.6, length(cell_types)), y=0.7+2*c(1:length(cell_types)), label=cell_types, size=3.5) +
    annotate("text", x=10, y=1.3, label="Mb", size=3) +
    annotate("text",x=scales_total_x, y=scales_total_y, label=scales_total_label, size=2.5) +
    guides(color="none") +
    scale_color_manual(values = unique(segments$color)) +
    labs(x = chrom, title = snp_id) +
    theme(panel.background = element_rect(fill="white", colour = "white")) +
    theme(
      axis.title.y=element_blank(),
      axis.text=element_blank(),
      axis.ticks=element_blank()) +
    theme(plot.title = element_text(hjust=0.5))
  # dev.off()

}
