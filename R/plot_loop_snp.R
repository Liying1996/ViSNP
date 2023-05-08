#' Plot 3D-interating SNPs.
#'
#' @param snp Required.
#' @param input_type Optional. The assembly version of the input SNP. "rsID", "hg19" and "hg38" can be selected. Default is "rsID".
#' @param pop Required. Population.
#' @param output_assembly Optional. The UCSC Assembly versions of chromosome coordinates of outputs. "hg19" or "hg38" can be selected. Default is "hg19".
#' @param show_cells Optional. The number of cell types shown in the figure. Default is 3 (Top 3 cell types). The names of cell types are available as well.
#'
#' @return A plot.
#' @export
#'
#' @examples
#' plot_loop_snp(snp="rs10", pop="AFR")
#' plot_loop_snp(snp="rs10", pop="AFR", output_assembly='hg38')
#' plot_loop_snp(snp="rs10", pop = pop, show_cells = c("Ventricle_Right", "Spleen"))
#' plot_loop_snp(snp='rs10040658', pop = "AFR", show_cells = c("GM12878", "Lung"))


plot_loop_snp <- function(snp, input_type="rsID", pop, output_assembly="hg19", show_cells=3){

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

  loop_snps <- get_loop_snp(snp_id)
  loop_snps <- loop_snps[complete.cases(loop_snps),]

  loop_snps <- loop_snps[loop_snps$type == pop, ]

  cell_types <- unique(loop_snps$cell)
  if (class(show_cells)=="numeric"){
    if (length(cell_types) > show_cells){
      cell_types <- cell_types[1:show_cells]
    }
  }

  if (class(show_cells)=="character"){
    if (sum(show_cells %in% cell_types)==length(show_cells)){
      loop_snps <- loop_snps[loop_snps$cell %in% show_cells, ]
    }else{
      return(message("Invalid cell type!"))
    }
  }

  cell_types <- unique(loop_snps$cell)
  chrom <- unique(loop_snps$chr)
  if (output_assembly=="hg19"){
    snp_loc <- get_snp_loc(snp_id, "hg19")
  }else{
    snp_loc <- get_snp_loc(snp_id, "hg38")
  }
  snp_pos <- as.numeric(str_split_fixed(snp_loc, ":", 2)[2])

  colors_default <- pal_igv(alpha = 0.8)(9)
  c_index <- 2

  scales_total_label <- c()
  scales_total_x <- c()
  scales_total_y <- c()
  cell_y <- c()
  cell_labels <- c()
  snp_text_x <- c()
  snp_text_y <- c()
  snp_text_label <- c()

  for (i in 1:length(cell_types)){
    seg_color <- colors_default[1]
    snp_B_id <- loop_snps[loop_snps$cell == cell_types[i], ]$SNP_B

    start_cha <- str_split_fixed(loop_snps[loop_snps$cell == cell_types[i], ]$start, ":" , 2)[,2]

    if (length(start_cha)==1){
      start1 <- as.numeric(str_split_fixed(start_cha, "-", 2)[1])
      start2 <- as.numeric(str_split_fixed(start_cha, "-", 2)[2])
      end_cha <- str_split_fixed(loop_snps[loop_snps$cell == cell_types[i], ]$end, ":" , 2)[,2]
      end1 <- as.numeric(str_split_fixed(end_cha, "-", 2)[1])
      end2 <- as.numeric(str_split_fixed(end_cha, "-", 2)[2])
    }else{
      start1 <- as.numeric(str_split_fixed(start_cha, "-", 2)[,1])
      start2 <- as.numeric(str_split_fixed(start_cha, "-", 2)[,2])
      end_cha <- str_split_fixed(loop_snps[loop_snps$cell == cell_types[i], ]$end, ":" , 2)[,2]
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


    snp_B_pos <- c()
    for (j in snp_B_id){
      if (output_assembly=="hg19"){
        snp_B_loc <- get_snp_loc(j, 'hg19')
        snp_B_pos <-  c(snp_B_pos, as.numeric(str_split_fixed(snp_B_loc, ":", 2)[2]))
      }else{
        snp_B_loc <- get_snp_loc(j, 'hg38')
        snp_B_pos <-  c(snp_B_pos, as.numeric(str_split_fixed(snp_B_loc, ":", 2)[2]))
      }
    }

    anchor1_start <- (start1 - min(start1))/(max(end2) - min(start1)) * 10
    anchor1_end <- (start2 - min(start1))/(max(end2) - min(start1)) * 10
    anchor2_start <- (end1 - min(start1))/(max(end2) - min(start1)) * 10
    anchor2_end <- (end2 - min(start1))/(max(end2) - min(start1)) * 10
    snp_input <- (snp_pos - min(start1))/(max(end2) - min(start1)) * 10

    snp_B <- (snp_B_pos - min(start1))/(max(end2) - min(start1)) * 10


    if (i == 1){
      snp_input_total <- tibble(
        x = snp_input,
        y = 2 * i
      )

      points_B <- tibble(
        x = snp_B,
        y = rep(2 * i, length(snp_B))
      )

      arrows <- tibble(
        x1 = rep(snp_input, length(snp_B)),
        x2 = snp_B,
        y1 = rep(2 * i + 0.03, length(snp_B)),
        y2 = rep(2 * i + 0.03, length(snp_B))
      )

      segments <- tibble(
        x1 = c(anchor1_start, anchor2_start),
        x2 = c(anchor1_end, anchor2_end),
        y1 = rep(2 * i, length(c(anchor1_start, anchor2_start))),
        y2 = rep(2 * i, length(c(anchor1_start, anchor2_start))),
        color = rep(seg_color, 2)
      )
    }else{
      snp_input_total <- add_row(snp_input_total, tibble(
        x = snp_input,
        y = 2 * i
      ))
      points_B <- add_row(points_B, tibble(
        x = snp_B,
        y = rep(2 * i, length(snp_B))
      ))

      arrows <- add_row(arrows, tibble(
        x1 = rep(snp_input, length(snp_B)),
        x2 = snp_B,
        y1 = rep(2 * i + 0.03, length(snp_B)),
        y2 = rep(2 * i + 0.03, length(snp_B))
      ))

      segments <- add_row(segments, tibble(
        x1 = c(anchor1_start, anchor2_start),
        x2 = c(anchor1_end, anchor2_end),
        y1 = rep(2 * i, length(c(anchor1_start, anchor2_start))),
        y2 = rep(2 * i, length(c(anchor1_start, anchor2_start))),
        color = rep(seg_color, 2)
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

    snp_B2 <- c()
    snp_B_id2 <- c()
    for (x in 1:length(snp_B)){
      if (x > 1){
        if ((snp_B[x] - snp_B[x-1] < 0.2) & (snp_B_id[x] != snp_B_id[x-1])){
            snp_B2 <- c(snp_B2, snp_B[x] + 0.15)
            snp_B_id2 <- c(snp_B_id2, snp_B_id[x])
          }else if (snp_B_id[x] != snp_B_id[x-1]){
            snp_B2 <- c(snp_B2, snp_B[x])
            snp_B_id2 <- c(snp_B_id2, snp_B_id[x])
          }
      }else{
        snp_B2 <- snp_B[1]
        snp_B_id2 <- c(snp_B_id2, snp_B_id[1])
        }
    }

    snp_text_x = c(snp_text_x, snp_input, snp_B2)
    snp_text_y = c(snp_text_y, rep(i * 2 - 0.3, length(snp_B2)+1))
    snp_text_label <- c(snp_text_label, snp_id, snp_B_id2)
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

  segments <- segments[!duplicated(segments),]
  tmp_df <- data.frame(snp_text_x=snp_text_x, snp_text_y=snp_text_y, snp_text_label=snp_text_label)
  tmp_df <- tmp_df[!duplicated(tmp_df),]

  g <- ggplot() + geom_segment(data = sequence, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.7) +
    geom_segment(data=segments, aes(x=x1, y=y1, xend=x2, yend=y2, color=color), alpha=0.7, size=2) +
    geom_point(data = snp_input_total, aes(x=x, y=y), shape=18, color="#F75000", size=3) +
    geom_point(data = points_B, aes(x=x, y=y), shape=18, color="#00A600", size=3) +
    geom_point(data = seq_label, aes(x=x, y=y), shape="I", size=1.5) +
    geom_curve(data = arrows, aes(x=x1, y=y1, xend=x2, yend=y2), arrow = arrow(length = unit(0.01, "npc")), size = 0.5, color="purple", curvature=0.8) +
    annotate("text", x=tmp_df$snp_text_x, y=tmp_df$snp_text_y,label=tmp_df$snp_text_label, size=2, color="black", angle=45, alpha=0.6) +
    annotate("text", x=rep(0.6, length(cell_types)), y=0.7+2*c(1:length(cell_types)), label=cell_types, size=3.5) +
    annotate("text", x=10, y=1.3, label="Mb", size=3) +
    annotate("text",x=scales_total_x, y=scales_total_y, label=scales_total_label, size=2.5) +
    guides(color="none") +
    labs(x = chrom, title = snp_id) +
    theme(panel.background = element_rect(fill="white", colour = "white")) +
    theme(
      axis.title.y=element_blank(),
      axis.text=element_blank(),
      axis.ticks=element_blank()) +
    theme(plot.title = element_text(hjust=0.5))

  return(g)

}
