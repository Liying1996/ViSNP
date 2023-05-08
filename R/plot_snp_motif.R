#' Return the disrupted motif sequences.
#'
#' @param snp Required.
#' @param input_type Optional. The assembly version of the input SNP. "rsID", "hg19" and "hg38" can be selected. Default is "rsID".
#'
#' @return A plot.
#' @export
#'
#' @examples
#' plot_snp_motif(snp="rs10040658")


plot_snp_motif <- function(snp, input_type="rsID", tf="all"){

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

  motifs <- get_snp_motif(snp_id)
  if (tf[1]!="all"){
    if (sum(tf %in% motifs$Motif) != length(tf)){
      return(message("Invalid TF name!"))
    }else{
      motifs <- motifs[motifs$Motif%in%tf, ]
    }
  }

  seq_x = c()
  seq_y = c()
  seq_label = c()
  seq_col <- c()
  seq_n_total <- c()

  axis_x <- c()
  axis_y <- c()
  axis_label <- c()


  for (i in 1:nrow(motifs)){
    seqs <- strsplit(motifs$Sequence[i], "")[[1]]
    sour <- motifs$Source[i]
    motif_name <- motifs$Motif[i]
    alter_pos <- as.numeric(motifs$Altered_Pos[i])

    seq_n <- length(seqs)
    seq_n_total <- c(seq_n_total, seq_n)
    seq_x <- c(seq_x, seq(0, seq_n - 1))
    seq_y <- c(seq_y, rep(i, seq_n))
    seq_label <- c(seq_label, seqs)
    c <- rep('blue', seq_n)
    c[alter_pos+1] <- 'red'
    seq_col <- c(seq_col, c)

    axis_x <- c(axis_x,  seq(0, seq_n - 1))
    axis_y <- c(axis_y, rep(i-0.25, seq_n))
    axis_label <- c(axis_label, seq(0, seq_n - 1))

    if (i == 1){
      seq_axis = tibble(
        x = seq(0, seq_n-1),
        y = rep(i - 0.18, seq_n)
      )
    }else{
      seq_axis = add_row(seq_axis, tibble(
        x = seq(0, seq_n-1),
        y = rep(1 * i - 0.18, seq_n)
      ))
    }
  }

  segments <-  tibble(
    x1 = rep(0, nrow(motifs)),
    x2 = seq_n_total-1,
    y1 = seq(1, nrow(motifs)) - 0.2,
    y2 = seq(1, nrow(motifs)) - 0.2
  )

  motif_names <- motifs$Motif
  motif_strand <- motifs$Strand
  motif_x <- rep(max(segments$x2)*0.15, nrow(motifs))
  motif_y <- seq(1, nrow(motifs)) + 0.35


  g <- ggplot() +
    geom_segment(data=segments, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.7) +
    annotate("text", x=seq_x, y=seq_y , label=seq_label, size=6, color=seq_col) +
    annotate("text", x=motif_x, y=motif_y , label=paste(motif_names, "(", motif_strand, ")", sep=" "), size=5) +
    annotate("text", x=axis_x, y=axis_y, label=axis_label, size=3.5) +
    geom_point(data = seq_axis, aes(x=x, y=y), shape="I", size=3) +
    theme(panel.background = element_rect(fill="white", colour = "white")) +
    labs(title = snp_id) +
    theme(
      axis.title=element_blank(),
      axis.text=element_blank(),
      axis.ticks=element_blank()) +
    theme(plot.title = element_text(hjust=0.5, size=18))

  return(g)
}
