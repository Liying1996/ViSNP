#' Return barplots of consequences.
#'
#' @param data Required. The annotation results from VEP.
#' @param data_source Optional. The source of the input data. "API" or "Upload" can be selected. Default is "API" (Input annotated data from get_batch_vep() function).
#' @param show_num Optional. The number of affected genes shown on the plot. Default is 7.
#' @param colors Optional. Default is `ggsci::pal_igv(alpha = 0.8)(show_num)`.
#'
#' @return Barplots.
#' @export
#'
#' @examples
#' plot_batch_consequence(data, show_num = 10, data_source="API", olors = rainbow(10))


plot_batch_consequence <- function(data, show_num=7, data_source="API", colors="default"){
    if (colors[1] != "default"){
      if (length(colors) != show_num){
        return(message("Please use the same number of colors as show_lines!"))
      }
    }else{
      colors <- pal_igv(alpha = 0.8)(show_num)
    }

    if (data_source=="Upload"){
      new_data <- data[, c("Uploaded_variation", "Consequence")]
      new_data <- unique(new_data)
    }else{
      new_data <- data[,c('id', 'most_severe_consequence')]
      new_data <- unique(new_data)
      colnames(new_data) <- c('id', 'Consequence')
    }
    counts <- data.frame(table(new_data$Consequence))
    colnames(counts) <- c("Consequence", "Freq")
    counts <- counts[order(counts$Freq, decreasing = T), ]
    if (nrow(counts) > show_num){
      counts <- counts[1:show_num,]
    }
    counts <- counts[order(counts$Freq, decreasing = F), ]
    counts$Consequence <- gsub("_", " ", counts$Consequence)
    counts$Consequence <- factor(counts$Consequence, levels = counts$Consequence)
    g1 <- ggplot(counts, aes(x = Consequence, y = Freq, fill = Consequence)) + geom_bar(stat = "identity") +
        coord_flip() +
        guides(fill = "none") +
        labs(y = "Frequency") +
        theme_snp() +
        theme(axis.text.x = element_text(vjust=0.5), axis.text = element_text(size=12)) +
        scale_fill_manual(values=colors) +
        scale_x_discrete(labels=function(x) str_wrap(x, width=10))
    print(g1)


    # impacts <- data.frame(table(new_data$IMPACT))
    # colnames(impacts) <- c("IMPACT", "Freq")
    # impacts <- impacts[order(impacts$Freq, decreasing = T), ]
    # g2 <- ggplot(impacts, aes(x = IMPACT, y = Freq, fill = IMPACT)) +
    #     geom_bar(stat = "identity", color = "white") +
    #     coord_flip() +
    #     guides(fill = "none") +
    #     labs(y = "Frequency") +
    #     theme_snp() +
    #     theme(axis.text.x = element_text(vjust=0.5), axis.text = element_text(size=10)) +
    #     scale_fill_uchicago()
    # print(g2)


}
