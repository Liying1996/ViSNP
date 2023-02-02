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
#' pdf("affected_genes.pdf", height = 7, width = 8)
#' plot_affect_gene(data, plot_type = "merged")
#' dev.off()




plot_affect_gene <- function(data, plot_type="merged", show_num=7){
    if (!plot_type %in% c("all", "gene", "snp", "feature", "merged")){
        return(message("Invalid parameter! Please use one parameter in 'all','gene', 'snp', 'feature', 'merged'."))
    }

    new_data <- data[, c("Uploaded_variation", "Gene", "Feature_type")]
    new_data <- unique(new_data)
    genes <- data.frame(table(new_data$Gene))
    colnames(genes) <- c("Gene", "Freq")
    unaffect_num <- genes[genes$Gene == "-", 2]
    affect_num <- sum(genes[genes$Gene != "-", 2])

    genes_df <- data.frame(c("Affect", "Unaffect"), c(affect_num, unaffect_num))
    colnames(genes_df) <- c("type", "freq")
    g1 <- ggplot(data = genes_df, aes(x = type, y = freq, fill = type)) + geom_bar(stat = "identity") +
        scale_fill_jco() +
        guides(fill = guide_legend(title=NULL)) +
        labs(x = "Affected Genes Frequency", y = "Frequency") +
        theme_snp() +
        theme(plot.title = element_text(hjust = 0.5, size = 12)) +
        labs(title = "Summay of affected genes")


    genes <- genes[genes$Gene != "-", ]
    genes <- genes[order(genes$Freq, decreasing = T), ]
    genes$Gene <- factor(genes$Gene, levels=rev(genes$Gene))
    genes <- genes[1:show_num,]
    g2 <- ggplot(genes, aes(x = Gene, y = Freq, fill = Gene)) + geom_bar(stat = "identity") +
        coord_flip() +
        guides(fill = "none") +
        labs(y = "Frequency") +
        theme_snp()  +
        labs(title = "Top affected genes") +
        theme(axis.text.x = element_text(vjust=0.5), axis.text = element_text(size=10)) +
        theme(plot.title = element_text(hjust = 0.5, size = 12)) +
        scale_fill_nejm()

    # The Number of affected genes for each SNP
    # new_data <- data[, c("Uploaded_variation", "Gene", "Feature_type", "Biotype")]
    # new_data <- unique(new_data)
    genes <- data.frame(table(new_data$Uploaded_variation))
    colnames(genes) <- c("Variation", "Freq")
    genes <- genes[order(genes$Variation, decreasing = T), ]
    genes <- genes[1:show_num,]
    genes <- genes[order(genes$Freq, decreasing = F), ]
    genes$Variation <- factor(genes$Variation, levels=genes$Variation)

    g3 <- ggplot(genes, aes(x = Variation, y = Freq, fill = Variation, 20)) + geom_bar(stat = "identity") +
        coord_flip() +
        guides(fill = "none") +
        labs(y = "Frequency", title = "Affected genes of each SNP") +
        theme_snp() +
        theme(axis.text.x = element_text(vjust=0.5), axis.text = element_text(size=10)) +
        theme(plot.title = element_text(hjust = 0.5, size = 12))  +
        scale_fill_lancet()


    # Affected feature Type
    feature <- data.frame(table(new_data$Feature_type))
    colnames(feature) <- c("Feature", "Freq")
    feature <- feature[feature$Feature != "-",]
    feature <- feature[order(feature$Feature, decreasing = T), ]
    g4 <- ggplot(feature, aes(x = Feature, y = Freq, fill = Feature, 20)) + geom_bar(stat = "identity") +
        guides(fill = "none") +
        labs(y = "Frequency", title = "Affected feature types") +
        theme_snp() +
        theme(axis.text.x = element_text(vjust=0.5), axis.text = element_text(size=10)) +
        theme(plot.title = element_text(hjust = 0.5, size = 12)) +
        scale_fill_igv()

    if (plot_type=="all"){
        print(g1)
        print(g2)
        print(g3)
        print(g4)
    }else if (plot_type=="gene"){
        print(g1)
        print(g2)
    }else if (plot_type=="snp"){
        print(g2)
    }else if (plot_type=="feature"){
        print(g3)
    }else{
        grid.arrange(g1, g2, g3, g4, ncol=2)
    }
}
