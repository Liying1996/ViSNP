#' Return barplots of affected genes of SNPs.
#'
#' @param data Required. The annotation results from VEP.
#' @param data_source Optional. The source of the inpit data. "API" or "Upload" can be selected. Default is "API" (Input annotated data from get_batch_vep() function).
#' @param gene Optional. The version of gene names. "Ensembl" and "Symbol" can be selected.  Default is "Ensembl".
#' @param plot_type Optional."gene", "snp", "all" and "merged" can be selected. Default is "rsID".
#' @param show_num Optional. The number of affected genes shown on the plot. Default is 7.
#' @param go_enrichment Optional. Whether to do the GO enrichment analysis of genes. Default is FALSE. (Users can do the enrichment analysis by analyze_gwas_enrich() as well)
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' data <- test_upload
#' pdf("affected_genes.pdf", height = 3, width = 8)
#' plot_affect_gene(data, data_source="Upload", gene = "Symbol", plot_type = "merged")
#' dev.off()


plot_affect_gene <- function(data, data_source="API", gene = "Ensembl", plot_type="merged", show_num=7, go_enrichment=FALSE){
    if (!plot_type %in% c("all", "gene", "snp", "merged")){
        return(message("Invalid parameter! Please use one parameter in 'all','gene', 'snp', 'feature', 'merged'."))
    }

    if (data_source=="Upload"){
      new_data <- data[, c("Existing_variation", "Gene")]
      new_data <- unique(new_data)
      colnames(new_data) <- c('SNP', 'Gene')
      unaffect_num <- sum(new_data$Gene == '-')
      affect_num <- sum(new_data$Gene != '-')

    }else{
      genes <- c()
      snps <- c()
      for (i in 1:nrow(data)){
        curr <- unique(data$transcript_consequences[[i]])
        if ('gene_id' %in% colnames(curr)){
          curr_genes <- unique(curr$gene_id)
        }else{curr_genes <- '-'}
        genes <- c(genes, curr_genes)
        snps <- c(snps, rep(data$id[i], length(curr_genes)))
      }
      new_data <- data.frame(SNP=snps, Gene=genes)
      new_data <- unique(new_data)
      unaffect_num <- sum(new_data$Gene == '-')
      affect_num <- sum(new_data$Gene != '-')
    }

    genes_df <- data.frame(c("Affect", "Unaffect"), c(affect_num, unaffect_num))
    colnames(genes_df) <- c("type", "freq")
    g1 <- ggplot(data = genes_df, aes(x = type, y = freq, fill = type)) + geom_bar(stat = "identity") +
        scale_fill_jco() +
        guides(fill = guide_legend(title=NULL)) +
        labs(x = "", y = "Frequency") +
        theme_snp() +
        theme(axis.title = element_text(size = 11), axis.text = element_text(size = 9)) +
        theme(plot.title = element_text(hjust = 0.5, size = 11)) +
        labs(title = "Summay of affected genes")

    new_genes <- new_data[new_data$Gene != "-", ]
    new_genes <- data.frame(table(new_genes$Gene))
    colnames(new_genes) <- c('Gene', 'Freq')

    if (gene == "Symbol"){
      gene_symbols <- bitr(new_genes$Gene, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
      # new_genes <- new_genes[new_genes$Gene %in% gene_symbols$ENSEMBL, ]
      tmp_data <- merge(new_genes, gene_symbols, by.x="Gene", by.y="ENSEMBL")
      new_genes <- tmp_data[,c('SYMBOL', 'Freq')]
      colnames(new_genes) <- c('Gene', 'Freq')

    }

    new_genes <- new_genes[order(new_genes$Freq, decreasing = T), ]
    new_genes <- unique(new_genes)
    new_genes <- new_genes[1:show_num,]
    new_genes$Gene <- factor(new_genes$Gene, levels=rev(new_genes$Gene))

    g2 <- ggplot(new_genes, aes(x = Gene, y = Freq, fill = Gene)) + geom_bar(stat = "identity") +
        coord_flip() +
        guides(fill = "none") +
        labs(y = "Frequency") +
        theme_snp()  +
        labs(title = "Top affected genes") +
        theme(axis.text.x = element_text(vjust=0.5), axis.text = element_text(size=10)) +
        theme(plot.title = element_text(hjust = 0.5, size = 12)) +
        scale_fill_nejm()

    # The Number of affected genes for each SNP
    new_SNPs <- new_data[new_data$Gene != "-", ]
    new_SNPs <- data.frame(table(new_SNPs$SNP))
    colnames(new_SNPs) <- c("SNP", "Freq")
    new_SNPs$SNP <- str_split_fixed(new_SNPs$SNP, ',', 2)[,1]
    new_SNPs <- new_SNPs[order(new_SNPs$SNP, decreasing = T), ]
    new_SNPs <- new_SNPs[1:show_num,]
    new_SNPs <- new_SNPs[order(new_SNPs$Freq, decreasing = F), ]
    new_SNPs$SNP <- factor(new_SNPs$SNP, levels=new_SNPs$SNP)

    g3 <- ggplot(new_SNPs, aes(x = SNP, y = Freq, fill = SNP, 20)) + geom_bar(stat = "identity") +
        coord_flip() +
        guides(fill = "none") +
        labs(y = "Frequency", title = "Affected genes of each SNP") +
        theme_snp() +
        theme(axis.text.x = element_text(vjust=0.5), axis.text = element_text(size=10)) +
        theme(plot.title = element_text(hjust = 0.5, size = 12))  +
        scale_fill_lancet()

    # Affected feature Type
    # feature <- data.frame(table(new_data$Feature_type))
    # colnames(feature) <- c("Feature", "Freq")
    # feature <- feature[feature$Feature != "-",]
    # feature <- feature[order(feature$Feature, decreasing = T), ]
    # g4 <- ggplot(feature, aes(x = Feature, y = Freq, fill = Feature, 20)) + geom_bar(stat = "identity") +
    #     guides(fill = "none") +
    #     labs(y = "Frequency", title = "Affected feature types") +
    #     theme_snp() +
    #     theme(axis.text.x = element_text(vjust=0.5), axis.text = element_text(size=10)) +
    #     theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    #     scale_fill_igv()

    if (plot_type=="all"){
        print(g1)
        print(g2)
        print(g3)
        # print(g4)
    }else if (plot_type=="gene"){
        print(g1)
        print(g2)
    }else if (plot_type=="snp"){
        print(g2)
    }else if (plot_type=="feature"){
        print(g3)
    }else{
        grid.arrange(g1, g2, g3, nrow=1)
    }

    if (go_enrichment){
      snp_genes <- unique(new_data[new_data$Gene != "-", 2])
      go_plot <- analyze_gene_go(snp_genes, output_type="dotplot")
      print(go_plot)
    }
}
