#' Return a barplot of ti/tv ratio.
#'
#' @param data Required. The annotation results from VEP.
#'
#' @return A barplot.
#' @export
#'
#' @examples
#' plot_batch_titv(data)



plot_batch_titv <- function(data){
    allele <- str_split_fixed(unique(data$Uploaded_variation), "_", 3)[,3]
    allele_df <- data.frame(str_split_fixed(allele, "/", 2))
    colnames(allele_df) <- c("ref", "alt")
    ti <- (allele_df$ref=="C"&allele_df$alt=="T")|(allele_df$ref=="T"&allele_df$alt=="C")|(allele_df$ref=="G"&allele_df$alt=="A")|(allele_df$ref=="A"&allele_df$alt=="G")

    tv <- (allele_df$ref=="C"&allele_df$alt=="A")|(allele_df$ref=="A"&allele_df$alt=="C")|(allele_df$ref=="G"&allele_df$alt=="T")|(allele_df$ref=="T"&allele_df$alt=="G")|(allele_df$ref=="C"&allele_df$alt=="A")|(allele_df$ref=="A"&allele_df$alt=="C")|(allele_df$ref=="G"&allele_df$alt=="T")|(allele_df$ref=="T"&allele_df$alt=="G")

    titv_tatio <- signif(sum(ti)/sum(tv),2)
    titv_df <- data.frame(c("Ti", "Tv"), c(sum(ti), sum(tv)))
    colnames(titv_df) <- c("type", "freq")
    g <- ggplot(data = titv_df, aes(x = type, y = freq, fill = type)) + geom_bar(stat = "identity") +
        scale_fill_nejm() +
        guides(fill = guide_legend(title=NULL)) +
        labs(x = "", y = "Frequency") +
        theme_snp() +
        annotate("text", x = Inf, y = Inf, label = paste("Ti/Tv Ratio: ", titv_tatio), size = 5, vjust = 1.5, hjust = 1.2) +
        scale_x_discrete(labels = function(x) str_wrap(x, width=10))
    print(g)
}
