#' Return boxplots of allele frequency.
#'
#' @param data Required. The annotation results from VEP.
#' @param version Optional."1000G" or "gnomAD" can be selected. Default is "1000G".
#'
#' @return A boxplot.
#' @export
#'
#' @examples
#' plot_batch_AF(data, version = "1000G")


plot_batch_AF <- function(data, version="1000G"){
    if (!version %in% c("1000G", "gnomAD")){
        return(message("Invalid parameter! Please use one parameter in '1000G', 'GnomAD'."))
        }

    g1000 <- c("AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "AA_AF", "EA_AF")
    gnomAD <- c("gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF","gnomAD_OTH_AF", "gnomAD_SAS_AF")

    af_version <- g1000
    if (version=="gnomAD"){
        af_version <- gnomAD
    }

    new_data <- data[, c("Uploaded_variation", af_version)]
    new_data <- unique(new_data)

    if (version=="gnomAD"){
        annotate_ratio <- sum(new_data$gnomAD_AF != "-") / nrow(new_data)
    }else{
        annotate_ratio <- sum(new_data$AF != "-") / nrow(new_data)
    }

    if (annotate_ratio < 0.2){
        print("Warning: the percentage of annotated SNPs of this version is less than 20% !")
    }

    if (version=="1000G"){
        melt_data <- melt(new_data, id.vars = c("Uploaded_variation"), measure.vars = af_version, variable.name = "group", value.name = "AF")
        melt_data$AF <- as.numeric(melt_data$AF)
        g <- ggplot(melt_data, aes(x = group, y = AF, fill = group)) + geom_boxplot(color = "gray") +
            scale_fill_jco() +
            labs(x = "Population", y  = "Allele Frequency", title = "Allele Frequency in different population") +
            theme_snp() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            guides(fill = guide_legend(NULL))
    }else{
        melt_data <- melt(new_data, id.vars = c("Uploaded_variation"), measure.vars = af_version, variable.name = "group", value.name = "gnomAD_AF")
        melt_data$gnomAD_AF <- as.numeric(melt_data$gnomAD_AF)
        g <- ggplot(melt_data, aes(x = group, y = gnomAD_AF, fill = group)) + geom_boxplot(color = "gray") +
            scale_fill_jco() +
            labs(x = "Population", y  = "Allele Frequency", title = "Allele Frequency in different population") +
            theme_snp() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
            guides(fill = guide_legend(NULL))
    }
    print(g)
}


