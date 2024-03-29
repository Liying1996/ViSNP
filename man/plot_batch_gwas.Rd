% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_batch_gwas.R
\name{plot_batch_gwas}
\alias{plot_batch_gwas}
\title{Return barplots of associated GWAS phenotypes of SNPs.}
\usage{
plot_batch_gwas(snps_loc, assembly = "hg38", show_num = 5, enrichment = FALSE)
}
\arguments{
\item{snps_loc}{Required. snps_loc Required. The locations of SNPs, the format should be "chr1:1014863".}

\item{assembly}{Optional. The assembly version of the input SNPs. "hg19" and "hg38" can be selected. Default is "hg38". Default is "hg38".}

\item{show_num}{Optional. The number of associated phenotypes shown on the plot. Default is 5.}

\item{enrichment}{Optional. Whether to do the enrichment analysis of GWAS studies. Default is FALSE. (Users can do the enrichment analysis by analyze_gwas_enrich() as well)}
}
\value{
A barplot.
}
\description{
Return barplots of associated GWAS phenotypes of SNPs.
}
\examples{
plot_batch_gwas(data)
}
