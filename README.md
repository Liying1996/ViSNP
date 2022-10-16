## Vi-SNP: An R-package to query and visualize SNPs


## 1. **Introduction**

Single-nucleotide polymorphism (SNP) is an important topic for most genetic association studies. This package can query, analyze, visualize and link SNPs with functional information (e.g., genes, regulatory elements, phenotypes).

## 2. **Installation**

``` bash
devtools::install_github("Liying1996/ViSNP")
```

## 3. **Functions of the package**

Overview...

## 4. **Details of functions**

### **4.1 For single SNP**

#### ***query_snp(snp, input_type="rsID", eqtl_tissue="Whole_Blood")***

`query_snp` output the basic and functional information of a SNP.

Input_type: "rsID", "hg19" or "hg38". If "hg19" or "hg38" was selected, the ***snp*** should be "chr22:19712094". (All genomic locations should be 1-based.) eqtl_tissue: Tissue ID of the tissue of GTEx.

```{r}
info_table <- query_snp("rs1891906")
info_table <- query_snp("chr1:950243", input_type="hg19")
info_table <- query_snp("chr1:1014863 ", input_type="hg38")
```

Output:

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/1_query_snp1.png)

print:

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/1_query_snp2.png)

#### ***plot_physcores(snps, colors="default", input_type="rsID")***

Return a table and a line chart of the physcores of +/-10bp of SNPs.

`snps` rsIDs or genomic locations. colors: Colors of lines. If not specified, default is to use the colors of `rainbow()`.

`input_type`: The assembly version of the input SNP. "rsID", "hg19" and "hg38" can be selected. Default is "rsID".

```{r}
snps <- c("rs2769466", "rs4661232")
physcores_table <- plot_physcores(snps, colors = c("blue", "orange"))
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/2_plot_physcores.png)

#### ***get_snp_rsID(snp_loc, assembly="hg19")***

Convert the chromosome coordinate of a SNP to rsID. 

```{r}
get_snp_rsID("chr22:19712094") # hg19
get_snp_rsID("chr22:19724571", assembly="hg38") 
```

#### ***get_snp_loc(snp, assembly="hg19")***

Convert SNP rsID to genomic locations.

```{r}
get_snp_loc("rs1059196")
```

#### ***get_snp_ccre(snp, input_type="rsID", output_assembly="hg19")***

Return the cCREs (candidate cis-regulatory elements) that SNPs overlapped.

```{r}
snp <- "rs1059196"
cCRE_intersect <- get_snp_ccre(snp)
# get_snp_ccre("chr22:19724571", input_type = "hg38", output_assembly = "hg38")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/3_get_snp_ccre.png)

#### ***get_snp_gwas(snp, input_type="rsID")***

Return the associated phenotypes and related research (from GWAS catelog) of the SNP.

```{r}
gwas_info <- get_snp_gwas("rs1891906")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/4_get_snp_gwas.png)

#### ***get_snp_eqtl(snp, input_type="rsID", eqtl_tissue = "Whole_Blood")***

Return GTEx significant single tissue eQTLs.

```{r}
gtex_data <- get_snp_eqtl("rs7524908", eqtl_tissue = "Whole_Blood")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/5_get_snp_eqtl.png)

#### ***get_snp_motif(snp, input_type="rsID")***

Return SNP-disrupted motifs.

```{r}
get_snp_motif("chr1:109676139", "hg38")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/6_get_snp_motif.png)

#### ***screenshot_GenomeBroswer(snp, input_type="rsID", filename, window_size=20)***

Return the Genome Browser screenshot of a SNP.

```{r}
screenshot_GenomeBroswer("rs456", filename="rs456.pdf")
screenshot_GenomeBroswer("chr22:19712094", input_type = "hg19", "~/test.pdf")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/7_genomebrowser.png)

#### ***get_LD(snps, pop="CEU", r2d="r2", assembly = "hg19", plot_heatmap=TRUE)***

Generates pairwise linkage disequilibrium statistics.

`snps` Required. between 1 - 10 variants, using an rsID or chromosome coordinate.

`pop` Optional. Default is "CEU". list_pop() return all the options.

`r2d` Optional. Default is "r2". Either "r2" for LD R2 (R-squared) or "d" for LD D' can be selected.

`assembly` Optional. Output UCSC Assembly versions. "hg19" or "hg38" can be selected. Default is "hg19".

`plot_heatmap` Optional. Whether to output the heatmap of the LD block. Default is TRUE.

```{r}
ld_block <- get_LD(c("chr7:24962419", "rs114", "rs127", "rs7805287", "rs60676332", "rs10239961"), assembly = "hg19")
```

Table:

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/8_ld_table.png)

Heatmap:

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/8_ld_heatmap.png)

#### ***plot_snp_circos(snp_info_table, window_size=2e5, output_assembly = "hg19", savefile="circos.pdf")***

Return the circos plot.

`snp_info_table` Required. The output data.frame of `query_snp()`.

`output_assembly` Optional. The version of chromosome coordinates of output circos plots. Default is "hg19".

`window_size` Optional. Window size around the SNP displayed. Default is 2e5. `savefile` Optional. The filename of the circos plot.

```{r}
info_table <- query_snp("rs1059196")
plot_snp_circos(info_table, savefile="circos.pdf")

info_table <- query_snp("rs1891906")
plot_snp_circos(snp_info_table=info_table, window_size = 1e5, savefile="~/test/circos2.pdf")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/9_circos.png)

#### ***hg19tohg38(snp_loc)***

Convert hg19 to hg38 coordinate.

#### ***hg38tohg19(snp_loc)***

Convert hg38 to hg19 coordinate.

#### ***get_genes_loc(genes, assembly="hg19")***

Return the genomic locations of gene symbols.

```{r}
get_genes_loc(get_genes_loc(c('GP1BB', 'TBX1')))
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/10_get_genes_loc.png)

#### **4.2 For batch SNPs**

***---Preprocess---***

Ensembl Variant Effect Predictor (VEP) is one of the most widely used Variant Annotation tools for the analysis, annotation, and prioritization of genomic variants in coding and non-coding regions. If there are many SNPs, it's better to annotate with VEP first. Below is a visualization of the VEP results.

Please annotate SNPs with the following parameters first with VEP:

```
vep --cache --dir_cache ~/SNP_visualize/  \
    --fasta ~/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    -i test_region.txt \
    -o output_region.txt \
    --var_synonyms \
    --af --af_1kg --af_esp \
    --everything \
    --force_overwrite
```

Then use the summary_vep.py to preprocess the VEP results:

```
python addition/summary.py -i VEP_annotation.txt -o summary.txt
```

Load data in R:

```{r}
data <- read.csv('file.txt', header=T, sep="\t") # Your preprocessed VEP annotation file
```

#### ***plot_batch_consequence(data, show_lines=7, colors="default")***

Return barplots of consequences.

`data` Required. The annotation results from VEP.

`show_num` Optional. The number of affected genes shown on the plot. Default is 7.

`colors` Optional. Default is `ggsci::pal_igv(alpha = 0.8)(show_num)`.

```{r}
plot_consequence(data)
# plot_consequence(data, show_num = 10, colors = rainbow(10))
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/plot_consequence.png)

#### ***plot_affect_gene(data, plot_type="all", show_num=7)***

Return barplots of affected genes of SNPs.

`data` Required. The annotation results from VEP.

`plot_type` Optional."gene", "snp", "all" and "merged" can be selected. Default is "rsID".

`show_num` Optional. The number of affected genes shown on the plot. Default is 7.

```{r}
pdf("affected_genes.pdf", height = 7, width = 8)
plot_affect_gene(data, plot_type = "merged")
dev.off()
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/plot_affect_gene.png)

#### ***plot_batch_titv(data)***

Return a barplot of ti/tv ratio.

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/titv.png)

#### ***plot_batch_AF(data, version="1000G")***

Return boxplots of allele frequency.

`data` Required. The annotation results from VEP.

`version` Optional."1000G" or "gnomAD" can be selected. Default is "1000G".

```{r}
plot_batch_AF(data, version = "1000G")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/AF.png)

#### ***plot_batch_cCRE(data, assembly = "hg38", show_unclassified=FALSE)***

Return barplots of cCREs that SNPs overlapped.

`data` Required. The annotation results from VEP.

`assembly` Optional. "hg19" and hg38" can be selected. Default is "hg38".

`show_unclassified` Optional. Whether to display SNPs without overlapping cCREs.

```{r}
plot_cCRE(data, plot_type="bar", assembly = "hg38", show_unclassified=FALSE)
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/batch_cCREs.png)


#### ***plot_batch_gwas(data, assembly="hg38", show_num=5)***

Return barplots of associated GWAS phenotypes of SNPs.

`data` Required. The annotation results from VEP.

`assembly` Optional. The assembly version of the input SNPs. "hg19" and "hg38" can be selected. Default is "hg38". Default is "hg38".

`show_num` Optional. The number of associated phenotypes shown on the plot. Default is 5.

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/batch_gwas.png)

### 5.**Data**

These datasets can be used directly.

cCRE_data_hg19: ENCODE SCREEN GRCh38-cCREs data.

cCRE_data_hg38: ENCODE SCREEN GRCh38-cCREs data.

gwas_data: GWAS catalog associations data.

chrom_info_hg19: chromosome infomation (hg19).

chrom_info_hg38: chromosome infomation (hg38).

