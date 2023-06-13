## Vi-SNP: An R-package to query and visualize SNPs comprehensively


## 1. **Introduction**

Single-nucleotide polymorphism (SNP) is an important topic for most genetic association studies. This package can query, analyze, visualize and link SNPs with functional information (e.g., genes, regulatory elements, phenotypes) from multiple databases (3DSNPv2, GTEx, SCREEN and GWAS catalog).

## 2. **Installation**

``` bash
devtools::install_github("Liying1996/ViSNP")
```

## 3. **Functions of the package**

Workflow.
![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/workflow.jpg)

## 4. **Details of functions**

### **4.1 For single SNP**

#### ***query_snp(snp, input_type="rsID", eqtl_tissue="Whole_Blood", verbose=FALSE)***

`query_snp` output the basic and functional information of a SNP.

`input_type`: "rsID", "hg19" or "hg38". If "hg19" or "hg38" was selected, the ***snp*** should be "chr22:19712094". (All genomic locations should be 1-based.) 
`eqtl_tissue`: Tissue ID of the tissue of GTEx. All the tissues: https://gtexportal.org/api/v2/redoc#tag/Static-Association-Endpoints/operation/get_significant_single_tissue_eqtls_by_location_api_v2_association_singleTissueEqtlByLocation_get tissueSiteDetailId.

`verbose`: Whether to show the intermediate process. Default is FALSE.




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

`snps` rsIDs or genomic locations. 

`colors`: Colors of lines. If not specified, default is to use the colors of `rainbow()`.

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
get_snp_loc("rs10040658", 'hg38')
```

#### ***get_snp_ccre(snp, input_type="rsID", output_assembly="hg19")***

Return the ENCODE SCREEN (https://screen.encodeproject.org/) cCREs (candidate cis-regulatory elements) that SNPs overlapped.

```{r}
snp <- "rs1059196"
cCRE_intersect <- get_snp_ccre(snp)
# get_snp_ccre("chr22:19724571", input_type = "hg38", output_assembly = "hg38")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/3_get_snp_ccre.png)

#### ***get_snp_gwas(snp, input_type="rsID", output_type="core")***

Return the associated phenotypes and related research (from GWAS catalog) of the SNP.  

`output_type`: "core" or "full" can be selected. The output of "core" is subset of "full". Default is "core".

```{r}
gwas_info <- get_snp_gwas("rs1891906")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/4_get_snp_gwas.png)

#### ***get_snp_eqtl(snp, input_type="rsID", eqtl_tissue = "Whole_Blood")***

Return GTEx significant single tissue eQTLs.

```
gtex_data <- get_snp_eqtl("rs10040658", eqtl_tissue = "Whole_Blood")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/5_get_snp_eqtl.png)

#### ***get_snp_motif(snp, input_type="rsID")***

Return SNP-disrupted motifs.

```{r}
get_snp_motif("chr1:109676139", "hg38")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/6_get_snp_motif.png)



***plot_snp_motif(snp, input_type="rsID", tf='all')***

```
plot_snp_motif(snp="rs10040658")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/6_plot_snp_motif.jpg)

```
plot_snp_motif(snp="rs10040658", tf=c('AP2alpha', 'CHCH_01'))
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/6_plot_snp_motif2.jpg)

#### ***get_loop_gene(snp, input_type="rsID",  output_type="core")***

Return 3D-interacting genes of a SNP.

`output_type`: "core" or "full" can be selected. The output of "core" is subset of "full". Default is "core".

```{R}
loop_genes <- get_loop_gene("rs1059196", output_type="full")
```
![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/11_loop_genes.jpg)

***plot_loop_gene(snp, input_type="rsID", output_assembly="hg19", show_cells=3, curvature=0.8)***

`output_assembly` Optional. The UCSC Assembly versions of chromosome coordinates of outputs. "hg19" or "hg38" can be selected. Default is "hg19".

`show_cells` Optional. The number of cell types shown in the figure. Default is 3 (Top 3 cell types). The names of cell types are available as well.

`curvature` Optional. The "curvature" parameter of `geom_curve`. Default is "0.8".

```
plot_loop_gene(snp="rs10040658", output_assembly='hg19', show_cells=c("GM12878", "K562"))
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/11_plot_loop_gene2.jpg)

```
plot_loop_gene(snp="rs4942486", output_assembly='hg19', show_cells=c("VentricleLeft", "Caki2", "HepG2"), curvature=-0.3)
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/11_plot_loop_gene.jpg)

#### ***get_loop_snp(snp, input_type="rsID",  output_type="core")***

Return 3D-interacting SNPs of a SNP.

`output_type`: "core" or "full" can be selected. The output of "core" is subset of "full". Default is "core".

```{R}
loop_snps <- get_loop_snp("rs1059196")
```
![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/11_loop_snps.jpg)

***plot_loop_snp(snp, input_type="rsID", pop, output_assembly="hg19", show_cells=3)***

`pop` Required. Population.

`output_assembly` Optional. The UCSC Assembly versions of chromosome coordinates of outputs. "hg19" or "hg38" can be selected. Default is "hg19".

`show_cells` Optional. The number of cell types shown in the figure. Default is 3 (Top 3 cell types). The names of cell types are available as well.

`curvature` Optional. The "curvature" parameter of `geom_curve`. Default is "0.8".

```
plot_loop_snp(snp='rs10040658', pop = 'AFR', show_cells = c("GM12878", "Lung", "Caki2"))

```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/11_plot_loop_snp2.jpg)

```
plot_loop_snp(snp='rs10', pop = 'AFR', show_cells = c("Bladder", "Spleen"))
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/11_plot_loop_snp.jpg)

#### ***screenshot_GenomeBroswer(snp, input_type="rsID", filename, window_size=20)***

Return the Genome Browser screenshot of a SNP.

```{r}
screenshot_GenomeBroswer("rs456", filename="rs456.pdf")
screenshot_GenomeBroswer("chr22:19712094", input_type = "hg19", "~/test.pdf")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/7_genomebrowser.png)

#### ***get_LD(snps, pop="CEU", r2d="r2", assembly = "hg19", plot_heatmap=TRUE)***

Generates pairwise linkage disequilibrium statistics.

**It should be noted that `LDlinkR` is embeded in this function, so a personal access token is required in In order to access the LDlinkAPI via `LDlinkR`.  For details, please see https://rdrr.io/cran/LDlinkR/f/vignettes/LDlinkR.Rmd or you can follow the steps below:**

1. Make a one-time request for your personal access token from a web browser at https://ldlink.nci.nih.gov/?tab=apiaccess.
2. Once registered, your personal access token will be emailed to you. It is a string of 12 random letters and numbers.
3. Then,  from the R console, do (to open the .Renviron file):

```
usethis::edit_r_environ() # Open the .Renviron file

```

4. Next, add a line that looks like this in the .Renviron file:

```
LDLINK_TOKEN=YourTokenHere123
```

5. Finally, save and close the `.Renviron` file. Restart R, as environment variables are only loaded from `.Renviron` at the start of a new R session. Now, you can check by entering:

```
Sys.getenv("LDLINK_TOKEN")
```



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

```
info_table <- query_snp("rs10040658")
plot_snp_circos(snp_info_table=info_table, savefile="~/test/circos3.pdf")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/9_circos2.jpg)



***plot_snp_network(snp, input_type = "rsID", cells = "all", pop, tissue_type="Whole_Blood", colors="default")***

`colors` Optional. Vectors with 6 colors. "default" is `brewer.pal(6, "Set1")`.

```
plot_snp_network(snp="rs10040658", pop='EAS', cells=c("A549", "K562", "Caki2", "GM12878", "LNCAP"))
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/12_plot_networks2.jpg)

```
plot_snp_network(snp="rs1059196", pop='EAS')
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/12_plot_networks.jpg)

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

##### [Annotated via VEP API]

***get_batch_vep(snps, input_type='rsID')***

Ensembl Variant Effect Predictor (VEP) is one of the most widely used Variant Annotation tools for the analysis, annotation, and prioritization of genomic variants in coding and non-coding regions. 

`snps` Required.
`input_type` Optional. The type of the input SNP. "rsID", hg19" or "hg38" can be selected. Default is "rsID".

```
vep_anno <- get_batch_vep(c("rs56116432", "rs10040658"), input_type = "rsID")

# get_batch_vep(snps=c("chr9:136131429", "chr5:139051015"), input_type = 'hg19')

# get_batch_vep(snps=c("chr9:133256042", "chr5:139671430"), input_type = 'hg38')
```

##### [Annotated via `query_snp`] 

***get_batch_SNP_info(snps, input_type="rsID", eqtl_tissue="Whole_Blood")***

Return information of a batch of SNPs (Annotated via `query_snp`. This method is relatively slow due to the need to access multiple databases via APIs.

`snp` Required.
`input_type` Optional. The assembly version of the input SNP. "rsID", "hg19" and "hg38" can be selected. Default is "rsID".
`eqtl_tissue` Optional. Tissue ID of the tissue of interest. Default is "Whole_Blood".

```
batch_info_table <- get_batch_SNP_info(snps = c("rs1891906", "rs10"), input_type="rsID", eqtl_tissue="Whole_Blood")

# batch_info_table <- get_batch_SNP_info(snps = c("chr1:1014863", "chr7:92754574"), input_type="hg38", eqtl_tissue="Whole_Blood")
```

##### [Annotated by VEP (Linux commands - Uploaded by users)]

***---Preprocess---***

 If there are many SNPs, it's better to annotate with VEP command lines. Below is the preprocessor of  Linux command lines.

Please annotate SNPs with the following parameters first with VEP:

```
# Manually downloading caches
# cd ~/SNP_visualize/
# curl -O https://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_vep_105_GRCh38.tar.gz
# tar xzf homo_sapiens_vep_109_GRCh38.tar.gz

vep --cache --dir_cache ~/SNP_visualize/  \
    --fasta ~/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    -i test.vcf \  # Your VCF file
    -o VEP_annotation.txt \
    --var_synonyms \
    --af --af_1kg --af_esp \
    --everything \
    --force_overwrite \
    --offline
```

Then use the addition/summary_vep.py to preprocess the VEP results:

```
python addition/summary_vep.py -i VEP_annotation.txt -o summary.txt
```

Load data in R:

```{r}
# eg.
data <- read.csv('summary.txt', header=T, sep="\t") # Your preprocessed VEP annotation file
```

**[Test Data]**

```
# SNP list 
test_snps <- test_snps # The dataset can be used directly.
snp_list <- test_snps[,1]

# Annotated data from get_batch_vep() for the following analyses and visualization
test_api <- get_batch_vep(snp_list, input_type = "rsID")

# Uploaded example data from VEP (Linux command)
test_upload <- test_upload  # The dataset can be used directly.
```

**[Functions for analyses and visualization] **

#### ***plot_batch_consequence(data, show_num=7, data_source="API", colors="default")***

Return barplots of consequences.

`data` Required. The annotation results from VEP.

`show_num` Optional. The number of affected genes shown on the plot. Default is 7.

`data_source` Optional. The source of the input data. "API" or "Upload" can be selected. Default is "API" (Input annotated data from get_batch_vep() function).

`colors` Optional. Default is `ggsci::pal_igv(alpha = 0.8)(show_num)`.

```{r}
plot_batch_consequence(test_upload, data_source="Upload")
# plot_batch_consequence(test_api, show_num = 10, data_source="API", colors = rainbow(10))
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/plot_consequence.png)

***analyze_gene_go(snp_genes, output="dotplot")***

 Return the go enrichment of genes.

`snp_genes` Required. The annotated genes of input SNPs.
`output_type` Optional. Show dotplot/barplot of GO enrichment results. Default is dotplot.

```
genes <- unique(test_upload$Gene)
analyze_gene_go(genes)
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/analyze_go.png)

#### ***plot_affect_gene(data, data_source="API", gene="Ensembl", plot_type="merged", show_num=7, go_enrichment=FALSE)***

Return barplots of affected genes of SNPs.

`data` Required. The annotation results from VEP.

`data_source` Optional. The source of the inpit data. "API" or "Upload" can be selected. Default is "API" (Input annotated data from get_batch_vep() function).

`gene` Optional. The version of gene names. "Ensembl" and "Symbol" can be selected.  Default is "Ensembl".

`plot_type` Optional."gene", "snp", "all" and "merged" can be selected. Default is "merged".

`show_num` Optional. The number of affected genes shown on the plot. Default is 7.

`go_enrichment` Optional. Whether to do the GO enrichment analysis of genes. Default is FALSE. (Users can do the enrichment analysis by analyze_gwas_enrich() as well)



```{r}
pdf("affected_genes.pdf", height = 7, width = 8)
plot_affect_gene(data=test_upload, data_source="Upload", gene = "Symbol", go_enrichment=FALSE) # If go_enrichment is TRUE, the dotplot generated by analyze_gene_go() is returned as well
# plot_affect_gene(data=test_api, data_source="API", go_enrichment=TRUE) 
dev.off()
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/plot_affect_gene.jpg)

#### ***plot_batch_titv(data)***

Return a barplot of ti/tv ratio.

`data` The annotation results from VEP. (Uploaded VEP version only)

```
plot_batch_titv(test_upload)
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/titv.png)

#### ***plot_batch_AF(data, version="1000G")***

Return boxplots of allele frequency.

`data` Required. The annotation results from VEP. (Uploaded VEP version only)

`version` Optional."1000G" or "gnomAD" can be selected. Default is "1000G".

```{r}
plot_batch_AF(test_upload, version = "1000G")
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/AF.png)



***analyze_ccre_enrich(snps_loc, assembly = "hg38", show_p=FALSE)***

**R package `bedtoolsr` is embeded in this function and it uses bedtools to overlap, so when using `analyze_ccre_enrich`and `plot_batch_cCREs`, the user must configure the path of bedtools, for example:**

```
options(bedtools.path = "~/anaconda3/bin/")
```

Simulation enrichment  is used to test whether the input SNPs are enriched in cCREs. 

Method:

1). Randomly select the same number of SNPs as input SNPs from 1000Genome as control ( 1 million SNPs randomly selected from 1000 Genome were used and embeded in Vi-SNP in order to speed up the test), iterate 500 times;

2). Calculate the number of input and control SNPs enriched in each cCRE;

3). For each cCRE, calculate the Z-score. Z-score = (the number of input SNPs enriched in the cCRE - mean number of control SNPs enriched in the cCRE) / (standard deviation of control SNPs enriched in the cCRE);

4). Calculate P-values. P = 2 * pnorm(-Zscore)

`snps_loc` Required. The genomic locations of input SNPs, the format should be like: "chr1:1014863".
`assembly` Optional. "hg19" and hg38" can be selected. Default is "hg38".
`show_p` Optional. Show p-values or "\*\*\*" (significant) on output plot. Default is TRUE. "\*\*\*" : P <= 1e-7; "\*\*": 1e-3 <= P < 1e-7;  "\*": 0.05 <= p < 1e-3; "NS": p > 0.05.

```
snps_loc <- unique(test_upload$Location)
analyze_ccre_enrich(snps_loc, assembly='hg38')

# Use API data
snps_loc <- unique(paste('chr', test_api$seq_region_name, ':', test_api$start, sep = ''))
snps_loc <- snps_loc[grep(pattern="CHR",snps_loc, invert = TRUE)]
analyze_ccre_enrich(snps_loc, assembly='hg38')
```

If p = TRUE:

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/analyze_ccre_1.png)



If p = FALSE:

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/analyze_ccre_2.png)



#### ***plot_batch_cCREs(snps_Loc, assembly = "hg38", show_unclassified=FALSE)***

Return barplots of cCREs that SNPs overlapped.

`data` Required. The annotation results from VEP.

`assembly` Optional. "hg19" and hg38" can be selected. Default is "hg38".

`show_unclassified` Optional. Whether to display SNPs without overlapping cCREs.

`enrichment` Optional. Whether to do the enrichment analysis of cCREs. Default is FALSE. (Users can do the enrichment analysis by analyze_ccre_enrich() as well)
`show_p` Optional.Show p-values or "***" (significant) on enrichment plot. Default is FALSE.



```{r}
plot_batch_cCREs(snps_loc, assembly = "hg38", show_unclassified=FALSE, enrichment=FALSE, show_p=FALSE)

# plot_batch_cCREs(snps_loc, assembly = "hg38", show_unclassified=FALSE, enrichment=TRUE, show_p=FALSE)
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/batch_cCREs.png)



***analyze_gwas_enrich(snps_loc, assembly='hg38')***

Return the enrichment results of GWAS studies.

`snps_loc` Required. The locations of SNPs, the format should be "chr1:1014863".
`assembly` Optional. "hg19" and hg38" can be selected. Default is "hg38".

```
snps_loc <- test_snps[,2][450:1000]
analyze_gwas_enrich(snps_loc)
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/analyze_gwas.png)

#### ***plot_batch_gwas(snps_loc, assembly="hg38", show_num=5,  enrichment=FALSE)***

Return barplots of associated GWAS phenotypes of SNPs.

`data` Required. The annotation results from VEP.

`assembly` Optional. The assembly version of the input SNPs. "hg19" and "hg38" can be selected. Default is "hg38". Default is "hg38".

`show_num` Optional. The number of associated phenotypes shown on the plot. Default is 5.

`enrichment` Optional. Whether to do the enrichment analysis of cCREs. Default is FALSE. (Users can do the enrichment analysis by analyze_ccre_enrich() as well)
`show_p` Optional.Show p-values or "***" (significant) on enrichment plot. Default is FALSE.

```
plot_batch_gwas(snps_loc, assembly = "hg38", show_num=5, enrichment=FALSE)
```

![](https://github.com/Liying1996/ViSNP/blob/master/example_figs/batch_gwas.png)

### 5.**Data**

These datasets can be used directly.

cCRE_data_hg19: ENCODE SCREEN GRCh19-cCREs data.

cCRE_data_hg38: ENCODE SCREEN GRCh38-cCREs data.

gwas_data: GWAS catalog associations data.

chrom_info_hg19: chromosome infomation (hg19).

chrom_info_hg38: chromosome infomation (hg38).

control_samples: The 1000,000 control SNPs for simulation enrichment randomly downsampled from 1000 Genome Project.

test_snps: example SNPs list.

test_upload: example annotation file from VEP linux command and preprossed by addition/summary_vep.py.

