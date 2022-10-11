#' Return the Genome Browser screenshot of a SNP.
#'
#' @param snp Required.
#' @param input_type Optional. The type of the input SNP. "rsID", hg19" or "hg38" can be selected. Default is "rsID".
#' @param filename Optional. The filename of the output screenshot.
#' @param window_size Optional. The window size around the SNP.
#'
#' @return A PDF screenshot.
#' @export
#'
#' @examples
#' screenshot_GenomeBroswer("chr22:19712094", input_type = "hg19", "~/test.pdf")
#'



screenshot_GenomeBroswer <- function(snp, input_type="rsID", filename="screenshot.pdf", window_size=20){
  if (input_type=="rsID"){
    snp <- get_snp_loc(snp, assembly = "hg38")
  }

  chrom <- str_split_fixed(snp, ":", 2)[,1]
  pos <- as.numeric(str_split_fixed(snp, ":", 2)[,2])
  pos1 <- pos + window_size
  pos2 <- pos - window_size

  if (input_type!="hg19"){
    url_page <- paste("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=", chrom, "%3A", pos1, "%2D", pos2, "&hgsid=1280864071_rnlawaLfUQQ1i4quzzQsWWEDlzt1&hgt.psOutput=on", sep="")}else{
    url_page <- paste("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=", chrom, "%3A", pos1, "%2D", pos2, "&hgsid=1280864071_rnlawaLfUQQ1i4quzzQsWWEDlzt1&hgt.psOutput=on", sep="")
  }
  web <- readLines(url_page, encoding = "UTF-8")
  namerow <- grep("the current browser graphic in PDF", web)
  pdf_url <- paste("https://genome.ucsc.edu/", gsub(pattern = "\\.\\.\\/", replacement = "", x = str_split_fixed(web[namerow], '"', 3)[,2]), sep="")
  download(pdf_url, filename, mode="wb")
}
