#' Theme

#' @return
#' @export
#'
#' @examples
#'

theme_snp <- function(){
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16), axis.title = element_text(size=14),
          axis.text = element_text(size = 13)) +
        theme(legend.text = element_text(size=12))
}
