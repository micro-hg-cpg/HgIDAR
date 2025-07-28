deconvulation_matrix = function(x, ratios_inv, samples , injections ){

  Samples <- c(1:samples)
  Injection <- c(1:injections)
  tidy.peaks <- reshape2::melt(x, id=c("sample", "injection"))
  tidy.peaks <- tidy.peaks %>%
    dplyr::arrange(injection)%>%
    dplyr::arrange(sample)
  first.results <- data.frame()
  matrix.deconvolution.results <- data.frame()

  for (i in Samples) {
    for (j in Injection) {
      a <- x %>%
        dplyr::filter(sample == i & injection == j)
      col.peaks <- ncol(a)
      mat_dec <- ratios_inv %*% t(as.matrix(a[,3:col.peaks]))
      first.results <- rbind(first.results, mat_dec)
    }
  }
  matrix.deconvolution.results <- cbind(tidy.peaks, first.results)
  return(matrix.deconvolution.results)
}

