#' Check the quality of your integrations using the relative standard deviation of your peak's ratios
#' @description
#' A short description...
#'
#' @param x data frame with your integrated peaks
#' @param form IHg, MMHg, all
#' @param cutoff Numeric value. cutoff used for the plot. Default = 10%
#' @importFrom dplyr %>%
#' @return table with the area peaks RSD
#' @export
#'
#' @examples
#' rsd_peaks(peaks_1, form = 'MMHg')
#' rsd_peaks(peaks_2, form = 'all', cutoff = 2)
#' rsd_peaks(peaks_3, form = 'IHg', 5)
rsd_peaks=function(x, form = c('IHg', 'MMHg', 'all'), cutoff = 10){
  form <- match.arg(form)
  peaks1 <- x[, colSums(x != 0) > 0]

  if (form == 'IHg') {

    peaks_IHg <- peaks1 %>%
      dplyr::select(1:2 |contains("IHg"))
#cambiar los nombres de las funciones.
    a <- peaks_IHg%>%
      dplyr::group_by(sample)%>%
      dplyr::summarise_all(list(mean, stats::sd, TIGERr::compute_RSD()))%>%
      dplyr::select(-c('injection_fn1', 'injection_fn2', 'injection_fn3'))
    b <- peaks_IHg%>%
      dplyr::left_join(a)

  }
  if (form == 'MMHg'){
    peaks_MMHg <- peaks1 %>%
      dplyr::select(1:2 |contains("MMHg"))

    a <- peaks_MMHg%>%
      dplyr::group_by(sample)%>%
      dplyr::summarise_all(list(mean, stats::sd, TIGERr::compute_RSD()))%>%
      dplyr::select(-c('injection_fn1', 'injection_fn2', 'injection_fn3'))
    b <- peaks_MMHg%>%
      dplyr::left_join(a)

  }
  if (form == 'all') {

    a <- peaks1%>%
      dplyr::group_by(sample)%>%
      dplyr::summarise_all(list(mean, stats::sd, TIGERr::compute_RSD()))%>%
      dplyr::select(-c('injection_fn1', 'injection_fn2', 'injection_fn3'))
    b <- peaks1%>%
      dplyr::left_join(a)

  }

  return(b)
}

