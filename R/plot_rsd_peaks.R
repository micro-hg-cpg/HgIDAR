#' Check the quality of your integrations using the relative standard deviation of your peak's ratios
#' @description
#' A short description...
#'
#' @param x data frame with your integrated peaks
#' @param form IHg, MMHg, all
#' @param cutoff Numeric value. cutoff used for the plot. Default = 10%
#' @importFrom dplyr %>%
#' @return plot of the area peaks RSD
#' @export
#'
#' @examples
#' plot_rsd_peaks(peaks_1, form = 'MMHg')
#' plot_rsd_peaks(peaks_2, form = 'all', cutoff = 2)
#' plot_rsd_peaks(peaks_3, form = 'IHg', 5)
#'
plot_rsd_peaks=function(x, form = c('IHg', 'MMHg', 'all'), cutoff = 10){
  form <- match.arg(form)
  peaks1 <- x[, colSums(x != 0) > 0]
  peakscols1 <- ncol(peaks1)


  if (form == 'IHg') {

      peaks_IHg <- peaks1 %>%
        dplyr::select(1:2 | contains("IHg"))

      a <- peaks_IHg%>%
        dplyr::group_by(sample)%>%
        dplyr::summarise_all(rtemisalpha::rsd)%>%
        dplyr::select(-c(injection))

      b <- reshape2::melt(a, id.vars = 'sample')

      plot <- ggplot2::ggplot(b)+
        ggplot2::geom_point(ggplot2::aes(sample, value))+
        ggplot2::geom_hline(yintercept=cutoff)+
        ggplot2::scale_x_continuous(name="sample")+
        ggplot2::scale_y_sqrt()+
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_line(color = "gray"))+
        ggplot2::ylab("RSD peaks")+
        ggplot2::ggtitle("IHg")+
        ggplot2::facet_wrap(~variable)

    }
  if (form == 'MMHg'){
      peaks_MMHg <- peaks1 %>%
        dplyr::select(1:2 |contains("MMHg"))

      a <- peaks_MMHg%>%
        dplyr::group_by(sample)%>%
        dplyr::summarise_all(rtemisalpha::rsd)%>%
        dplyr::select(-c(injection))

      b <- reshape2::melt(a, id.vars = 'sample')

      plot <- ggplot2::ggplot(b)+
        ggplot2::geom_point(ggplot2::aes(sample, value))+
        ggplot2::geom_hline(yintercept=cutoff)+
        ggplot2::scale_x_continuous(name="sample")+
        ggplot2::scale_y_sqrt()+
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_line(color = "gray"))+
        ggplot2::ylab("RSD peaks")+
        ggplot2::ggtitle("MMHg")+
        ggplot2::facet_wrap(~variable)

    }
  if (form == 'all') {

      a <- peaks1%>%
        dplyr::group_by(sample)%>%
        dplyr::summarise_all(rtemisalpha::rsd)%>%
        dplyr::select(-c(injection))

      b <- reshape2::melt(a, id.vars = 'sample')

      plot <- ggplot2::ggplot(b)+
        ggplot2::geom_point(ggplot2::aes(sample, value))+
        ggplot2::geom_hline(yintercept=cutoff)+
        ggplot2::scale_x_continuous(name="sample")+
        ggplot2::scale_y_sqrt()+
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_line(color = "gray"))+
        ggplot2::ylab("RSD peaks")+
        ggplot2::ggtitle("IHg and MMHg")+
        ggplot2::facet_wrap(~variable)

  }

  return(plot)
}

