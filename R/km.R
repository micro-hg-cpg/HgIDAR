#' Calculations of Hg the methylation rates
#'
#' @param x data frame with your integrated peaks
#' @param y data frame with your masses
#' @param z data frame with your abundances
#' @param samples Numeric vector. Number of samples
#' @param injections Numeric vector. Number of injections per sample Default = 3
#' @param type Type of sample: liquid or solid
#' @param incub_IHg IHg isotope used for the incubation
#' @param incub_MMHg MMHg isotope used for the incubation
#' @param quant_MMHg MMHg isotope used for the quantification
#' @param quant_IHg IHg isotope used for the quantification
#' @param nat isotope used as natural for the incubation
#' @param spike_concentration_MMHg Concentration of the organic Hg quantification spike. Units: ng/g
#' @param spike_concentration_IHg Concentration of the inorganic Hg quantification spike. Units: ng/g
#' @importFrom dplyr %>%
#' @return Table with the sample mean Hg methylation rates and the sd.
#' @export
#'
#' @examples
#' km(peaks_3,
#'   masses_3,
#'   abundances_3,
#'   samples = 96,
#'   injections =3,
#'   type = 'solid',
#'   incub_IHg=199,
#'   incub_MMHg=201,
#'   quant_MMHg=202,
#'   quant_IHg=198,
#'   nat=200,
#'   spike_concentration_MMHg=0.22,
#'   spike_concentration_IHg=5)
#'  km(peaks_1,
#'   masses_1,
#'   abundances_1,
#'   samples = 24,
#'   injections = 3,
#'   type = 'liquid',
#'   incub_IHg= 199,
#'   incub_MMHg=201,
#'   quant_MMHg=202,
#'   quant_IHg=198,
#'   nat=200,
#'   spike_concentration_MMHg=0.22,
#'   spike_concentration_IHg=5)
km=function(x, y, z, samples, injections = 3, type = c('liquid', 'solid'), incub_IHg, incub_MMHg, quant_MMHg, quant_IHg, spike_concentration_MMHg, nat, spike_concentration_IHg ){
  type <- match.arg(type)

  resultsIHg <- IHg(x, y, z, samples = samples , injections = injections, type=type,incub_IHg=incub_IHg,incub_MMHg=incub_MMHg, quant_MMHg=quant_MMHg, quant_IHg=quant_IHg,nat=nat, spike_concentration_IHg = spike_concentration_IHg)
  resultsMMHg <- MMHg(x, y, z, samples = samples, injections = injections, type=type, incub_IHg=incub_IHg,incub_MMHg=incub_MMHg, quant_MMHg=quant_MMHg, quant_IHg=quant_IHg,nat=nat, spike_concentration_MMHg = spike_concentration_MMHg)


  concentration <- merge(resultsIHg, resultsMMHg, by = c('sample', 'injection'))

  isotopes <- list('196' , '198', '199',
                   '200', '201', '202' ,
                   '204')

  for (key in isotopes){
    key_numeric <- as.numeric(key)
    if(key_numeric == incub_IHg){
      print(paste("IHg incubation spike:", key))
      w <- concentration%>%
        dplyr::select(1:2 | contains(key))
    }
  }

  p <- w%>%
    dplyr::mutate("totalconcentration" = w[,3]+w[,4])%>%
    dplyr::mutate("methylation" = (w[,4]/totalconcentration))%>%
    dplyr::select(-c('totalconcentration'))


  final_output <- data.frame(p[,1],
                             p[,2],
                             p[,3],
                             p[,4],
                             p[,5])

  colnames(final_output)  <- c('sample',
                               'injection',
                               'IHg',
                               'MMHg',
                               'km')

  return(final_output)

}
