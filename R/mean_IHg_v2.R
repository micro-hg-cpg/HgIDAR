#' Calculates concentration of the inorganic (IHg) isotopes in sediments
#' @description
#' The mean concentration of IHg and their standard deviations of the different isotopes are returned in a table. Mean and sd are calculated with the injection values.
#'
#' @param x data frame with your integrated peaks
#' @param y data frame with your masses
#' @param z data frame with your abundances
#' @param samples Numeric vector. Number of samples
#' @param injections Numeric vector. Number of injections per sample Default = 3
#' @param type Type of sample: liquid or solid
#' @param quant_IHg IHg isotope used for the quantification
#' @param spike_concentration_IHg Concentration of the organic Hg quantification spike. Units: ng/g
#' @return Table with the mean sample concentrations and the standard deviations for each isotope.
#' @importFrom dplyr %>%
#' @export
#'
#'
#' @examples
#' mean_IHg(
#' x=peaks_3,
#'   y=masses_3,
#'   z=abundances_3,
#'   samples = 96,
#'   type = 'solid',
#' quant_IHg = 198,
#' spike_concentration_IHg = 5)
#' mean_IHg(peaks_1, masses_1, abundances_1, 24, 3, type = 'liquid' , 198, 5)
#' mean_IHg(peaks_2, masses_2, abundances_2, 26, 3, type = 'solid' ,199, 27.27)
mean_IHg=function(x, y, z, samples, injections = 3, type = c('liquid', 'solid'),   quant_IHg, spike_concentration_IHg ){
  type <- match.arg(type)
  #set the sample interval
  Samples <- c(1:samples)
  #set the injections that you made
  Injection <- c(1:injections)
  #select the peaks of inorganic
  peaks <- x %>%
    dplyr::select(1:2 |contains("IHg"))
  #remove columns with all 0 from your data, in case that you have any.
  peaks1 <- peaks[, colSums(peaks != 0) > 0]
  #calculate the inverse  of the ratios matrix obtained from the abundances
  ratios_inv <- abun_treat(z)
  #calculate the molecular weight of each of your isotopes
  M.w <- M.w(z)
  #select the molecular mass of the quantification spike
  mw1 <- as.data.frame(M.w)%>%
      tibble::rownames_to_column()
  mw2 <- mw1[grepl(paste(quant_IHg), mw1$rowname),]
  #calculate the added mols of spike and join it to the other masses used.
  mol_added <- (spike_concentration_IHg*y$spike_IHg)/mw2[1,2]
  masses <- cbind(y, mol_added)
  #create a minidataframe with the molecular weighs and abundance of your mayor isotope.
  #This is located in the excel in the top left corner.
  mw3 <-  tibble::rownames_to_column(as.data.frame(M.w))
  mayor.isotope <- as.data.frame(apply(z, 2, max))%>%
      tibble::rownames_to_column()
  mayor.iso.mw <- merge(mw3, mayor.isotope, by = 'rowname', all.x = T)%>%
      tibble::column_to_rownames(var = 'rowname')
  colnames(mayor.iso.mw) <- c('M.w', 'mayor.isotope')
  #made the first step of the calculations, the deconvulation matrix.
  matrix.deconvolution.results <- deconvulation_matrix(peaks1, ratios_inv, samples = samples, injections = injections)
  #calculation of the concentration in the extract or in the water if its liquid.
  liquid <- five_isotopes_2(decon = matrix.deconvolution.results, masses,
                            quant = quant_IHg, samples = samples, injections = injections,
                            mayor.iso.mw = mayor.iso.mw)
  concentration <- data.frame("sample"=liquid$sample,
                                "injection"=liquid$injection,
                                "isotope_IHg"=liquid$variable,
                                "concentration_IHg"=liquid$conc_ext)
  o <- tidyr::pivot_wider(data=concentration, names_from = isotope_IHg, values_from = concentration_IHg)
  p <- o%>%
      dplyr::group_by(sample)%>%
      dplyr::summarise_all(list(mean, stats::sd))%>%
      dplyr::select(-c('injection_fn1', 'injection_fn2'))


  if (type == 'solid') {
    solid <- data.frame()
    for (i in Samples) {
      ext_mass <- masses%>%
        dplyr::filter(sample ==i)
      for (j in Injection) {
        h <- liquid %>%
          dplyr::filter(sample == i & injection == j)
        h <- h%>%
          dplyr::mutate("conxaci"=conc_ext*ext_mass$acid)
        h<- h%>%
          dplyr::mutate("conc_sed"=conxaci/ext_mass$sed)
        solid <- rbind(solid, h)

        concentration <- data.frame("sample"=solid$sample,
                                    "injection"=solid$injection,
                                    "isotope_IHg"=solid$variable,
                                    "concentration_IHg"=solid$conc_sed)

        o <- tidyr::pivot_wider(data=concentration, names_from = isotope_IHg, values_from = concentration_IHg)
        p <- o%>%
          dplyr::group_by(sample)%>%
          dplyr::summarise_all(list(mean, stats::sd))%>%
          dplyr::select(-c('injection_fn1', 'injection_fn2'))
      }
    }
  }

    final_output <- data.frame(p[,1],
                               p[,2], p[,7],
                               p[,3], p[,8],
                               p[,4], p[,9],
                               p[,5], p[,10],
                               p[,6], p[,11])
    colnames(final_output)  <- gsub('fn1', 'mean', colnames(final_output))
    colnames(final_output)  <- gsub('fn2', 'sd', colnames(final_output))


  return(final_output)

}

