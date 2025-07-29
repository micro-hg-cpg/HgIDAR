#' Calculates concentration of the inorganic (MMHg) isotopes in sediments
#' @description
#' The mean concentration of MMHg and their standard deviations of the different isotopes are returned in a table. Mean and sd are calculated with the injection values.
#'
#' @param x data frame with your integrated peaks
#' @param y data frame with your masses
#' @param z data frame with your abundances
#' @param samples Numeric vector. Number of samples
#' @param injections Numeric vector. Number of injections per sample Default = 3
#' @param type Type of sample: liquid or solid
#' @param incub_IHg IHg isotope used for the incubation
#' @param incub_MMHg MMHg isotope used for the incubation
#' @param quant_IHg IHg isotope used for the quantification
#' @param quant_MMHg MMHg isotope used for the quantification
#' @param nat isotope used as natural for the incubation
#' @param spike_concentration_MMHg Concentration of the organic Hg quantification spike. Units: ng/g
#' @return Table with the mean sample concentrations and the standard deviations for each isotope.
#' @importFrom dplyr %>%
#' @export
#' @examples
#' MMHg(
#' x=peaks_3,
#'   y=masses_3,
#'   z=abundances_3,
#'   samples = 96,
#'   type = 'solid',
#' incub_IHg = 199,
#' incub_MMHg=201,
#' quant_MMHg=202,
#' quant_IHg=198,
#' nat=200,
#' spike_concentration_MMHg = 0.5)
#' MMHg(peaks_1, masses_1, abundances_1, 24, 3, type = 'liquid' ,199,201,198,202,200, 0.22)
#' MMHg(peaks_2, masses_2, abundances_2, 26, 3, type = 'solid' ,204,200,199,201,202,0.23)


MMHg=function(x, y, z, samples, injections = 3, type = c('liquid', 'solid') ,incub_IHg,
              incub_MMHg=NULL,quant_IHg,quant_MMHg,nat, spike_concentration_MMHg ){
  type <- match.arg(type)
  #set the sample interval
  Samples <- c(1:samples)
  #set the injections that you made
  Injection <- c(1:injections)
  #select the peaks of inorganic
  peaks <- x %>%
    dplyr::select(1:2 |contains("MMHg"))
  #remove columns with all 0 from your data, in case that you have any.
  peaks1 <- peaks[, colSums(peaks != 0) > 0]
  #calculate the inverse  of the ratios matrix obtained from the abundances
  ratios_inv <- abun_treat(z)
  #calculate the molecular weight of each of your isotopes
  M.w <- M.w(z)
  #select the molecular mass of the quantification spike
  mw1 <- as.data.frame(M.w)%>%
    tibble::rownames_to_column()
  mw2 <- mw1[grepl(paste(quant_MMHg), mw1$rowname),]
  #calculate the added mols of spike and join it to the other masses used.
  mol_added <- (spike_concentration_MMHg*y$spike_MMHg)/mw2[1,2]
  masses <- cbind(y, mol_added)
  #create a minidataframe with the molecular weighs and abundance of your mayor isotope.
  #This is located in the excel in the top left corner.
  mw3 <-  tibble::rownames_to_column(as.data.frame(M.w))%>%
    tidyr::separate(rowname, c('no','isotope'), sep = '_')%>%
    dplyr::select(-no)
  mayor.isotope <- mayor.isotope.fn(z,incub_IHg=incub_IHg,incub_MMHg=incub_MMHg,
                                    quant_IHg=quant_IHg, quant_MMHg=quant_MMHg, nat=nat)
  mayor.iso.mw <- merge(mw3, mayor.isotope, by = 'isotope', all.x = T)%>%
    tibble::column_to_rownames(var = 'isotope')
  colnames(mayor.iso.mw) <- c('M.w', 'mayor.isotope')
  #made the first step of the calculations, the deconvulation matrix.
  matrix.deconvolution.results <- deconvulation_matrix(peaks1, ratios_inv, samples = samples, injections = injections)
  #calculation of the concentration in the extract or in the water if its liquid.
  liquid <- five_isotopes_2(decon = matrix.deconvolution.results, masses,
                            quant = quant_MMHg, samples = samples, injections = injections,
                            mayor.iso.mw = mayor.iso.mw)
  concentration <- data.frame("sample"=liquid$sample,
                              "injection"=liquid$injection,
                              "isotope_MMHg"=liquid$variable,
                              "concentration_MMHg"=liquid$conc_ext)
  p <- tidyr::pivot_wider(data=concentration, names_from = isotope_MMHg, values_from = concentration_MMHg)


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
                                    "isotope_MMHg"=solid$variable,
                                    "concentration_MMHg"=solid$conc_sed)

        p <- tidyr::pivot_wider(data=concentration, names_from = isotope_MMHg, values_from = concentration_MMHg)
      }
    }
  }

  return(p)

}

