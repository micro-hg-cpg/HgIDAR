five_isotopes_2=function(decon, masses,  quant, samples, injections, mayor.iso.mw ){
  Samples <- c(1:samples)
  Injection <- c(1:injections)
  third.results <- data.frame()
  fourth.results <- data.frame()
  fifth.results <- data.frame()
  sixth.results <- data.frame()
  concentration_extract <- data.frame()


  for (i in Samples) {
    for (j in Injection) {
      b <- decon %>%
        dplyr::filter(sample == i & injection == j)
      b <- b%>%
        dplyr::mutate("MD.A"=V1/mayor.iso.mw$mayor.isotope)
      third.results<- rbind(third.results, b)
      c <- third.results %>%
        dplyr::filter(sample == i & injection == j)

      c2 <- c[grepl(paste(quant), c$variable),]
      c<- c%>%dplyr::mutate("fac.MD.A"= MD.A/c2[1,6])

      fourth.results<- rbind(fourth.results, c)
      d <- masses%>%
        dplyr::filter(sample == i)
      e <- fourth.results%>%
        dplyr::filter(sample==i)
      e <- e %>%
        dplyr::mutate("facxadd"=fac.MD.A*d$mol_added)
      fifth.results <- rbind(fifth.results, e)
      f <- fifth.results %>%
        dplyr::filter(sample == i & injection == j)
      f <- f%>%
        cbind(mayor.iso.mw$M.w)%>%
        dplyr::mutate("facxaddxmw"=facxadd*mayor.iso.mw$M.w)
      sixth.results<- rbind(sixth.results, f)
      ext_mass <- masses%>%
        dplyr::filter(sample ==i)
      g <- sixth.results %>%
        dplyr::filter(sample == i & injection == j)
      g <- g%>%
        dplyr::mutate("conc_ext"=facxaddxmw/ext_mass$extract)
      concentration_extract<- rbind(concentration_extract, g)
    }
  }
  return(concentration_extract)
}
