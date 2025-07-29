mayor.isotope.fn=function(z,incub_IHg,
                    incub_MMHg=NULL, quant_IHg, quant_MMHg, nat){
  if (is.null(incub_MMHg)) {
    IHg_incub <- z%>%
      dplyr::select(1,contains(paste(incub_IHg)))
    colnames(IHg_incub) <- c('isotope','abun')
    IHg_incub <- IHg_incub%>%
      dplyr::filter(isotope==paste(incub_IHg))


    IHg_quant <- z%>%
      dplyr::select(1,contains(paste(quant_IHg)))
    colnames(IHg_quant) <- c('isotope','abun')
    IHg_quant <- IHg_quant%>%
      dplyr::filter(isotope==paste(quant_IHg))


    natural <- z%>%
      dplyr::select(1,contains(paste(nat)))
    colnames(natural) <- c('isotope','abun')
    natural <- natural%>%
      dplyr::filter(isotope==paste(nat))


    MMHg_quant <- z%>%
      dplyr::select(1,contains(paste(quant_MMHg)))
    colnames(MMHg_quant) <- c('isotope','abun')
    MMHg_quant<-MMHg_quant%>%
      dplyr::filter(isotope==paste(quant_MMHg))

    mi1 <- rbind(IHg_incub, IHg_quant,
                 natural, MMHg_quant)

    rm(IHg_incub, IHg_quant,
       natural, MMHg_quant)

    mi2<- mi1%>%
      dplyr::arrange(isotope)
  } else {
    IHg_incub <- z%>%
      dplyr::select(1,contains(paste(incub_IHg)))
    colnames(IHg_incub) <- c('isotope','abun')
    IHg_incub <- IHg_incub%>%
      dplyr::filter(isotope==paste(incub_IHg))

    MMHg_incub <- z%>%
      dplyr::select(1,contains(paste(incub_MMHg)))
    colnames(MMHg_incub) <- c('isotope','abun')
    MMHg_incub <- MMHg_incub%>%
      dplyr::filter(isotope==paste(incub_MMHg))


    IHg_quant <- z%>%
      dplyr::select(1,contains(paste(quant_IHg)))
    colnames(IHg_quant) <- c('isotope','abun')
    IHg_quant <- IHg_quant%>%
      dplyr::filter(isotope==paste(quant_IHg))


    natural <- z%>%
      dplyr::select(1,contains(paste(nat)))
    colnames(natural) <- c('isotope','abun')
    natural <- natural%>%
      dplyr::filter(isotope==paste(nat))


    MMHg_quant <- z%>%
      dplyr::select(1,contains(paste(quant_MMHg)))
    colnames(MMHg_quant) <- c('isotope','abun')
    MMHg_quant<-MMHg_quant%>%
      dplyr::filter(isotope==paste(quant_MMHg))

    mi1 <- rbind(IHg_incub, MMHg_incub, IHg_quant,
                 natural, MMHg_quant)

    rm(IHg_incub, MMHg_incub, IHg_quant,
       natural, MMHg_quant)

    mi2<- mi1%>%
      dplyr::arrange(isotope)

  }

  return(mi2)

}


