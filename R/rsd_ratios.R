#' Check the quality of your integrations using the relative standard deviation of your peak's ratios
#' @description
#' A short description...
#'
#' @param x data frame with your integrated peaks
#' @param form IHg, MMHg, all
#' @param samples Numeric vector. Number of samples
#' @param cutoff Numeric value. cutoff used for the plot. Default = 10%
#' @param incub_IHg IHg isotope used for the incubation
#' @param incub_MMHg MMHg isotope used for the incubation
#' @param quant_MMHg MMHg isotope used for the quantification
#' @param quant_IHg IHg isotope used for the quantification
#' @param nat isotope used as natural for the incubation

#' @importFrom dplyr %>%
#' @return a table and one plot per sample of the ones that are above the cutoff.
#' @export
#'
#' @examples
#' rsd_ratios(peaks_1,
#' form = 'MMHg',
#' sample = 24,
#' incub_IHg = 199,
#' incub_MMHg = 201,
#' quant_IHg = 198,
#' quant_MMHg = 202,
#' nat = 200)
#' rsd_ratios(peaks_2,
#' form = 'all',
#' sample = 26,
#' incub_IHg = 200,
#' incub_MMHg = 204,
#' quant_IHg = 199,
#' quant_MMHg = 201,
#' nat = 202)
#' rsd_ratios(peaks_3, form = 'MMHg', 96, 5, 199, 201, 198, 202, 200)
rsd_ratios=function(x, form = c('IHg', 'MMHg', 'all'), samples, cutoff = 10, incub_IHg,
                         incub_MMHg, quant_IHg, quant_MMHg, nat){
  form <- match.arg(form)
  peaks1 <- x[, colSums(x != 0) > 0]
  Samples <- c(1:samples)

    if (form == 'IHg') {
      peaks_IHg <- peaks1 %>%
        dplyr::select(1:2 |contains("IHg"))
      rp <- data.frame()
      Peaks_ratios <- data.frame()

      for (i in Samples) {
        rp <- peaks_IHg%>%
          dplyr::filter(sample == i)

        IHg_added_incub <- rp%>%
          dplyr::select(contains(paste(incub_IHg)))
          colnames(IHg_added_incub) <- 'IHg_added_incub'

        IHg_formed <- rp%>%
          dplyr::select(contains(paste(incub_MMHg)))
        colnames(IHg_formed) <- 'IHg_formed'

        IHg_added_quant <- rp%>%
          dplyr::select(contains(paste(quant_IHg)))
        colnames(IHg_added_quant) <- 'IHg_added_quant'

        natural <- rp%>%
          dplyr::select(contains(paste(nat)))
        colnames(natural) <- 'natural'

        natural2 <- rp%>%
          dplyr::select(contains(paste(quant_MMHg)))
        colnames(natural2) <- 'natural2'

        rp1 <- cbind(rp, IHg_added_incub, IHg_formed, IHg_added_quant,
                     natural, natural2)
        rm(IHg_added_incub, IHg_formed, IHg_added_quant,
           natural, natural2)

        rp2<- rp1%>%
          dplyr::mutate("Ratio_incub_quant" = IHg_added_incub/IHg_added_quant)%>%
          dplyr::mutate("RSD_ratio_incub_quant" = rtemisalpha::rsd(`Ratio_incub_quant`))%>%
          dplyr::mutate("Ratio_formed_quant"= IHg_formed/IHg_added_quant)%>%
          dplyr::mutate("RSD_formed_quant"= rtemisalpha::rsd( `Ratio_formed_quant`))%>%
          dplyr::mutate("Ratio_formed_nat"= IHg_formed/natural)%>%
          dplyr::mutate("RSD_formed_nat"= rtemisalpha::rsd( `Ratio_formed_nat`))%>%
          dplyr::mutate("Ratio_natural"= natural2/natural)
        Peaks_ratios <- rbind(Peaks_ratios, rp2)
      }
      Peaks_ratios <- Peaks_ratios%>%
        dplyr::select(-c(IHg_added_incub, IHg_formed, IHg_added_quant,
                  natural, natural2))

      # RSD <- Peaks_ratios %>%
      #   dplyr::filter(RSD_ratio_incub_quant > cutoff)

    }

    if (form == 'MMHg'){
      peaks_MMHg <- peaks1 %>%
        dplyr::select(1:2 |contains("MMHg"))
      rp <- data.frame()
      Peaks_ratios <- data.frame()

      for (i in Samples) {
        rp <- peaks_MMHg%>%
          dplyr::filter(sample == i)

        MMHg_formed <- rp%>%
          dplyr::select(contains(paste(incub_IHg)))
        colnames(MMHg_formed) <- 'MMHg_formed'

        MMHg_added_incub <- rp%>%
          dplyr::select(contains(paste(incub_MMHg)))
        colnames(MMHg_added_incub) <- 'MMHg_added_incub'

        natural2 <- rp%>%
          dplyr::select(contains(paste(quant_IHg)))
        colnames(natural2) <- 'natural2'

        natural <- rp%>%
          dplyr::select(contains(paste(nat)))
        colnames(natural) <- 'natural'

        MMHg_added_quant <- rp%>%
          dplyr::select(contains(paste(quant_MMHg)))
        colnames(MMHg_added_quant) <- 'MMHg_added_quant'

        rp1 <- cbind(rp, MMHg_formed, MMHg_added_incub, MMHg_added_quant,
                     natural, natural2)

        rm(MMHg_formed, MMHg_added_incub, MMHg_added_quant,natural, natural2)

        rp2<- rp1%>%
          dplyr::mutate("Ratio_incub_quant" = MMHg_added_incub/MMHg_added_quant)%>%
          dplyr::mutate("RSD_ratio_incub_quant" = rtemisalpha::rsd(`Ratio_incub_quant`))%>%
          dplyr::mutate("Ratio_formed_quant"= MMHg_formed/MMHg_added_quant)%>%
          dplyr::mutate("RSD_formed_quant"= rtemisalpha::rsd( `Ratio_formed_quant`))%>%
          dplyr::mutate("Ratio_formed_nat"= MMHg_formed/natural)%>%
          dplyr::mutate("RSD_formed_nat"= rtemisalpha::rsd( `Ratio_formed_nat`))%>%
          dplyr::mutate("Ratio_natural"= natural2/natural)
        Peaks_ratios <- rbind(Peaks_ratios, rp2)
      }
      Peaks_ratios <- Peaks_ratios%>%
        dplyr::select(-c(MMHg_formed, MMHg_added_incub, natural2,
                  natural, MMHg_added_quant))
      # RSD <- Peaks_ratios %>%
      #   dplyr::filter(RSD_ratio_incub_quant > cutoff)

    }

    if (form == 'all') {
      peaks_MMHg <- peaks1 %>%
        dplyr::select(1:2 |contains("MMHg"))
      rp.mm <- data.frame()
      Peaks_ratios.mm <- data.frame()

      for (i in Samples) {
        rp.mm <- peaks_MMHg%>%
          dplyr::filter(sample == i)
        MMHg_formed <- rp.mm%>%
          dplyr::select(contains(paste(incub_IHg)))
        colnames(MMHg_formed) <- 'MMHg_formed'

        MMHg_added_incub <- rp.mm%>%
          dplyr::select(contains(paste(incub_MMHg)))
        colnames(MMHg_added_incub) <- 'MMHg_added_incub'

        natural2 <- rp.mm%>%
          dplyr::select(contains(paste(quant_IHg)))
        colnames(natural2) <- 'natural2'

        natural <- rp.mm%>%
          dplyr::select(contains(paste(nat)))
        colnames(natural) <- 'natural'

        MMHg_added_quant <- rp.mm%>%
          dplyr::select(contains(paste(quant_MMHg)))
        colnames(MMHg_added_quant) <- 'MMHg_added_quant'

        rp.mm1 <- cbind(rp.mm, MMHg_formed, MMHg_added_incub, MMHg_added_quant,
                     natural, natural2)
        rm(MMHg_formed, MMHg_added_incub, MMHg_added_quant,natural, natural2)
        rp.mm2<- rp.mm1%>%
          dplyr::mutate("Ratio_incub_quant.mm" = MMHg_added_incub/MMHg_added_quant)%>%
          dplyr::mutate("RSD_ratio_incub_quant.mm" = rtemisalpha::rsd(`Ratio_incub_quant.mm`))%>%
          dplyr::mutate("Ratio_formed_quant.mm"= MMHg_formed/MMHg_added_quant)%>%
          dplyr::mutate("RSD_formed_quant.mm"= rtemisalpha::rsd(`Ratio_formed_quant.mm`))%>%
          dplyr::mutate("Ratio_formed_nat"= MMHg_formed/natural)%>%
          dplyr::mutate("RSD_formed_nat"= rtemisalpha::rsd( `Ratio_formed_nat`))%>%
          dplyr::mutate("Ratio_natural.mm"= natural2/natural)
        Peaks_ratios.mm <- rbind(Peaks_ratios.mm, rp.mm2)
      }
      Peaks_ratios.mm <- Peaks_ratios.mm%>%
        dplyr::select(-c(natural, MMHg_formed, MMHg_added_incub, natural2,
                         MMHg_added_quant))

      peaks_IHg <- peaks1 %>%
        dplyr::select(1:2 |contains("IHg"))
      rp.in <- data.frame()
      Peaks_ratios.in <- data.frame()


      for (i in Samples) {
        rp.in <- peaks_IHg%>%
          dplyr::filter(sample == i)

        IHg_added_incub <- rp.in%>%
          dplyr::select(contains(paste(incub_IHg)))
        colnames(IHg_added_incub) <- 'IHg_added_incub'

        IHg_formed <- rp.in%>%
          dplyr::select(contains(paste(incub_MMHg)))
        colnames(IHg_formed) <- 'IHg_formed'

        IHg_added_quant <- rp.in%>%
          dplyr::select(contains(paste(quant_IHg)))
        colnames(IHg_added_quant) <- 'IHg_added_quant'

        natural <- rp.in%>%
          dplyr::select(contains(paste(nat)))
        colnames(natural) <- 'natural'

        natural2 <- rp.in%>%
          dplyr::select(contains(paste(quant_MMHg)))
        colnames(natural2) <- 'natural2'

        rp.in1 <- cbind(rp.in, IHg_added_incub, IHg_formed, IHg_added_quant,
                     natural, natural2)
        rm(IHg_added_incub, IHg_formed, IHg_added_quant,
              natural, natural2)
        rp.in2<- rp.in1%>%
          dplyr::mutate("Ratio_incub_quant.in" = IHg_added_incub/IHg_added_quant)%>%
          dplyr::mutate("RSD_ratio_incub_quant.in" = rtemisalpha::rsd(`Ratio_incub_quant.in`))%>%
          dplyr::mutate("Ratio_formed_quant.in"= IHg_formed/IHg_added_quant)%>%
          dplyr::mutate("RSD_formed_quant.in"= rtemisalpha::rsd( `Ratio_formed_quant.in`))%>%
          dplyr::mutate("Ratio_formed_nat"= IHg_formed/natural)%>%
          dplyr::mutate("RSD_formed_nat"= rtemisalpha::rsd( `Ratio_formed_nat`))%>%
          dplyr::mutate("Ratio_natural.in"= natural2/natural)
        Peaks_ratios.in <- rbind(Peaks_ratios.in, rp.in2)
      }
      Peaks_ratios.in <- Peaks_ratios.in%>%
        dplyr::select(-c(IHg_added_incub, IHg_formed, IHg_added_quant,
                         natural, natural2))
      Peaks_ratios <- merge(Peaks_ratios.in, Peaks_ratios.mm, by= c('sample','injection'), all = T)



    }


  return(Peaks_ratios)
}


