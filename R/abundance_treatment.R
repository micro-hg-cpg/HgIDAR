abun_treat=function(z) {
    # Convert first column to row names
    abun.1 <- tibble::column_to_rownames(z, var = 'isotope')

    #Remove the isotopes that we wont use to calculate the inverse of the ratios.
    abun.2 <- as.data.frame(t(abun.1))%>%
      tibble::rownames_to_column()%>%
      tidyr::separate(rowname, c('isotope', 'number'), sep = '_')
    numrows <- nrow(abun.2)
    if (numrows == 5 ) {
      abun.3 <- abun.2%>%
        dplyr::select(1:2 |contains(paste(abun.2[1,2]))|contains(paste(abun.2[2,2]))|
                        contains(paste(abun.2[3,2]))|contains(paste(abun.2[4,2]))|
                        contains(paste(abun.2[5,2])))%>%
        tidyr::unite('isotope', c(isotope, number), sep = '_')

      abun.4 <- as.data.frame(t(abun.3))
      colnames(abun.4) <- abun.4[1,]
      abun.4 <- abun.4[-1,]
      abun.5 <- as.data.frame(lapply(abun.4, function(x) as.numeric(as.character(x))))

      # Construct a 5x5 matrix for ratios_matrix when numcols1 != numcols2
      ratios_matrix <- c(
        abun.5[,1]/abun.5[1,1],
        abun.5[,2]/abun.5[2,2],
        abun.5[,3]/abun.5[3,3],
        abun.5[,4]/abun.5[4,4],
        abun.5[,5]/abun.5[5,5]
      )
      ratios_matrix <- matrix(ratios_matrix, nrow = 5, ncol = 5, byrow = TRUE)
      ratios_matrix <- t(ratios_matrix)
    } else {
      abun.3 <- abun.2%>%
        dplyr::select(1:2 |contains(paste(abun.2[1,2]))|contains(paste(abun.2[2,2]))|
                        contains(paste(abun.2[3,2]))|contains(paste(abun.2[4,2])))%>%
        tidyr::unite('isotope', c(isotope, number), sep = '_')

      abun.4 <- as.data.frame(t(abun.3))
      colnames(abun.4) <- abun.4[1,]
      abun.4 <- abun.4[-1,]
      abun.5 <- as.data.frame(lapply(abun.4, function(x) as.numeric(as.character(x))))

      # Construct a 4x4 matrix for ratios_matrix when numcols1 != numcols2
      ratios_matrix <- c(
        abun.5[,1]/abun.5[1,1],
        abun.5[,2]/abun.5[2,2],
        abun.5[,3]/abun.5[3,3],
        abun.5[,4]/abun.5[4,4]
      )
      ratios_matrix <- matrix(ratios_matrix, nrow = 4, ncol = 4, byrow = TRUE)

    }

    # Check if ratios_matrix is square and contains no NaN or Inf values before solving
    if (nrow(ratios_matrix) != ncol(ratios_matrix)) {
      stop("ratios_matrix must be square.")
    }
    if (any(is.nan(ratios_matrix)) || any(is.infinite(ratios_matrix))) {
      stop("ratios_matrix contains NaN or Inf values.")
    }

    # Calculate inverse of the ratios_matrix
    ratios_inv <- tryCatch({
      solve(ratios_matrix)
    }, error = function(e) {
      stop("Error in solving the ratios_matrix: ", e$message)
    })

    return(as.matrix(ratios_inv))
  }
