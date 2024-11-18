M.w = function(z){
  abun.1 <- tibble::column_to_rownames(z, var = 'isotope')
  abun.2 <- abun.1[, colSums(abun.1 != 0) > 0]
  numcols1 <- ncol(abun.1)
  numcols2 <- ncol(abun.2)

  abun_frac <- as.matrix(abun.1[1:7,1:numcols1]/100)

  molecular.masses <- c(195.9658144,
                        197.9667523,
                        198.9682623,
                        199.9683093,
                        200.9702853,
                        201.9706253,
                        203.9734753)
  mol.mass<- t(as.matrix(molecular.masses))

  M.w <- t(mol.mass%*%abun_frac)

  return(M.w)
}
