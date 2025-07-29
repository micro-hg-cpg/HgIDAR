summarise_isotopes <- function(data) {
  final_output=data %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(dplyr::across(tidyselect::contains("m"),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE)),
                     .names = "{.fn}_{.col}")) %>%
    dplyr::ungroup()


  return(final_output)
}
