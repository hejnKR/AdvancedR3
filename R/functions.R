#' Descriptive stats
#'
#' @param df data frame containing a metabolite column.
#'
#' @return a data.frame/tibble summarizing the mean and sd for each observation in the metabolite column.
descriptive_stats <- function(df) {
  df %>%
    dplyr::group_by(metabolite) %>%
    dplyr::summarise(dplyr::across(value, list(mean = mean, sd = sd))) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, digits = 1)))
}
