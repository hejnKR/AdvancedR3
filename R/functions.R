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

#' Distribution plot of metabolites
#'
#' @param df with value and metabolite column
#'
#' @return histograms for for metabolite values
#'
plot_distribution <- function(df) {
  df %>% ggplot2::ggplot(aes(x = value)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(vars(metabolite), scales = "free")
}

#' column values to snake case
#'
#' @param df cols to make to snake case in df
#'
#' @return df
#'
column_values_to_snake_case <- function(df, cols) {
  df %>%
    dplyr::mutate(dplyr::across({{ cols }}, snakecase::to_snake_case))
}

#' pivot metabolite to wider
#'
#' @param df containing columns metabolite, value, and mean
#'
#' @return df wide format
#'
metabolites_to_wider <- function(df) {
  df %>%
    tidyr::pivot_wider(
      names_from = metabolite,
      values_from = value,
      values_fn = mean,
      names_prefix = "metabolite_"
    )
}

#' Create recipe spec to pre-process data
#'
#' @param df wide df
#' @param metabolite_variable column of metabolite of intrest
#'
#' @return returns recipe with specifications
#'
create_recipe_spec <- function(df, metabolite_variable) {
  recipes::recipe(df) %>%
    recipes::update_role({{ metabolite_variable }}, age, gender, new_role = "predictor") %>%
    recipes::update_role(class, new_role = "outcome") %>%
    recipes::step_normalize(tidyselect::starts_with("metabolite_"))
}
