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

#' Create a workflow object of the model and the transformations
#'
#' @param model_specs the model specs
#' @param recipe_specs the recipe specs
#'
#' @return a workflow objects
#'
create_model_workflow <- function(model_specs, recipe_specs) {
  workflows::workflow() %>%
    workflows::add_model(model_specs) %>%
    workflows::add_recipe(recipe_specs)
}

#' tidyi model results output
#'
#' @param workflow_fitted_model the model workflow object that has been fitted
#'
#' @return data frame
#'
tidy_model_output <- function(workflow_fitted_model) {
  workflow_fitted_model %>%
    workflows::extract_fit_parsnip() %>%
    broom::tidy(exponentiate = TRUE)
}

#' converte long format to wide split by metabolite
#'
#' @param df long format df with metabolite column
#'
#' @return list of wide form df for each metabolite
#'
split_by_metabolite <- function(df) {
  df %>%
    column_values_to_snake_case(metabolite) %>%
    dplyr::group_split(metabolite) %>%
    purrr::map(metabolites_to_wider)
}

#' Generated results for the glm model for each metabolite
#'
#' @param df long format with metabolite column
#'
#' @return df
#'
generate_model_results <- function(df) {
  create_model_workflow(
    parsnip::logistic_reg() %>%
      parsnip::set_engine("glm"),
    df %>%
      create_recipe_spec(tidyselect::starts_with("metabolite_"))
  ) %>%
    parsnip::fit(df) %>%
    tidy_model_output()
}

#' add orginal metabolite names to model results
#'
#' @param model_results df with statictic model results
#' @param df original data
#'
#' @return table
#'
add_original_metabolite_names <- function(model_results, df) {
  df %>%
    dplyr::select(metabolite) %>%
    dplyr::mutate(term = metabolite) %>%
    column_values_to_snake_case(term) %>%
    dplyr::mutate(term = stringr::str_c("metabolite_", term)) %>%
    dplyr::distinct(term, metabolite) %>%
    dplyr::right_join(model_results, by = "term")
}

#' claculate etsimate and make nice
#'
#' @param df with lipidomics data
#'
#' @return
#' @export
#'
#' @examples
calculate_estimates <- function(df){
    df %>%
        split_by_metabolite() %>%
        purrr::map(generate_model_results) %>%
        purrr::list_rbind() %>%
        dplyr::filter(stringr::str_detect(term, "metabolite_")) %>%
        add_original_metabolite_names(df)
}
