
lag_maker <- function(df, vars, lags) {
  
  map_lag <- lags %>% purrr::map( ~ purrr::partial(lag, n = .x) )
  
  return(
    df %>% dplyr::mutate(across({{ vars }}, .fns = map_lag, .names = "{.col}_lag{.fn}"))
  )
  
}