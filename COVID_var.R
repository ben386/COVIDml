

covid_var <- function(data, lags, var_type = "const", covidindex) {
  
  simplevar <- vars::VAR(data %>% dplyr::select(where(is.numeric)) %>% as.matrix,
                         p = lags,
                         type = var_type)
  
  s0 <- abs(data[covidindex, 2:ncol(data)] %>% as.matrix) %*% (1/matrixStats::colSds(data[,2:ncol(data)] %>% as.matrix, na.rm = TRUE)) * (1/(ncol(data)-1))
  
  beta0 <- purrr::map(coefficients(simplevar), ~ { .x %>% as_tibble() %>% dplyr::slice( c(n(), (1:(n() - 1))) ) %>% dplyr::select(1) })
  
  sigma0 <- diag(t(residuals(simplevar)) %*% residuals(simplevar)) * (1/nrow(data))
  
  theta <- c(s0, do.call(rbind.data.frame, beta0) %>% as.matrix, sigma0)
  
  
  gls <- optim(theta, covidnormal.lik, method = "Nelder-Mead", hessian = F, data = data, lagsize = lags, index = covidindex)
  
  
  
  df <- data
  
  for ( i in 1:length(covidindex) ) {
    
    df[covidindex[i], 2:ncol(data)] <- df[covidindex[i], 2:ncol(data)]/gls$par[i]
    
  }
  
  
  
  # not the best way to retrieve model results by re-weight the data using the hyper-parameter estimates
  # will have to manually compute the residuals, sigma, test statistics etc. 
  var_weighted <- vars::VAR(df %>% dplyr::select(where(is.numeric)) %>% as.matrix, 
                            p = lags, 
                            type = var_type)
  
  resid <- residuals(simplevar) 
  
  for ( i in 1:length(covidindex) ) {
    
    resid[covidindex[i] - lags, 2:ncol(resid)] <- resid[covidindex[i] - lags, 2:ncol(resid)] * gls$par[i]
    
  }
  
  
  
  
  result <- list("Beta" = broom::tidy(var_weighted), # standard test statistics are invalid here
                 "In-Sample Fit" = purrr::map(var_weighted$varresult, ~ { broom::glance(.x) }), # not technically correct
                 "Unweighted residuals" = resid,
                 "Weighted VAR" = var_weighted
                 )
  
  return(result)
  
  
  
}
