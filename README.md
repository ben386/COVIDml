Estimating a VAR via maximum likelihood estimation with COVID observations. This is following the methodology from Lenza & Primicri (2020)

The way to estimate the VAR is by reading all the scripts/functions and executing the function "covid_var()" in the "covid_var.R" script. 

The "example.R" script gives a basic overview.

NOTE: this is just basic, preliminary code for applying the ML estimation method. More work is needed to be done with properly computing residuals and the parameter estimates in the error covariance matrix. Usual test statistics are no longer valid and need to be adjusted accordingly (if we are assuming the errors are scaled exogenously because of the pandemic). 
