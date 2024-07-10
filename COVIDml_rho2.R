
setwd('E:/5. Research - Common/10. Modelling/Teams/Ben/02. Models')
source('lag_maker.R')



covidnormal.lik <- function(theta, data, lagsize, index){
  
  # Setting the indexing and the data matrices
  
  # Sorting dates and the matrix of endogenous variables
  date <- data[,1]
  data <- data[,2:ncol(data)]
  
  # Number of endogenous variables
  nvar <- ncol(data)
  
  # Number of observations 
  nobs <- nrow(data)
  
  # Need to index the COVID observations
  covid_dates <- date[index]
  
  # Setting lag length
  #lagsize <- vars::VARselect(data)$selection[3]
  
  # Making lags for endogenous variables
  lagmatrix <- data%>%
    as.data.frame%>%
    lag_maker(everything(), 1:lagsize)%>%
    dplyr::select(contains("lag"))%>%
    dplyr::select(matches(paste0(c(1:lagsize),
                                 '$',
                                 collapse = '|')))%>%
    dplyr::slice((lagsize + 1):n())%>%
    dplyr::relocate(ends_with(as.character(1:lagsize)))%>%
    dplyr::mutate(constant = 1, .after = last_col())%>%
    dplyr::select(constant, everything())
  
  
  
  
  # Setting the model parameters
  # COVID hyper-parameters
  S = theta[1:length(index)]
  rho <- theta[(length(index) + 1)]
  
  # Slope coefficient parameters (vectorised)
  beta <- theta[(length(index) + 2):(( nvar*((lagsize*nvar) + 1) ) + (length(index) + 1))]
  
  # Error variances of the endogenous variables
  sigma <- theta[(( nvar*((lagsize*nvar) + 1) ) + (length(index) + 2)):(( nvar*((lagsize*nvar) + 1) ) + (length(index) + 1) + nvar)]
  
  
  
  
  
  # Log-likelihood components
  # 1.
  loghyper <- -(nvar)*( sum(log(S)) )
  
  for ( i in (which(as.Date(date) == covid_dates[length(index)]) + 1):nobs ) {
    
    loghyper <- loghyper + log( 1 + (S[length(index)] - 1)*rho^(i - which(as.Date(date) == covid_dates[length(index)])) )
    
  }
  
  
  
  # 2.
  logerror <- -((nobs-lagsize)/2)*(log((prod(sigma))))
  
  
  
  # Starting to construct the sum of squared errors component in the likelihood function
  logsumsq <- 0
  
  for ( i in 1:nrow(lagmatrix) ) {
    
    if ( which(as.Date(date) == date[(lagsize + i)]) < which(as.Date(date) == covid_dates[1]) ) {
      
      residual <- as.numeric(data[(lagsize + i),]) - ( kronecker(diag(nvar), t(as.numeric(as.matrix(lagmatrix[i, ])))) %*% beta )
      
      wss <- t(residual) %*% solve(diag(sigma)) %*% residual
      
    } else {}
    
    logsumsq <- logsumsq + wss
    
  }
  
  
  for ( i in which(date %in% covid_dates) ) {
    
    residual <- as.numeric(data[(i),]) - ( kronecker(diag(nvar), t(as.numeric(as.matrix(lagmatrix[(i - lagsize), ])))) %*% beta )
    
    wss <- t(residual) %*% solve( ( theta[(i - which(date %in% covid_dates)[1] + 1)]*diag(sigma) ) ) %*% residual
    
    
    logsumsq <- logsumsq + wss
    
  }
  
  
  for ( i in (which(as.Date(date) == covid_dates[length(index)]) + 1):nobs ) {
    
    residual <- as.numeric(data[(i),]) - ( kronecker(diag(nvar), t(as.numeric(as.matrix(lagmatrix[(i - lagsize), ])))) %*% beta )
    
    wss <- t(residual) %*% solve( ( ( 1 + (S[length(index)] - 1)*rho^(i - which(as.Date(date) == covid_dates[length(index)])) )*diag(sigma) ) ) %*% residual
    
    
    logsumsq <- logsumsq + wss
    
  }
  
  # 3.
  logsumsq <- -(0.5)*( logsumsq )
  
  
  # The log-likelihood function
  logl <- loghyper + logerror + logsumsq
  
  return(-logl)
  
}

