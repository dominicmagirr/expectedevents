#' Power under non-constant hazard ratio.
#' 
#' \code{get_power} returns the power of a group sequential test assuming a non-constant hazard ratio.
#' \param{events} Number of events at each analysis.
#' \param{av_hrs} Average hazard ratio at each analysis.
#' \param{upper_bound} The upper stopping boundary for the (standardized) test statistic.
#' \param{lower_bound} The upper stopping boundary for the (standardized) test statistic. Default is NULL.
#' \param{recruitment} List of recruitment information. 
#'   Containing \enumerate{
#'                 \item Sample size on control, \code{n_0} 
#'                 \item Sample size on experimental, \code{n_1} 
#'                 \item Recruitment period, \code{r_period}
#'                 \item Recruitment parameter for power model, \code{k} 
#'               }
#' \return The cumulative power at each analysis.
#' @export


get_power = function(events,
                     av_hrs,
                     upper_bound,
                     lower_bound = NULL,
                     recruitment){  
  
  R = recruitment$n_1 / recruitment$n_0
  
  mu = -log(av_hrs) * sqrt(events * R / (1 + R) ^ 2)
  
  n_analyses = length(events)
  
  m1 = matrix(events, 
              nrow = n_analyses, 
              ncol = n_analyses)
  
  
  sigma = sqrt(pmin(m1, t(m1)) / pmax(m1, t(m1)))
  
  if (is.null(lower_bound)){
    
    power_so_far = numeric(n_analyses)
    
    for (i in 1:n_analyses){
      
      power_so_far[i] = 1 - mvtnorm::pmvnorm(lower = rep(-Inf, i),
                                             upper = upper_bound[1:i],
                                             mean = mu[1:i],
                                             sigma = sigma[1:i, 1:i])[1]
      
    }
    
    upper_stopping = c(power_so_far[1], diff(power_so_far))
    
    return(list(upper_stopping = upper_stopping,
                power_so_far = power_so_far))
    
  }
  else {
    
    upper_stopping = numeric(n_analyses)
    lower_stopping = numeric(n_analyses)
    
    upper_stopping[1] = 1 - pnorm(upper_bound[1], 
                                  mu[1], 
                                  sqrt(sigma[1,1]))
    
    lower_stopping[1] = pnorm(lower_bound[1], 
                              mu[1], 
                              sqrt(sigma[1,1]))
    
    
    for (i in 2:n_analyses){
      
      upper_stopping[i] = mvtnorm::pmvnorm(lower = c(lower_bound[1:(i-1)], upper_bound[i]),
                                           upper = c(upper_bound[1:(i-1)], Inf),
                                           mean = mu[1:i],
                                           sigma = sigma[1:i, 1:i])[1]
     
      lower_stopping[i] = mvtnorm::pmvnorm(lower = c(lower_bound[1:(i-1)], -Inf),
                                           upper = c(upper_bound[1:(i-1)], lower_bound[i]),
                                           mean = mu[1:i],
                                           sigma = sigma[1:i, 1:i])[1] 
    }
    
    power_so_far = cumsum(upper_stopping)
  }
  
  list(upper_bound = upper_bound,
       lower_bound = lower_bound,
       events = events,
       av_hrs = av_hrs,
       upper_stopping = upper_stopping,
       power_so_far = power_so_far,
       lower_stopping = lower_stopping)
  
}

