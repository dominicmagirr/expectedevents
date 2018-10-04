########################################
## conditional probability of an event
## occuring before patient time 'min_f' and
## before calendar time 'r_period + f_period',
## given that a patient is recruited at 'r_time'.
## Note: 'r_time' is a vector here

cond_p_event = function(r_time,
                        min_f,
                        r_period,
                        f_period,
                        change_points,
                        lambdas,
                        mu,
                        dropouts){
  
  
  f_time = pmin(min_f, r_period + f_period - r_time)
  
  event_prob = numeric(length(f_time))
  
  for (i in seq_along(event_prob)){
    
    event_prob[i] = 1 - surv_pieces(f_time[i], change_points, lambdas, mu, dropouts)
    
  }
  
  event_prob
  
}

########################################
## probability density of recruitment time

p_r = function(r_time, r_period, k) k * (r_time / r_period) ^ (k - 1) / r_period

########################################
## cond_p_event times p_r

p_event_r = function(r_time,
                     min_f,
                     r_period,
                     f_period,
                     change_points,
                     lambdas,
                     mu,
                     dropouts,
                     k){
  cond_p_event(r_time,
               min_f,
               r_period,
               f_period,
               change_points,
               lambdas,
               mu,
               dropouts) * p_r(r_time, r_period, k)
}
########################################
## probability of an event occuring 
## prior to patient-time c(change_points, Inf) AND
## before calendar time 'r_period + f_period'

p_event = function(r_period,
                   f_period,
                   change_points,
                   lambdas,
                   mu,
                   dropouts,
                   k){
  
  min_fs = c(change_points, Inf)
  
  integrals = numeric(length(min_fs))
  
  for (i in seq_along(integrals)){
    
    integrals[i] = integrate(p_event_r, 
                             lower = 0,
                             upper = r_period,
                             min_f = min_fs[i],
                             r_period = r_period,
                             f_period = f_period,
                             change_points = change_points,
                             lambdas = lambdas,
                             mu = mu,
                             dropouts = dropouts,
                             k = k)$value 
  }
  
  integrals
}


#' Expected number of events
#' 
#' \code{expected_events_single_arm} returns the expected number of events at a given data-cut-off time.
#' @param{dco} Time of data cut-off. 
#' @param{recruitment} List of recruitment information. 
#'   Containing \enumerate{
#'                 \item Sample size, \code{n} 
#'                 \item Recruitment period, \code{r_period}
#'                 \item Recruitment parameter for power model, \code{k} 
#'               }
#' @param{model} The piecewise hazard model.
#'   A list containing the \code{change_points} and \code{lambdas}.
#' @param{mu} The dropout rate. Assumed independent to event hazard. 
#' @param{dropouts} Is dropout the event of interest? Default is FALSE.
#' @param{total_only} By default, the function returns only the expected number of events.
#'   To see more detailed output set \code{total_only = FALSE}.
#' @return The expected number of events at \code{dco}.
#' @export


expected_events_single_arm = function(dco = 28,
                                      recruitment,
                                      model,
                                      mu = 0,
                                      dropouts = FALSE,
                                      total_only = TRUE){
  
  n = recruitment$n
  r_period = recruitment$r_period
  k = recruitment$k
  
  change_points = model$change_points
  lambdas = model$lambdas
    
  if (length(dco) > 1) stop("dco must be length 1")
  if (length(n) > 1) stop("n must be length 1")
  
  # assume dco is greater than r_period.
  # what if dco < r_period
  
  if (dco < r_period){
    
    n = n * (dco / r_period) ^ k
    r_period = dco
    
  }
  
  f_period = dco - r_period
  
  expected_events = n * p_event(r_period,
                                f_period,
                                change_points,
                                lambdas,
                                mu,
                                dropouts,
                                k)
  
  total_events = expected_events[length(expected_events)]
  
  if (total_only) return(total_events)
  
  
  events_per_period = diff(c(0, expected_events))
  
  prop_events = events_per_period / total_events
  
  list(total_events = total_events,
       events_per_period = cbind(c(change_points, Inf), events_per_period),
       prop_events = prop_events,
       n_recruited = n)
  
}



#' Expected number of events
#' 
#' \code{expected_events_two_arm} returns the expected number of events at a given data-cut-off time.
#' @param{dco} Time of data cut-off. 
#' @param{recruitment} List of recruitment information. 
#'   Containing \enumerate{
#'                 \item Sample size on control, \code{n_0} 
#'                 \item Sample size on experimental, \code{n_1} 
#'                 \item Recruitment period, \code{r_period}
#'                 \item Recruitment parameter for power model, \code{k} 
#'               }
#' @param{model} The piecewise hazard model.
#'   A list containing the \code{change_points} and \code{lambdas}.
#' @param{mu} The dropout rate. Assumed independent to event hazard. 
#' @param{dropouts} Is dropout the event of interest? Default is FALSE.
#' @param{total_only} By default, the function returns only the expected number of events.
#'   To see more detailed output set \code{total_only = FALSE}.
#' @return The expected number of events at \code{dco}.
#' @export

expected_events_two_arm = function(dco = 28,
                                   recruitment,
                                   model,
                                   mu = 0,
                                   dropouts = FALSE,
                                   total_only = TRUE){
  
  n_0 = recruitment$n_0
  n_1 = recruitment$n_1
  r_period = recruitment$r_period
  k = recruitment$k
  
  change_points = model$change_points
  lambdas_0 = model$lambdas_0
  lambdas_1 = model$lambdas_1
  
  if (length(dco) > 1) stop("dco must be length 1")
  if (length(n_0) > 1) stop("n_0 must be length 1")
  if (length(n_1) > 1) stop("n_1 must be length 1")
  
  # assume dco is greater than r_period.
  # what if dco < r_period
  
  if (dco < r_period){
    
    n_0 = n_0 * (dco / r_period) ^ k
    n_1 = n_1 * (dco / r_period) ^ k
    r_period = dco
    
  }
  
  f_period = dco - r_period
  
  expected_events_0 = n_0 * p_event(r_period,
                                    f_period,
                                    change_points,
                                    lambdas_0,
                                    mu,
                                    dropouts,
                                    k)
  
  expected_events_1 = n_1 * p_event(r_period,
                                    f_period,
                                    change_points,
                                    lambdas_1,
                                    mu,
                                    dropouts,
                                    k)
  
  
  total_events_0 = expected_events_0[length(expected_events_0)]
  total_events_1 = expected_events_1[length(expected_events_1)]
  total_events = total_events_0 + total_events_1
  
  if (total_only) return(c(total_events_0 = total_events_0,
                           total_events_1 = total_events_1,
                           total_events = total_events))
  
  
  events_per_period = cbind(events_0 = diff(c(0, expected_events_0)), 
                            events_1 = diff(c(0, expected_events_1)))
  
  prop_events = rowSums(events_per_period) / total_events
  
  hrs = lambdas_1 / lambdas_0
  
  av_hr = exp(sum(log(hrs) * prop_events))
  
  list(total_events_0 = total_events_0,
       total_events_1 = total_events_1,
       total_events = total_events,
       events_per_period = cbind(c(change_points, Inf), events_per_period),
       prop_events = prop_events,
       hrs = hrs,
       av_hr = av_hr,
       n_recruited_0 = n_0,
       n_recruited_1 = n_1,
       n_recruited = n_0 + n_1)
  
}



