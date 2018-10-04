
#' Expected time of dco
#' 
#' \code{expected_dco_single_arm} returns the expected number of events at a given data-cut-off time.
#' @param{total_events} Number of events
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
#' @return The expected calendar time of reaching \code{total_events}.
#' @export




expected_dco_single_arm = function(total_events,
                                   recruitment,
                                   model,
                                   mu,
                                   dropouts){
  
  if (length(total_events) > 1) stop("total_events must be length 1")
  
  find_dco = function(x){
    
    
    expected_events_single_arm(dco = x,
                               recruitment,
                               model,
                               mu,
                               dropouts,
                               total_only = TRUE) - total_events
    
  }
  
  
  uniroot(find_dco, c(0.001, 1000))$root 
  
}


#' Expected time of dco
#' 
#' \code{expected_dco_two_arm} returns the expected calendar time of data-cut-off for a given number of events.
#' @param{total_events} Number of events.
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
#' @return The expected calendar time of reaching \code{total_events}.
#' @export



expected_dco_two_arm = function(total_events,
                                recruitment,
                                model,
                                mu,
                                dropouts){
  
  if (length(total_events) > 1) stop("total_events must be length 1")
  
  find_dco = function(x){
    
    
    expected_events_two_arm(dco = x,
                            recruitment,
                            model,
                            mu,
                            dropouts,
                            total_only = TRUE)["total_events"] - total_events
    
  }

  
  uniroot(find_dco, c(0.001, 1000))$root 
  
}