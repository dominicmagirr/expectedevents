###############################################
## single arm
###############################################


expected_events_single_arm_mu = function(mu,
                                         prop_withdrawn,
                                         dco = 28,
                                         recruitment,
                                         model){
  
  n = recruitment$n
  
  expected_events_single_arm(dco = dco,
                             recruitment = recruitment,
                             model = model,
                             mu = mu,
                             dropouts = TRUE,
                             total_only = TRUE) / n - prop_withdrawn
  
}



#' Convert a proportion of dropouts to a dropout rate.
#' 
#' \code{prop_to_mu_single_arm} converts a proportion of patients who have withdrawn at data cut-off \code{dco}
#'   to the corresponding exponential dropout rate. This takes into acount the competing 
#'   risk of the event-of-interest.
#' @param{dco} Time of data cut-off. 
#' @param{recruitment} List of recruitment information. 
#'   Containing \enumerate{
#'                 \item Sample size, \code{n} 
#'                 \item Recruitment period, \code{r_period}
#'                 \item Recruitment parameter for power model, \code{k} 
#'               }
#' @param{model} The piecewise hazard model.
#'   A list containing the \code{change_points} and \code{lambdas}.
#' @return The dropout rate \code{mu} that would lead to \code{prop_withdrawn} at \code{dco}.
#' @export


prop_to_mu_single_arm = function(prop_withdrawn = 0.2,
                                 dco = 28,
                                 recruitment,
                                 model){
  
  uniroot(expected_events_single_arm_mu, 
          c(0.0001, 1000),
          prop_withdrawn = prop_withdrawn,
          dco = dco,
          recruitment = recruitment,
          model = model)$root
  
  
}


####################################################################
## two arm
####################################################################



expected_events_two_arm_mu = function(mu,
                                      prop_withdrawn,
                                      dco = 28,
                                      recruitment,
                                      model){
  
  n_0 = recruitment$n_0
  n_1 = recruitment$n_1
  
  expected_events_two_arm(dco = dco,
                          recruitment = recruitment,
                          model = model,
                          mu = mu,
                          dropouts = TRUE,
                          total_only = TRUE)["total_events"] / (n_0 + n_1) - prop_withdrawn
  
}



#' Convert a proportion of dropouts to a dropout rate.
#' 
#' \code{prop_to_mu_two_arm} converts a proportion of patients who have withdrawn at data cut-off \code{dco}
#'   to the corresponding exponential dropout rate. This takes into acount the competing 
#'   risk of the event-of-interest.
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
#' @return The dropout rate \code{mu} that would lead to \code{prop_withdrawn} at \code{dco}.
#' @export


prop_to_mu_two_arm = function(prop_withdrawn = 0.2,
                              dco = 28,
                              recruitment,
                              model){
  
  uniroot(expected_events_two_arm_mu, 
          c(0.0001, 1000),
          prop_withdrawn = prop_withdrawn,
          dco = dco,
          recruitment = recruitment,
          model = model)$root
  
  
}



