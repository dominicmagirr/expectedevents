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



