devtools::load_all()
##############################################################
##############################################################
get_t_pieces(t = 5, 
             ts = c(0, 0.001, 200))



surv_pieces(t = 5, 
            change_points = 10,
            lambdas = c(log(2) / 15, 
                        log(2) / 15), 
            mu = 0.12, 
            dropouts = FALSE)

################################################
#devtools::use_vignette("dropout_explainer")
################################################


## bring in recruitment rate...

recruitment_1 = list(n = 100, r_period = 12, k = 1)
model_1 = list(change_points = 4, lambdas = c(log(2) / 15, log(2) / 15))


expected_events_single_arm(dco = 30,
                           recruitment = recruitment_1,
                           model = model_1,
                           mu = 0,
                           dropouts = FALSE,
                           total_only = FALSE)


dcos = seq(0.001, 30, length.out = 100)

events_dco = purrr::map(dcos,
                        expected_events_single_arm,
                        recruitment = recruitment_1,
                        model = model_1,
                        mu = 0.01,
                        total_only = FALSE)


events_dco_n = unlist(purrr::map(events_dco, function(x) x$n_recruited))
events_dco_e = unlist(purrr::map(events_dco, function(x) x$total_events))

dropouts_dco = purrr::map(dcos,
                          expected_events_single_arm,
                          recruitment = recruitment_1,
                          model = model_1,
                          mu = 0.01,
                          dropouts = TRUE)


plot(dcos, events_dco_n, type = 'l', col = 1, xlab = "time", ylab = "count")
points(dcos, events_dco_e, type = 'l', col = 2)
points(dcos, unlist(dropouts_dco), type = 'l', col = 3)


legend("topleft", c("recruited", "events", "withdrawn"), lty = c(1,1,1), col = c(1,2,3))

## suppose I wanted 20% dropout at 30 months

prop_to_mu_single_arm(prop_withdrawn = 0.2,
                      dco = 30,
                      recruitment = recruitment_1,
                      model = model_1)

##################################
## two arm

## I would need to make mu 0.0161

recruitment_2 = list(n_0 = 100, n_1 = 100, r_period = 12, k = 1)
model_2 = list(change_points = 4, lambdas_0 = c(log(2) / 15, log(2) / 15), lambdas_1 = c(log(2) / 15, 0.7 * log(2) / 15))


expected_events_two_arm(dco = 30,
                        recruitment = recruitment_2,
                        model = model_2,
                        mu = 0.0153,
                        dropouts = FALSE,
                        total_only = FALSE)


prop_to_mu_two_arm(prop_withdrawn = 0.2,
                   dco = 30,
                   recruitment = recruitment_2,
                   model = model_2)

events_dco = purrr::map(dcos,
                        expected_events_two_arm,
                        recruitment = recruitment_2,
                        model = model_2,
                        mu = 0.0153,
                        total_only = FALSE)


events_dco_n = unlist(purrr::map(events_dco, function(x) x$n_recruited))
events_dco_e = unlist(purrr::map(events_dco, function(x) x$total_events))

dropouts_dco = purrr::map(dcos,
                          expected_events_two_arm,
                          recruitment = recruitment_2,
                          model = model_2,
                          mu = 0.0153,
                          dropouts = TRUE)

dropouts_dco_e = unlist(purrr::map(dropouts_dco, function(x) x["total_events"]))

plot(dcos, events_dco_n, type = 'l', col = 1, xlab = "time", ylab = "count")
points(dcos, events_dco_e, type = 'l', col = 2)
points(dcos, dropouts_dco_e, type = 'l', col = 3)


legend("topleft", c("recruited", "events", "withdrawn"), lty = c(1,1,1), col = c(1,2,3))

devtools::use_vignette("expectedevents")
devtools::use_package("purrr")
