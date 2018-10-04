#' Time-at-risk intervals.
#' 
#' \code{get_t_pieces} returns the time at risk between intervals \code{ts} up to time \code{t}.
#' @param{t} Total time.
#' @param{ts} The breakpoints for intervals, starting at 0, ending at Inf.
#' @return A vector of length one less than \code{ts}, containing time at risk between consecutive breakpoints in \code{ts}.

get_t_pieces = function(t, ts){
  
  if (length(t) > 1) stop("t must be length 1")
  
  pmax(0, pmin(t - ts[-length(ts)], diff(ts)))
} 

surv_pieces_simple = function(t, change_points, lambdas){
  
  if (length(t) > 1) stop("t must be length 1")
  if (length(change_points) != length(lambdas) - 1) stop("require one event rate per time period")
  ts = c(0, change_points, Inf) 
  
  exp(-sum(get_t_pieces(t, ts) * lambdas))
  
}

#' Probability of event-free at time t, under piecewise constant hazard.
#' 
#' \code{surv_pieces} returns the probability of being event-free at time \code{t} assuming piecewise constant hazard.
#' @param{t} time
#' @param{change_points} The breakpoints for the piecewise hazard function.
#' @param{lambdas} A vector of length one more than \code{change_points}. 
#'   The hazard between consecutive breakpoints in \code{c(0, change_points, Inf)}.
#' @param{mu} The dropout hazard (a constant). Assumed independent to event hazard. 
#' @param{dropouts} Is dropout the event of interest? Default is FALSE.
#' @return The probability of being free of the event-of-interest at time \code{t}.
#' @export


surv_pieces = function(t, change_points, lambdas, mu = 0, dropouts = FALSE){
  
  if (length(t) > 1) stop("t must be length 1")
  if (length(change_points) != length(lambdas) - 1) stop("require one event rate per time period")
  ts = c(0, change_points, Inf) 
  
  both_surv = surv_pieces_simple(t, 
                                 change_points = change_points,
                                 lambdas = lambdas + mu)
  
  ## term 1: pr(T > C)
  
  term_1 = mu / (lambdas + mu)
  
  ## or, if dropout is the event of interest, pr(C > T) 
  
  if (dropouts){
    term_1 = lambdas / (lambdas + mu) 
  }
  
  ## term 2: pr(no event by t_{i-1})
  
  term_2 = unlist(purrr::map(ts[-length(ts)], 
                             surv_pieces_simple, 
                             change_points = change_points, 
                             lambdas = lambdas + mu))
  
  ## term 3: pr(event by t_i | no event by t_{i-1})
  
  term_3 = 1 - exp(-(mu + lambdas) * get_t_pieces(t, ts))
  
  
  other_happens_first = sum(term_1 * term_2 * term_3)
  
  return(both_surv + other_happens_first)
  
}