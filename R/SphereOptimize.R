#' Converting spherical coordinates to Cartesian coordinates
#' @description The function from.Sphere convert a list of angles representing a point
#' on a unit sphere to the corresponding Cartesian coordinates.
#'
#' @param theta A list of angles. The first item should be between 0 to pi, and the following items
#' should be between 0 to 2*pi.
#'
#' @return A vector of the corresponding Cartesian coordinates.
#'
#' @examples
#' from.Sphere(c(pi/3, pi/4, pi/5))

from.Sphere = function(theta){
  s = rep(0,(length(theta)+1))
  for(k in length(theta)){
    theta[k] = .regularize.Theta(theta[k])
  }
  for(i in 1:length(s)){
    s[i] = .med.From.Sphere(theta, i)
  }
  return(s)
}



#' Converting Cartesian coordinates to spherical coordinates
#' @description The function to.Sphere convert a list of Cartesian coordinates representing a point
#' on a unit sphere to the corresponding spherical coordinates.
#'
#' @param s A list of Cartesian coordinates.
#'
#' @return A vector of the corresponding angles in spherical coordinating system.
#'
#' @examples
#' s = from.Sphere(c(pi/3, pi/4, pi/5))
#' theta = to.Sphere(s)
#' theta = round(theta, 5)
#' theta == round(c(pi/3, pi/4, pi/5), 5)

to.Sphere = function(s){
  theta = rep(0, length(s)-1)
  theta[1] = acos(s[1])
  for(i in 2:length(theta)){
    theta[i] = acos(s[i] /sqrt(sum(s[i:length(s)]^2)))
  }
  if(s[length(s)] < 0){
    theta[length(theta)] = 2*pi - theta[length(theta)]
  }
  na_id = which(is.na(theta))
  theta[na_id] = 0
  return(theta)
}

#' Conducting optimization on a unit sphere
#' @description The function SphereOptimize conducts optimization on a unit sphere.
#' If the size of neighbor near the initial value is specified, the L-BFGS-B opitmization
#' algorithm will be called. Otherwise this function searches the whole unit sphere using
#' Nelder-Mead algorithm by default. Other optimization methods are allowed.
#'
#' @param par Initial values for the parameters to be optimized over. Must be
#' in Cartesian coordinates and on a unit sphere.
#' @param fn A function to be minimized (or maximized).
#' @param neighbor Radius of neighbor to search for the optimal results. If not specified, this
#' function will search for the whole unit sphere.
#' @param ... Extra arguments that can be passed to optim().
#'
#' @return A list compose three items.
#' \itemize{
#' \item par The optimal restuls found.
#' \item value The value of fn corresponding to par.
#' \item method The optimization algorithm used.
#' }
#'
#' @examples
#' fn = function(s){
#'     return(sum(s^3))
#' }
#'
#' s = c(sqrt(2)/2,sqrt(2)/2)
#' k = SphereOptimize(s, fn, control = list(fnscale = -1))
#' k$value
#' k$par


SphereOptimize = function(par, fn, neighbor = NULL, ...){
  # check if initial value is on the unit sphere
  if(round(sum(par^2),10) != 1){
    stop('Initial value is not on a unit sphere.')
  }

  # create a temporary function taking spherical coordinates as inputs
  temp_fn = function(t){
    s = from.Sphere(t)
    res = fn(s)
    return(res)
  }
  # optimize over the whole unit sphere
  if(is.null(neighbor)){
    theta = to.Sphere(par)
    k = stats::optim(theta, temp_fn, ...)
    method = "Nelder-Mead"
  } else { # optimize over a neighbor of the initial value
    theta = to.Sphere(par)
    low_bd = theta - rep(neighbor, length(theta))
    up_bd = theta + rep(neighbor, length(theta))
    k = stats::optim(par, temp_fn, method = "L-BFGS-B", lower = low_bd, upper = up_bd, ...)
    method = "L-BFGS-B"
  }
  par = from.Sphere(k$par) # convert from spherical coordinates to Eculidean coordinates
  value = k$value
  return(list(
    par = par,
    value = value,
    method = method
  ))
}
