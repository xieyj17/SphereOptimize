#' @description convert spherical coordinates to Cartesian coordinates at each coordinate
.med.From.Sphere = function(theta, i){
  res = 1
  if(i <= length(theta)){
    temp_theta = theta[1:i]
    if(i == 1){
      res = cos(temp_theta)
    } else {
      for(p in 1:(i-1)){
        res = res * sin(temp_theta[p])
      }
      res = res * cos(temp_theta[i])
    }
  } else {
    temp_theta = theta
    for(p in 1:(i-2)){
      res = res * sin(temp_theta[p])
    }
    res = res * sin(temp_theta[i-1])
  }
  return(res)
}

#' @description make sure the angles are in correct range
.regularize.Theta = function(ang){
  while(ang > 2*pi){
    ang = ang - 2*pi
  }
  while(ang < 0){
    ang = ang + 2*pi
  }
  return(ang)
}
