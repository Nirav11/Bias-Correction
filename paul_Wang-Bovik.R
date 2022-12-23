###################################################################################### Wang_Bovik begins
#' Calculate the Wang-Bovik index
#' 
#' This function calculated the Wang-Bovik index of two series with identical length
#' 
#' @param y Dataset 1 (length=n)
#' @param x DAtaset 2 (length=n)
#'
#' @return Returns Wang-Bovik index and its three components (Mxy, Vxy, and Rxy).  
#' Mo et al 2014 equation 28.

#' @references
#' Mo,R., Ye, C., and Whitfield, P.H. 2014. Application potential of four nontraditional 
#' similarity metrics in hydrometeorology. Journal of Hydrometerology 15:1862-1880.
#'
#' @export 
#' @examples
#' Wang-Bovik (x,y)

Wang_Bovik<- function (y,x) {
  
  # Get the joint minimum of x and y
  xy_min = min(c(min(x), min(y)))
  
  # Manipulate x
  xm = mean(x)
  xd = x - xm
  sxxd = mean(xd^2)
  
  # Manipulate y
  ym = mean(y)
  yd = y - ym
  syyd = mean(yd^2)
  
  sxyd = mean(xd*yd)
  xn = xm - xy_min
  yn = ym - xy_min
  
  # Calulate the three components of the modified WB index defined in Eq.(28)
  mxy = 2 * xn * yn / (xn^2 + yn^2)
  vxy = 2 * sqrt(sxxd * syyd) / (sxxd + syyd)
  rxy = sxyd / sqrt(sxxd * syyd)
  
  # Calculate the modified WB index defined in Eq.(28)
  wb_index = mxy * vxy * rxy
  
  # wb_index can also be computed directly as follows
  #wb_index = 4 * xn * yn * sxyd / ((xn*xn + yn*yn) * (sxxd + syyd))
  #-----------------------------------------------------------------------------------------
  
  result <-list(wb_index,mxy,vxy,rxy)
  names(result) <-c("Wang-Bovik Index", "Mxy", "Vxy", "Rxy")
  
  return (result)
}

###################################################################################### Wang_Bovik ends