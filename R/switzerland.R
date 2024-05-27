#' Convert Swiss grid coordinates
#'
#' Convert coordinates from the Swiss coordinate system LV03 into global 
#' longitude and latidute values.
#'
#' @note
#' The Swiss coordinate system uses x to refer to the North-South axis, while 
#' y designates the East-West axis.
#'
#' @param x Coordinate on the Swiss grid in km.
#' @param y As *x*.
#' @return lat, lon Longitude and latitude in degrees.
#'
swiss_coords_to_lat_lon = function(x, y) {
  # Shift origin and convert units
  x_c = (y - 2e5) / 1e6
  y_c = (x - 6e5) / 1e6
  lon = 2.6779094 + 4.728982*y_c + 0.791484*y_c*x_c + 0.1306*y_c*x_c^2 - 
        0.0436*x_c^3
  lat = 16.9023892 + 3.238272*x_c - 0.270978*y_c^2 - 0.002528*x_c^2 -
        0.0447*y_c^2*x_c - 0.014*y_c^3
  conversion = 100./36.
  return(conversion * c(lat, lon))
}

