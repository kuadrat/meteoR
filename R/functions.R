#' Inverse relative distance between Earth and Sun
#'
#' @param DOY day of year
inverse_relative_distance = function(DOY) {
  1. + 0.033 * cos(terran_angular_velocity * DOY)
}

#' Get solar declination for given DOY
#'
#' The solar declination starts from 0 degree in the March equinox (March 20) 
#' and increases to its maximum value of the Earth's tilt angle of 23.4 
#' degrees in the Summer solstice (June 20/21). From there it falls again, 
#' crossing 0 during the September equinox (September 22/23) before 
#' approaching its minimum of -23.4 degree on the Winter solstice (December 
#' 21/22).
#' @param DOY day of year
#' @return declination Solar declination in radians
#'
solar_declination = function(DOY) {
  # 23.4 degree = 0.408 RAD
  # The DOY for March 20 is 31 + 28 + 20 = 79
  return(0.408 * sin(terran_angular_velocity * (DOY - 79)))
}

#' Calculate sunset hourangle for given DOY at given latitude.
#'
#' @param DOY day of year
#' @param latitude Geographical latitude in degrees
#' @return Hourangle in radians
#'
sunset_hourangle = function(DOY, latitude) {
  declination = solar_declination(DOY)
  return(acos(-tan(latitude * radian) * tan(declination)))
}

#' Extraterrestrial radiation at given DOY at given latitude.
#'
#' @param DOY day of year
#' @param lat Geographical latitude in degree
#' @return Average extraterrestrial radiation in J/s/m^2, sometimes called 
#'   srad.  
#'
extraterrestrial_radiation = function(DOY, lat) {
  d = solar_declination(DOY)
  hourangle = sunset_hourangle(DOY, lat)
  term1 = 1 / pi * solar_constant * inverse_relative_distance(DOY)
  term2 = (hourangle * sin(lat * radian) * sin(d) + 
           cos(lat * radian) * cos(d) * sin(hourangle))
  return(term1 * term2)
}

#' Calculate the day length at given DOY at given latitude.
#'
#' @param DOY day of year
#' @param latitude Geographical latitude in degrees
#' @return Day length hours
#'
day_length = function(DOY, latitude) {
  hourangle = sunset_hourangle(DOY, latitude)
  # Factor 2 because sunrise and sunset are symmetric around noon.
  # THe angular velocity of Earth's rotation around itself is roughly 15 deg 
  # per hour.
  return( 2 * hourangle / (15*radian) )
}

#' Convert irradiance to photorsynthetically active radiation
#'
#' Integrates the daily average irradiation over a full day 
#' and applies a fractional conversion factor.
#'
#' @param srad Average sunlight irradiance in J/s/m^2
#' @param fraction Portion of *srad* that is photosynthetically active
#' @return PAR in MJ/m^2
#'
srad_to_PAR = function(srad, fraction = 0.47) {
  # There are 86400 seconds in a day; 1e-6 is for unit conversion.
  return(srad * fraction * 86400 * 1e-6)
}

#' Portion of extraterrestrial radiation that reaches Earth's surface
#'
#' Implements the Angstrom-Prescott formula for estimation of terrestrial 
#' radiation based on extraterrestrial radiation and relative sunshine 
#' duration.  
#' The Rietveld coefficients a and b are taken from 
#' Rietveld, M. R. A New Method for Estimating the Regression Coefficients in 
#' the Formula Relating Solar Radiation to Sunshine. Agricultural Meteorology 
#' 1978, 19 (2), 243–252. https://doi.org/10.1016/0002-1571(78)90014-6.
#'
#' @param r_extraterrestrial Average extraterrestrial radiation in J/s/m^2.
#' @param relative_sunshine_duration Expressed as fraction of a day.
#' @return Average terrestrial radiation in W/m^2
#'
terrestrial_radiation = function(r_extraterrestrial, 
                                 relative_sunshine_duration) {
  rssd = relative_sunshine_duration
  a = 0.10 + 0.24 * rssd
  b = 0.78 - 0.44 * rssd
  return(1.07 * (a + b*rssd) * r_extraterrestrial)
}

#' FAO Hargreaves equation
#'
#' Estimate reference evapotranspiration based on temperature and radiation 
#' data. See equation 52 in Chapter 3 of "Crop evapotranspiration - 
#' Guidelines for computing crop water requirements - FAO Irrigation and 
#' drainage paper 56"
#'
#' @param T_mean Daily average temperature in °C
#' @param T_max Daily maximum temperature in °C
#' @param T_min Daily minimum temperature in °C
#' @param r_extraterrestrial Daily average extraterrestrial radiation in W/m^2
#'
hargreaves_et0 = function(T_mean, T_max, T_min, r_extraterrestrial) {
  # Integrate extraterrestrial radiation over a day (86400 seconds) and apply 
  # a unit conversion from MJ/m^2/d to mm/day (energy for evaporating x mm H2O).
  r_integrated = r_extraterrestrial * 86400 / 2.45e6 
  return(0.0023 * (T_mean + 17.8) * (T_max - T_min)^0.5 * r_integrated)
}

#' Original Hargreaves-Samani equation
#'
#' @param T_mean Daily average temperature in °C
#' @param r_terrestrial Daily average terrestrial radiation in W/m^2
#' @return Reference evapotranspiration in mm/day
#'
hargreaves_et0_terrestrial = function(T_mean, r_terrestrial) {
  # Convert radiation in W/m2 to mm evaporated H2O per day.
  r_integrated = r_terrestrial * 86400 / 2.45e6
  return(0.0135 * r_integrated * (T_mean + 17.8))
}

#' FAO Hargreaves equation, alternative implementation
#'
#' Estimate reference evapotranspiration based on temperature and radiation 
#' data. See equation 52 in Chapter 3 of "Crop evapotranspiration - 
#' Guidelines for computing crop water requirements - FAO Irrigation and 
#' drainage paper 56"
#'
#' @param T_mean Daily average temperature in °C
#' @param T_max Daily maximum temperature in °C
#' @param T_min Daily minimum temperature in °C
#' @param r_extraterrestrial Daily average extraterrestrial radiation in W/m^2
#'
hargreaves_et0_extraterrestrial = function(T_mean, T_max, T_min, 
                                           r_extraterrestrial) {
  r_terrestrial = r_extraterrestrial * 0.17 * (T_max - T_min)^0.5
  return(hargreaves_et0_terrestrial(T_mean, r_terrestrial))
}

#' Convert Watt/m2 to mm evaporated water per day
#'
#' @param radiation Radiation in W/m2
#' @return radiative evapotranspiration in mm H2O per day
#'
radiation_to_evapotranspiration_unit_conversion = function(radiation) {
  return(radiation * 86400 * 0.408e-6)
}

