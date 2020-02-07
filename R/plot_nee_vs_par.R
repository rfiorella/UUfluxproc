#' Title
#'
#' @param data Input dataframe. Must have all 9 terms (3 CO2, 3 water, 3 temperature) from the eddy covariance data product in the timeseries.
#'
#'
#' @return Multipanel plot showing all 9 terms (3 net exchange terms, 3 turbulent flux terms, 3 storage terms.)
#' @export
#'
#' @examples
plot_nee_vs_par <- function(data,nee.thres=50) {
  
  # list of required packages.
  require(ggplot2)
  require(dplyr)
  require(gridExtra)
  
  # check to make sure all required terms are present.
  # note: requires columns to have specific names! easiest
  # way to ensure this script works is to feed in an extended
  # data set calculated by extract_NEON_fluxes.
  
  if (!all(c("NEE","PAR","Year","DoY","Hour") %in% names(data))) {
    stop("Check to see if the correct data frame is being fed to this function!")
  }
  
  # change time to something posixct can handle.
  ddoy <- data$DoY + data$Hour/24
  data$Time <- as.POSIXct(paste0(data$Year,"-",ddoy),format="%Y-%j")
  
  # set up ggplot.
  p1 <- ggplot(data=data,aes(x=PAR,y=NEE)) +
    geom_point(alpha = .1) + 
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_y_continuous(limits=c(-nee.thres,nee.thres))
  
  print(p1)
  
}