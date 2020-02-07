#' Title
#'
#' @param data Input dataframe. Must have all 9 terms (3 CO2, 3 water, 3 temperature) from the eddy covariance data product in the timeseries.
#'
#'
#' @return Multipanel plot showing all 9 terms (3 net exchange terms, 3 turbulent flux terms, 3 storage terms.)
#' @export
#'
#' @examples
plot_met_terms <- function(data) {
  
  # list of required packages.
  require(ggplot2)
  require(dplyr)
  require(gridExtra)
  
  # check to make sure all required terms are present.
  # note: requires columns to have specific names! easiest
  # way to ensure this script works is to feed in an extended
  # data set calculated by extract_NEON_fluxes.
  
  if (!all(c("Tair","rH","VPD","PAR","Ustar","Year","DoY","Hour") %in% names(data))) {
    stop("Check to see if the correct data frame is being fed to this function!")
  }
  
  # change time to something posixct can handle.
  ddoy <- data$DoY + data$Hour/24
  data$Time <- as.POSIXct(paste0(data$Year,"-",ddoy),format="%Y-%j")
  
  # set up ggplot.
  p1 <- ggplot(data=data,aes(x=Time,y=Tair)) +
    geom_line() + 
    theme_bw() +
    scale_x_datetime(name="",expand=c(0,0))
  
  p2 <- ggplot(data=data,aes(x=Time,y=rH)) +
    geom_line() + 
    theme_bw() +
    scale_x_datetime(name="",expand=c(0,0))
  
  p3 <- ggplot(data=data,aes(x=Time,y=VPD)) +
    geom_line() + 
    theme_bw() +
    scale_x_datetime(name="",expand=c(0,0))
  
  p4 <- ggplot(data=data,aes(x=Time,y=PAR)) +
    geom_line() + 
    theme_bw() +
    scale_x_datetime(name="",expand=c(0,0))
  
  p5 <- ggplot(data=data,aes(x=Time,y=Ustar)) +
    geom_line() + 
    theme_bw() +
    scale_x_datetime(expand=c(0,0))
  
  grid.arrange(p1,p2,p3,p4,p5,nrow=5)
  
}