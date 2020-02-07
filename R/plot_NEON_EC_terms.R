#' Title
#'
#' @param data Input dataframe. Must have all 9 terms (3 CO2, 3 water, 3 temperature) from the eddy covariance data product in the timeseries.
#'
#'
#' @return Multipanel plot showing all 9 terms (3 net exchange terms, 3 turbulent flux terms, 3 storage terms.)
#' @export
#'
#' @examples
plot_NEON_EC_terms <- function(data) {
  
  # list of required packages.
  require(ggplot2)
  require(dplyr)
  require(gridExtra)
  
  # check to make sure all required terms are present.
  # note: requires columns to have specific names! easiest
  # way to ensure this script works is to feed in an extended
  # data set calculated by extract_NEON_fluxes.
  
  if (!all(c("NEE","LH","H","Fc","Fw","FT","Sc","Sw","ST","Year","DoY","Hour") %in% names(data))) {
    stop("Check to see if the correct data frame is being fed to this function!")
  }
  
  # change time to something posixct can handle.
  ddoy <- data$DoY + data$Hour/24
  data$Time <- as.POSIXct(paste0(data$Year,"-",ddoy),format="%Y-%j")
  
  # set up ggplot.
  p1 <- ggplot(data=data,aes(x=Time,y=NEE)) +
    geom_line() + 
    theme_bw()
  
  p2 <- ggplot(data=data,aes(x=Time,y=Fc)) +
    geom_line() + 
    theme_bw()
  
  p3 <- ggplot(data=data,aes(x=Time,y=Sc)) +
    geom_line() + 
    theme_bw()
  
  p4 <- ggplot(data=data,aes(x=Time,y=LH)) +
    geom_line(col="red") + 
    theme_bw() +
    scale_x_datetime(name="")
  
  p5 <- ggplot(data=data,aes(x=Time,y=Fw)) +
    geom_line(col="red") + 
    theme_bw() +
    scale_x_datetime(name="")
  
  p6 <- ggplot(data=data,aes(x=Time,y=Sw)) +
    geom_line(col="red") + 
    theme_bw() +
    scale_x_datetime(name="")
  
  p7 <- ggplot(data=data,aes(x=Time,y=H)) +
    geom_line(col="blue") + 
    theme_bw() +
    scale_x_datetime(name="")
  
  p8 <- ggplot(data=data,aes(x=Time,y=FT)) +
    geom_line(col="blue") + 
    theme_bw() +
    scale_x_datetime(name="")
  
  p9 <- ggplot(data=data,aes(x=Time,y=ST)) +
    geom_line(col="blue") + 
    theme_bw() +
    scale_x_datetime(name="")
  
  grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,nrow=9)
  
}