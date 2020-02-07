#' Title
#'
#' @param data Input dataframe. Must have all 9 terms (3 CO2, 3 water, 3 temperature) from the eddy covariance data product in the timeseries.
#'
#'
#' @return Multipanel plot showing all 9 terms (3 net exchange terms, 3 turbulent flux terms, 3 storage terms.)
#' @export
#'
#' @examples
plot_storage_vs_ustar <- function(data) {
  
  # list of required packages.
  require(ggplot2)
  require(dplyr)
  require(gridExtra)
  
  # check to make sure all required terms are present.
  # note: requires columns to have specific names! easiest
  # way to ensure this script works is to feed in an extended
  # data set calculated by extract_NEON_fluxes.
  
  if (!all(c("Sc","Sw","ST","Ustar","Year","DoY","Hour") %in% names(data))) {
    stop("Check to see if the correct data frame is being fed to this function!")
  }
  
  # change time to something posixct can handle.
  ddoy <- data$DoY + data$Hour/24
  data$Time <- as.POSIXct(paste0(data$Year,"-",ddoy),format="%Y-%j")
  
  # set up ggplot.
  p1 <- ggplot(data=data,aes(x=Ustar,y=Sc)) +
    geom_hex() + 
    theme_bw() +
    theme(legend.position = "bottom")
  
  p2 <- ggplot(data=data,aes(x=Ustar,y=Sw)) +
    geom_hex() + 
    theme_bw() +
    theme(legend.position = "bottom")
  
  p3 <- ggplot(data=data,aes(x=Ustar,y=ST)) +
    geom_hex() + 
    theme_bw() +
    theme(legend.position = "bottom")

  grid.arrange(p1,p2,p3,nrow=1)
  
}