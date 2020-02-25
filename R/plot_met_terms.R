#' plot_met_terms
#' 
#' Plots the meteological variable terms extracted from extract_NEON_fluxes
#'
#' @param data \code{data.frame} containing the meteorological variables that should be plotted. Requires
#' \code{data.frame} to be the expanded version from extract_NEON_fluxes.
#'
#'
#' @return Multipanel plot showing Tair, RH, VPD, PAR, and Ustar.
#' 
#' @author Rich Fiorella \email{rich.fiorella@@utah.edu}
#' 
#' @export
#'
#'
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
  time.origin <- as.POSIXct(paste(data$Year[1],data$DoY[1]+data$Hour[1]/24),format="%Y %j",tz=attr(data,"tzone"))

  data$Time <- time.origin + (data$DoY-1)*86400 + data$Hour*3600

  print(data$Time)
  
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