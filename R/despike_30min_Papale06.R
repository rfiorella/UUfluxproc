#' Title
#'
#' @param flux.df 
#' @param z 
#'
#' @return
#' @export
#'
#' @examples
despike_30min_Papale06 <- function(flux.df,z=4) {
  
  # remove spikes from a timeseries of 30 minute fluxes using the 
  # method described in section 2.2 of Papale et al. 2006, Biogeosceinces
  
  # # check to see if there are Rg, nee, LH, H columns. 
  # if (!(c("Rg","NEE","LH","H") %in% colnames(flux.df))) {
  #   print("missing variables in data frame.")
  # }
  
  #---------------------------------------------
  # calculate "d" value.
  d <- vector() # initialize d.
  
  # take double difference to calculate d.
  for (i in 1:nrow(flux.df)) {
    if (i == 1 | i == nrow(flux.df)) {
      d[i] <- NA
    } else {
      d[i] <- (flux.df$NEE[i]-flux.df$NEE[i-1])-(flux.df$NEE[i+1]-flux.df$NEE[i])
    }
  }
  
  #---------------------------------------------
  # calculate rolling median of NEE
  Md <- median(d,na.rm=TRUE)
  
  #----------------------------------------------
  # Calculate MAD values.
  
  MAD <- rollapply(abs(d-Md),3,median,na.rm=TRUE,fill=NA)
  
  #----------------------------------------------
  # Calculate upper/lower bounds
  
  upr.limit <- Md + (z*MAD/0.6745)
  lwr.limit <- Md - (z*MAD/0.6745)

  #----------------------------------------------
  # get logical vector of where d < lwr limit or d > upr limit
  
  spikes <- ((d < lwr.limit) | (d > upr.limit))
  
  print(length(spikes))
  print(nrow(flux.df))
  #----------------------------------------------
  # set spikes to NA, and return flux.df back to workspace
  print(paste(100*sum(spikes,na.rm=TRUE)/length(spikes)," % of NEE values appear to be spikes."))
  
  flux.df$NEE[spikes == TRUE] <- NA
  
  # for debugging, return a dataframe of all of these metrics.
  out.df <- data.frame(NEE=flux.df$NEE,d,Md,spikes,lwr.limit,upr.limit,MAD,
                       NEE3=rollapply(flux.df$NEE,3,median,na.rm=TRUE,fill=NA))
  
  return(out.df)
  
}