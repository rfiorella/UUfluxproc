#' extract_NEON_PMvars
#'
#' Extracts NEON eddy covariance data from HDF5 files, merges with relevant meteorolgical data,
#' and returns a \code{data.frame} for a specified time period. A specific year can be requested,
#' though the default is to return the entire period of record for that site. Current setup 
#' assumes that you have a local copy of eddy covariance files (DP4.00200.001) for the site of interest.
#' Meterological data *can* also be local, and a location specified by the \code{met.path} argument; if
#' no value is provided to \code{met.path}, meteorological data for the site for the requested period
#' will be retrieved from the NEON API. Clear spikes are removed using a moving "deviation from median" 
#' filter, described by Brock 1986 and shown to perform well in Starkenburg et al. 2016.
#' Default options return a \code{data.frame} that is immediately useable for flux partitioning using 
#' \code{REddyProc}, though users are encouraged to verify that the despiking worked appropriately for
#' their site. An option to save the returned \code{data.frame}, set by the arguments 
#' \code{write.to.file == TRUE} and \code{out.path}, is also available.
#' 
#' @section Notes on median deviation filter:
#' This function by default applies a "median absolute deviation" filter to the net flux variables
#' of sensible heat, latent heat, and net ecosystem exchange. The filter is applied to 
#' a moving window of data, currently 7 points, and each point in that window is compared to the median
#' of the window. If it differs by the median by a certain threshold, the data point is discarded. Due 
#' to the width of the filter, only spikes smaller than 3 points can be removed. As a result, there is 
#' a second weak filter that removes points that are more than 4sigma from the mean. The threshold values
#' in this filter are 10 umolCO2/m2/s for NEE, and 100 W/m2 for LH and H. All of these should be considered
#' provisional, and different thresholds may be suggested. If you need different thresholds, or suggest they
#' should be changed globally, contact me.
#'
#' @author Rich Fiorella \email{rich.fiorella@@utah.edu}
#' 
#' @param neon.site Four letter code indicating NEON site.
#' @param year Which year to process? If not specified, process all. 
#' @param flux.path Specify path to flux data you wish to process.
#' @param met.path  Specify path to meteorological data (if NULL - download from NEON API)
#' @param median.filter Filter half-hourly data using the Brock 86 median filter.
#' @param out.path If saving a file of results, where should it be saved?
#' @param write.to.file Write to csv file?
#'
#' @return A \code{data.frame} containing merged eddy covariance and meteorological variables 
#'         on a common time array. If \code{expanded = FALSE}, output will be exactly the input
#'         \code{data.frame}. required by \code{REddyProc}. If \code{expanded = TRUE}, 
#'         output \code{data.frame} includes: a) storage and turblent terms of net flux, b) expanded
#'         radiation variables, including outgoing SW, incoming and outgoing LW, and PAR.
#'         
#' @export
#' @import neonUtilities
#' @import rhdf5
#' @import xts
#' @import lubridate
#' @import tidyverse
extract_NEON_PMvars <- function(neon.site,
                                year=9999,
                                flux.path="~/Dropbox/NEON/DP4_00200_001",
                                met.path,
                                median.filter=TRUE,
                                write.to.file=TRUE,
                                out.path) {
  
  # stack flux data.
  fluxes <- neonUtilities::stackEddy(paste0(flux.path,"/",neon.site),level="dp04")
  
  fluxes.flat <- fluxes[[neon.site]] # flatten list structure.

  # extract required variables.
  fluxes.reduced <- fluxes.flat %>%
    select(timeBgn,timeEnd, # time variables
           data.fluxCo2.nsae.flux,data.fluxCo2.stor.flux,data.fluxCo2.turb.flux, # CO2 fluxes
           data.fluxH2o.nsae.flux,data.fluxH2o.stor.flux,data.fluxH2o.turb.flux, # H2O fluxes,
           data.fluxTemp.nsae.flux,data.fluxTemp.stor.flux,data.fluxTemp.turb.flux, # Temp fluxes,
           data.fluxMome.turb.veloFric) %>% # u*
    rename(nee=data.fluxCo2.nsae.flux,lhf=data.fluxH2o.nsae.flux,
           shf=data.fluxTemp.nsae.flux,ustar=data.fluxMome.turb.veloFric,
           Fc=data.fluxCo2.turb.flux,Fw=data.fluxH2o.turb.flux,
           Ft=data.fluxTemp.turb.flux,Sc=data.fluxCo2.stor.flux,
           Sw=data.fluxH2o.stor.flux,St=data.fluxTemp.stor.flux)
  
  # cut down further if not expanded
    fluxes.reduced <- fluxes.reduced %>%
      select(timeBgn,timeEnd,nee,lhf,shf,ustar)
  
  # convert "NEON" time to POSIXct
  fluxes.reduced$timeBgn <- as.POSIXct(fluxes.reduced$timeBgn,format="%Y-%m-%dT%H:%M:%S.%OSZ",tz="UTC")
  fluxes.reduced$timeEnd <- as.POSIXct(fluxes.reduced$timeEnd,format="%Y-%m-%dT%H:%M:%S.%OSZ",tz="UTC")

  # convert fluxes to xts. 
  fluxes.notime <- fluxes.reduced %>%
    select(-timeBgn,timeEnd)
  
  flux.xts <- as.xts(fluxes.notime,order.by=fluxes.reduced$timeBgn)
  
  # get lat/lon and timezone of station.
  slist <- list.files(path=paste0(flux.path,"/",neon.site,"/"),pattern=".h5",
                      full.names=TRUE,recursive=TRUE)
  
  attrs <- h5readAttributes(slist[[1]],neon.site)
  
  lat <- as.numeric(attrs$LatTow)
  lon <- as.numeric(attrs$LonTow)
  tzone <- as.character(attrs$ZoneTime)
  hgts <- as.numeric(attrs$DistZaxsLvlMeasTow)
  
  # get h2o profile
  h2oprof.tmp <- neonUtilities::stackEddy(paste0(flux.path,"/",neon.site),level="dp01",avg=30,var="rtioMoleWetH2o")
  
  # simplify h2o profile
  h2oprof <- h2oprof.tmp[[1]] %>%
    select(verticalPosition,timeBgn,data.h2oStor.rtioMoleWetH2o.mean) %>%
    rename(h2o = data.h2oStor.rtioMoleWetH2o.mean,time = timeBgn) %>%
    mutate(verticalPosition = as.numeric(verticalPosition)/10) %>%
    filter(!is.na(verticalPosition)) %>%
    mutate(height = hgts[verticalPosition]) %>%
    select(-verticalPosition) %>%
    pivot_wider(names_from=height,
                names_prefix="H2O_",
                names_sep="_",
                values_from=h2o)
  
  h2oprof$time <- as.POSIXct(h2oprof$time,format="%Y-%m-%dT%H:%M:%S.%OSZ",tz="UTC")
  
  # get co2 profile
  co2prof.tmp <- neonUtilities::stackEddy(paste0(flux.path,"/",neon.site),level="dp01",avg=30,var="rtioMoleDryCo2")
  
  # simplify co2 profile
  co2prof <- co2prof.tmp[[1]] %>%
    select(verticalPosition,timeBgn,data.co2Stor.rtioMoleDryCo2.mean) %>%
    rename(co2 = data.co2Stor.rtioMoleDryCo2.mean, time = timeBgn) %>%
    mutate(verticalPosition = as.numeric(verticalPosition)/10) %>%
    filter(!is.na(verticalPosition)) %>%
    mutate(height = hgts[verticalPosition]) %>%
    select(-verticalPosition) %>%
    pivot_wider(names_from=height,
                names_prefix="CO2_",
                values_from=co2)
  
  co2prof$time <- as.POSIXct(co2prof$time,format="%Y-%m-%dT%H:%M:%S.%OSZ",tz="UTC")

  #------------------------------------------------------------
  # load met data.
  #------------------------------------------------------------
  
  if (year > 2016) {
    
    # set start month.
    start.mon <- paste0(year,"-01")
    end.mon   <- paste0(year,"-12")
    
    # need, at a minimum, RH, Rg, and Tair.
    Rh.tmp <- neonUtilities::loadByProduct("DP1.00098.001",site=neon.site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
    Rg.tmp <- neonUtilities::loadByProduct("DP1.00023.001",site=neon.site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
    Ttrip.tmp <- neonUtilities::loadByProduct("DP1.00003.001",site=neon.site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
    Tsing.tmp <- neonUtilities::loadByProduct("DP1.00002.001",site=neon.site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
    ws.tmp <- neonUtilities::loadByProduct("DP1.00001.001",site=neon.site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
    bp.tmp <- neonUtilities::loadByProduct("DP1.00004.001",site=neon.site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
    soihf.tmp <- neonUtilities::loadByProduct("DP1.00040.001",site=neon.site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
    
  } else {
    stop("Year selected predates NEON operation.")
  }
  
  #====================================
  # pull out and simplify met variables.
  #------------------------------------
  # relative humidity
  Rh1 <- Rh.tmp$RH_30min %>% 
    select(RHMean,startDateTime,horizontalPosition,verticalPosition) %>%
    mutate(HOR.VER = paste(horizontalPosition,verticalPosition,sep=".")) %>%
    select(-horizontalPosition,-verticalPosition)
  Rh2 <- Rh.tmp$sensor_positions_00098 %>%
    select(HOR.VER,xOffset,yOffset,zOffset)
  
  Rh <- left_join(Rh1,Rh2,by="HOR.VER") %>%
    mutate(dist.from.tower = sqrt(xOffset^2 + yOffset^2)) %>%
    filter(dist.from.tower < 30) %>%
    select(-HOR.VER,-xOffset,-yOffset,-dist.from.tower) %>%
    rename(time = startDateTime) %>%
    pivot_wider(names_from=zOffset,
                names_prefix="RH",
                names_sep="_",
                values_from=RHMean)
  
  rm(Rh1, Rh2, Rh.tmp)
  
  # single-aspirated air temperature
  Tsing1 <- Tsing.tmp$SAAT_30min %>%
    select(tempSingleMean,startDateTime,horizontalPosition,verticalPosition) %>%
    mutate(HOR.VER = paste(horizontalPosition,verticalPosition,sep=".")) %>%
    select(-horizontalPosition,-verticalPosition)
  Tsing2 <- Tsing.tmp$sensor_positions_00002 %>%
    select(HOR.VER,xOffset,yOffset,zOffset)
  
  Tsing <- left_join(Tsing1,Tsing2,by="HOR.VER") %>%
    mutate(dist.from.tower = sqrt(xOffset^2 + yOffset^2)) %>%
    filter(dist.from.tower < 30) %>%
    select(-HOR.VER,-xOffset,-yOffset,-dist.from.tower) %>%
    rename(time = startDateTime) %>%
    pivot_wider(names_from=zOffset,
                names_prefix="singTemp",
                names_sep="_",
                values_from=tempSingleMean)
  
  rm(Tsing1, Tsing2, Tsing.tmp)
  
  # triple-aspirated air temperature
  Ttrip1 <- Ttrip.tmp$TAAT_30min %>%
    select(tempTripleMean,startDateTime,horizontalPosition,verticalPosition) %>%
    mutate(HOR.VER = paste(horizontalPosition,verticalPosition,sep=".")) %>%
    select(-horizontalPosition,-verticalPosition)
  Ttrip2 <- Ttrip.tmp$sensor_positions_00003 %>%
    select(HOR.VER,xOffset,yOffset,zOffset)
  
  Ttrip <- left_join(Ttrip1,Ttrip2,by="HOR.VER") %>%
    mutate(dist.from.tower = sqrt(xOffset^2 + yOffset^2)) %>%
    filter(dist.from.tower < 30) %>%
    select(-HOR.VER,-xOffset,-yOffset,-dist.from.tower) %>%
    rename(time = startDateTime) %>%
    pivot_wider(names_from=zOffset,
                names_prefix="tripTemp",
                names_sep="_",
                values_from=tempTripleMean)
  
  rm(Ttrip1, Ttrip2, Ttrip.tmp)
  
  # barometric pressure
  bp1 <- bp.tmp$BP_30min %>%
    select(staPresMean,startDateTime,horizontalPosition,verticalPosition) %>%
    mutate(HOR.VER = paste(horizontalPosition,verticalPosition,sep=".")) %>%
    select(-horizontalPosition,-verticalPosition)
  bp2 <- bp.tmp$sensor_positions_00004 %>%
    select(HOR.VER,xOffset,yOffset,zOffset)
  
  bp <- left_join(bp1,bp2,by="HOR.VER") %>%
    mutate(dist.from.tower = sqrt(xOffset^2 + yOffset^2)) %>%
    filter(dist.from.tower < 30) %>%
    select(-HOR.VER,-xOffset,-yOffset,-dist.from.tower) %>%
    rename(time = startDateTime) %>%
    pivot_wider(names_from=zOffset,
                names_prefix="pres",
                names_sep="_",
                values_from=staPresMean)
  
  rm(bp1, bp2, bp.tmp)
  
  # wind speed
  ws1 <- ws.tmp$`2DWSD_30min` %>%
    select(windSpeedMean,startDateTime,horizontalPosition,verticalPosition) %>%
    mutate(HOR.VER = paste(horizontalPosition,verticalPosition,sep=".")) %>%
    select(-horizontalPosition,-verticalPosition)
  ws2 <- ws.tmp$sensor_positions_00001 %>%
    select(HOR.VER,xOffset,yOffset,zOffset)
  
  ws <- left_join(ws1,ws2,by="HOR.VER") %>%
    mutate(dist.from.tower = sqrt(xOffset^2 + yOffset^2)) %>%
    filter(dist.from.tower < 30) %>%
    select(-HOR.VER,-xOffset,-yOffset,-dist.from.tower) %>%
    rename(time = startDateTime) %>%
    pivot_wider(names_from=zOffset,
                names_prefix="ws",
                names_sep="_",
                values_from=windSpeedMean)
  
  rm(ws1, ws2, ws.tmp)
  
  # radiation terms
  Rg1 <- Rg.tmp$SLRNR_30min %>%
    select(startDateTime,horizontalPosition,verticalPosition,inSWMean,outSWMean,outLWMean,inLWMean) %>%
    mutate(HOR.VER = paste(horizontalPosition,verticalPosition,sep=".")) %>%
    select(-horizontalPosition,-verticalPosition)
  Rg2 <- Rg.tmp$sensor_positions_00023 %>%
    select(HOR.VER,xOffset,yOffset,zOffset)
  
  Rg <- left_join(Rg1,Rg2,by="HOR.VER") %>%
    mutate(dist.from.tower = sqrt(xOffset^2 + yOffset^2)) %>%
    filter(dist.from.tower < 30) %>%
    select(-HOR.VER,-xOffset,-yOffset,-dist.from.tower) %>%
    mutate(Rnet = (inSWMean + inLWMean) - (outSWMean + outLWMean)) %>%
    rename(inSW = inSWMean, inLW = inLWMean, outSW = outSWMean, outLW = outLWMean, time = startDateTime) %>%
    pivot_wider(names_from=zOffset,
                values_from=c(inSW,inLW,outSW,outLW,Rnet))
  
  rm(Rg1,Rg2,Rg.tmp)
  
  # soil heat flux
  soihf1 <- soihf.tmp$SHF_30min %>%
    select(startDateTime,horizontalPosition,verticalPosition,SHFMean) %>%
    mutate(HOR.VER = paste(horizontalPosition,verticalPosition,sep=".")) %>%
    select(-horizontalPosition,-verticalPosition)
  soihf2 <- soihf.tmp$sensor_positions_00040 %>%
    select(HOR.VER,xOffset,yOffset,zOffset)
  
  soihf <- left_join(soihf1,soihf2,by="HOR.VER") %>%
    mutate(dist.from.tower = round(sqrt(xOffset^2 + yOffset^2),3)) %>%
    filter(dist.from.tower < 30) %>%
    filter(zOffset == -0.08) %>%
    select(-xOffset,-yOffset,-zOffset,-dist.from.tower) %>%
    rename(time = startDateTime) %>%
    pivot_wider(names_prefix="SHF",
                names_from=HOR.VER,
                values_from=SHFMean)
  
  rm(soihf1,soihf2,soihf.tmp)
  
  # change time vars to character, then to posix ct.
  Rg$time <- as.POSIXct(Rg$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  Rh$time <- as.POSIXct(Rh$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  Ttrip$time <- as.POSIXct(Ttrip$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  ws$time <- as.POSIXct(ws$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  bp$time <- as.POSIXct(bp$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  Tsing$time <- as.POSIXct(Tsing$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  soihf$time <- as.POSIXct(soihf$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")

  # convert data to xts objects and then merge.
  Rg.xts <- xts(Rg[,2:ncol(Rg)],order.by=Rg$time)
  Rh.xts <- xts(Rh[,2:ncol(Rh)],order.by=Rh$time)
  ws.xts <- xts(ws[,2:ncol(ws)],order.by=ws$time)
  bp.xts <- xts(bp[,2:ncol(bp)],order.by=bp$time)
  Ttrip.xts <- xts(Ttrip[,2:ncol(Ttrip)],order.by=Ttrip$time)
  Tsing.xts <- xts(Tsing[,2:ncol(Tsing)],order.by=Tsing$time)
  soihf.xts <- xts(soihf[,2:ncol(soihf)],order.by=soihf$time)
  h2oprof.xts <- xts(h2oprof[,2:ncol(h2oprof)],order.by=h2oprof$time)
  co2prof.xts <- xts(co2prof[,2:ncol(co2prof)],order.by=co2prof$time)

  minTime <- as.POSIXct(min(c(min(index(Rg.xts)),min(index(Rh.xts)),
                              min(index(Ttrip.xts)),min(index(flux.xts)),
                              min(index(bp.xts)),min(index(ws.xts)),
                              min(index(Tsing.xts)),min(index(soihf.xts)),
                              min(index(h2oprof.xts)),min(index(co2prof.xts)))))
  
  maxTime <- as.POSIXct(max(c(max(index(Rg.xts)),max(index(Rh.xts)),
                              max(index(Ttrip.xts)),max(index(flux.xts)),
                              max(index(bp.xts)),max(index(ws.xts)),
                              max(index(Tsing.xts)),max(index(soihf.xts)),
                              max(index(h2oprof.xts)),max(index(co2prof.xts)))))
  
  dummy.ts <- seq.POSIXt(minTime,maxTime,by=1800)
  dummy.data <- rep(NA,length(dummy.ts))
  
  dummy.xts <- xts(dummy.data,order.by=dummy.ts)
  
  all.data <- merge.xts(dummy.xts,flux.xts,Rg.xts,Rh.xts,Tsing.xts,Ttrip.xts,
                        bp.xts,ws.xts,soihf.xts,
                        h2oprof.xts,co2prof.xts)

  
  if (year != 9999) {
    all.data <- all.data[paste0(start.mon,"/",end.mon)]  
  }
  
  if (tzone == "PST") {
    tzone(all.data) <- "Etc/GMT+8"  
  } else if (tzone == "MST") {
    tzone(all.data) <- "Etc/GMT+7"
  } else if (tzone == "CST") {
    tzone(all.data) <- "Etc/GMT+6" 
  } else if (tzone == "EST") {
    tzone(all.data) <- "Etc/GMT+5"
  }
  
  data.out <- data.frame(time=index(all.data),coredata(all.data)) %>%
    select(-timeEnd,-dummy.xts) %>%
    rename(H = shf)
  
   #------------------------------------------------------------
  
  if (median.filter == TRUE) {
    
    #---------------- Filter NEE ----------------
    # put in median deviation filter from Brock 86 / Starkenburg 2016.
    
    NEE.filt <- rollapply(data.out$nee,7,median,na.rm=TRUE,fill=NA)
    NEE.logi <- abs(data.out$nee - NEE.filt) > 10
    
    data.out$nee[NEE.logi == TRUE] <- NA
    
    # remove points that are 5 sigma away from mean?
    NEE.mu <- mean(data.out$nee,na.rm=TRUE)
    NEE.sd <- sd(data.out$nee,na.rm=TRUE)

    NEE.oor <- (data.out$nee < NEE.mu-4*NEE.sd) | (data.out$nee > NEE.mu+4*NEE.sd)

    # set out of range values to missing:
    data.out$nee[NEE.oor == TRUE] <- NA

    #---------------- Filter LH -----------------
    # median deviation filter of Brock 86 / Starkenburg 16
    
    LH.filt <- zoo::rollapply(data.out$lhf,7,median,na.rm=TRUE,fill=NA)
    LH.logi <- abs(data.out$lhf - LH.filt) > 100 # this threshold has not been checked!
    
    data.out$lhf[LH.logi == TRUE] <- NA # replace with missing
    
    # remove points that are 4 sigma away from mean?
    LH.mu <- mean(data.out$lhf,na.rm=TRUE)
    LH.sd <- sd(data.out$lhf,na.rm=TRUE)
    
    LH.oor <- (data.out$lhf < LH.mu-4*LH.sd) | (data.out$lhf > LH.mu+4*LH.sd)
    
    # set out of range values to missing:
    data.out$lhf[LH.oor == TRUE] <- NA
    
    #---------------- Filter H ------------------
    # median deviation filter of Brock 86 / Starkenburg 16
    
    H.filt <- zoo::rollapply(data.out$H,7,median,na.rm=TRUE,fill=NA)
    H.logi <- abs(data.out$H - H.filt) > 100 # this threshold has not been checked!
    
    data.out$H[H.logi == TRUE] <- NA # replace with missing
    
    # remove points that are 4 sigma away from mean?
    H.mu <- mean(data.out$H,na.rm=TRUE)
    H.sd <- sd(data.out$H,na.rm=TRUE)
    
    H.oor <- (data.out$H < H.mu-4*H.sd) | (data.out$H > H.mu+4*H.sd)
    
    # set out of range values to missing:
    data.out$H[H.oor == TRUE] <- NA
  }
 
  attr(data.out,"lat") <- lat 
  attr(data.out,"lon") <- lon 
  attr(data.out,"tzone") <- tzone
  
  dfDigits <- function(x, digits = 3) {
    ## x is a data.frame
    for (col in colnames(x)[sapply(x, class) == 'numeric'])
      x[,col] <- round(x[,col], digits = digits)
    x
  }
  
  # round data frame to 3 dps
  data.out <- dfDigits(data.out)
  #------------------------------------------------------------
  # write out data file if requested.
  if (write.to.file == TRUE) {
    if (is.null(out.path)) {
      stop("Attempting to write to file, but no path specified!")
    } else {
      write.csv(data.out,
                  paste0(out.path,"/",Sys.Date(),"_",neon.site,"_",year,"_PMinv.csv"),
                  row.names=FALSE,quote=FALSE)
    }
  }
    
# return(data.out)
}  
  
  
  