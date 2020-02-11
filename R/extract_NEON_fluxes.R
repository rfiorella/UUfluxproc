#' extract_NEON_fluxes
#'
#' @param neon.site Four letter code indicating NEON site.
#' @param year Which year to process? If not specified, process all. 
#' @param flux.path Specify path to flux data you wish to process.
#' @param met.path  Specify path to meteorological data (if NULL - download from NEON API)
#' @param expanded 
#' @param median.filter Filter half-hourly data using the Brock 86 median filter.
#' @param filt.width Width of median filter. Old implementation.
#' @param out.path If saving a file of results, where should it be saved?
#' @param write.to.file Write to csv file?
#' @param fix.tz Convert from UTC to local time zone?
#'
#' @return
#' @export
#'
extract_NEON_fluxes <- function(neon.site,
                                year=9999,
                                flux.path="~/Dropbox/NEON/DP4_00200_001",
                                met.path,
                                expanded=FALSE,
                                median.filter=TRUE,
                                filt.width=3,
                                fix.tz=FALSE,
                                write.to.file=FALSE,
                                out.path) {
  
  # list required packages
  require(rhdf5)
  require(neonUtilities)
  require(xts)
  require(lubridate)
  require(tidyverse)
  
  # stack flux data.
  fluxes <- stackEddy(paste0(flux.path,"/",neon.site),level="dp04")
  
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
  if (!expanded) {
    fluxes.reduced <- fluxes.reduced %>%
      select(timeBgn,timeEnd,nee,lhf,shf,ustar)
  }
  
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
  print(slist)
  
  attrs <- h5readAttributes(slist[[1]],neon.site)
  
  lat <- as.numeric(attrs$LatTow)
  lon <- as.numeric(attrs$LonTow)
  tzone <- as.character(attrs$ZoneTime)
  #------------------------------------------------------------
  # load met data.
  #------------------------------------------------------------
  
  if (year == 9999) { # get all years w/ flux data.
    # need, at a minimum, RH, Rg, and Tair.
    Rh.tmp <- loadByProduct("DP1.00098.001",site=neon.site,startdate="2017-01",avg=30,check.size=F)
    Rg.tmp <- loadByProduct("DP1.00023.001",site=neon.site,startdate="2017-01",avg=30,check.size=F)
    Ta.tmp <- loadByProduct("DP1.00003.001",site=neon.site,startdate="2017-01",avg=30,check.size=F)
    
    if (expanded == TRUE) {
      PAR.tmp <- loadByProduct("DP1.00024.001",site=neon.site,startdate="2017-01",avg=30,check.size=F)
    }
  } else if (year > 2015) {
    
    # set start month.
    start.mon <- paste0(year,"-01")
    end.mon   <- paste0(year,"-12")
    
    # need, at a minimum, RH, Rg, and Tair.
    Rh.tmp <- loadByProduct("DP1.00098.001",site=neon.site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
    Rg.tmp <- loadByProduct("DP1.00023.001",site=neon.site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
    Ta.tmp <- loadByProduct("DP1.00003.001",site=neon.site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
    
    if (expanded == TRUE) {
      PAR.tmp <- loadByProduct("DP1.00024.001",site=neon.site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
    }
  } else {
    stop("Year selected predates NEON operation.")
  }
  
  # pull out met variables.
  Rh <- Rh.tmp$RH_30min %>%
    filter(as.numeric(horizontalPosition) == 0) %>%
    filter(as.numeric(verticalPosition) == max(as.numeric(verticalPosition))) %>%
    select(RHMean,startDateTime) 
  Ta <- Ta.tmp$TAAT_30min %>% 
    filter(as.numeric(horizontalPosition) == 0) %>% 
    filter(as.numeric(verticalPosition) == max(as.numeric(verticalPosition))) %>%
    select(tempTripleMean,startDateTime)
  
  
  if (expanded == TRUE) {
    Rg <- Rg.tmp$SLRNR_30min %>%
      filter(as.numeric(horizontalPosition) == 0) %>%
      filter(as.numeric(verticalPosition) == max(as.numeric(verticalPosition))) %>%
      select(inSWMean,outSWMean,inLWMean,outLWMean,startDateTime) 
    
    PAR <- PAR.tmp$PARPAR_30min %>%
      filter(as.numeric(horizontalPosition) == 0) %>%
      filter(as.numeric(verticalPosition) == max(as.numeric(verticalPosition))) %>%
      select(PARMean,startDateTime)  
    
  } else {
    Rg <- Rg.tmp$SLRNR_30min %>%
      filter(as.numeric(horizontalPosition) == 0) %>%
      filter(as.numeric(verticalPosition) == max(as.numeric(verticalPosition))) %>%
      select(inSWMean,startDateTime) 
  }

  names(Rh) <- c("Rh","time")
  names(Ta) <- c("Ta","time")

  if (expanded == TRUE) {
    names(Rg) <- c("SWdown","SWup","LWdown","LWup","time")
    names(PAR) <- c("PAR","time")
  } else {
    names(Rg) <- c("Rg","time")
  }
  
  # change time vars to character, then to posix ct.
  Rg$time <- as.POSIXct(Rg$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  Rh$time <- as.POSIXct(Rh$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  Ta$time <- as.POSIXct(Ta$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  
  if (expanded == TRUE) {
    PAR$time <- as.POSIXct(PAR$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  }
  
  # convert data to xts objects and then merge.
  if (expanded == TRUE) {
    Rg.xts <- xts(Rg[,1:4],order.by=Rg$time)
  }  else {
    Rg.xts <- xts(Rg$Rg,order.by=Rg$time)
  }
  Rh.xts <- xts(Rh$Rh,order.by=Rh$time)
  Ta.xts <- xts(Ta$Ta,order.by=Ta$time)
  
  if (expanded == TRUE) {
    PAR.xts <- xts(PAR$PAR,order.by=PAR$time)
  }

  if (expanded == TRUE) {
    names(Rg.xts) <- c("SWdown","SWup","LWdown","LWup")
    names(PAR.xts) <- "PAR"
  } else {
    names(Rg.xts) <- "Rg"
  }
  
  names(Rh.xts) <- "Rh"
  names(Ta.xts) <- "Ta"

  minTime <- as.POSIXct(min(c(min(index(Rg.xts)),min(index(Rh.xts)),min(index(Ta.xts)),min(index(flux.xts)))),origin="1970-01-01")
  maxTime <- as.POSIXct(max(c(max(index(Rg.xts)),max(index(Rh.xts)),max(index(Ta.xts)),max(index(flux.xts)))),origin="1970-01-01")
  
  dummy.ts <- seq.POSIXt(minTime,maxTime,by=1800)
  dummy.data <- rep(NA,length(dummy.ts))
  
  dummy.xts <- xts(dummy.data,order.by=dummy.ts)
  
  # create bound xts
  if (expanded == TRUE) {
    all.data <- merge.xts(dummy.xts,Rg.xts,Rh.xts,Ta.xts,PAR.xts,flux.xts)
  } else {
    all.data <- merge.xts(dummy.xts,Rg.xts,Rh.xts,Ta.xts,flux.xts)
  }
  
  # calculate vpd from Tair and Rh
  all.data$vpd <- (100-all.data$Rh)/100*ifelse(all.data$Ta < 0,
                      exp(23.33086-6111.72784/(all.data$Ta+273.15)+0.15215*log(all.data$Ta+273.15)), # vapor over ice
                      exp(53.67957-6743.769/(all.data$Ta+273.15)-4.8451*log(all.data$Ta+273.15))) # vapor over liquid
  
  if (year != 9999) {
    all.data <- all.data[paste0(start.mon,"/",end.mon)]  
  }
  
  if (tzone == "PST") {
    indexTZ(all.data) <- "Etc/GMT+8"  
  } else if (tzone == "MST") {
    indexTZ(all.data) <- "Etc/GMT+7"
  } else if (tzone == "CST") {
    indexTZ(all.data) <- "Etc/GMT+6" 
  } else if (tzone == "EST") {
    indexTZ(all.data) <- "Etc/GMT+5"
  }
  
   #------------------------------------------------------------
  # load more data if we're going for the expanded package.
  # convert to MPI required format.
  out.data <- coredata(all.data)
  time.tmp <- index(all.data)
  
  # pull out time variables as required.
  yr <- year(time.tmp)
  doy <- yday(time.tmp)
  hr <- as.numeric(hour(time.tmp)+minute(time.tmp)/60)
  
  if (expanded == TRUE) {
    head1 <- c("Year","DoY","Hour","NEE","Fc","Sc","LH","Fw","Sw","H","FT","ST","SWup",
               "SWdown","LWup","LWdown","Tair","Tsoil","rH","VPD","Ustar","PAR")
    
    data.out <- data.frame(yr,doy,hr,out.data[,"nee"],out.data[,"Fc"],out.data[,"Sc"],
                           out.data[,"lhf"],out.data[,"Fw"],out.data[,"Sw"],
                           out.data[,"shf"],out.data[,"Ft"],out.data[,"St"],
                           out.data[,"SWup"],out.data[,"SWdown"],out.data[,"LWup"],out.data[,"LWdown"],
                           out.data[,"Ta"],rep(NA,nrow(out.data)),
                           out.data[,"Rh"],out.data[,"vpd"],out.data[,"ustar"],out.data[,"PAR"])
  } else {
  
    head1 <- c("Year","DoY","Hour","NEE","LH","H","Rg","Tair","Tsoil","rH","VPD","Ustar")
    head2 <- c("--","--","--","umol-2-s","Wm-2","Wm-2","Wm-2","degC","degC","%","hPa","ms-1")
    data.out <- data.frame(yr,doy,hr,out.data[,"nee"],
                           out.data[,"lhf"],out.data[,"shf"],out.data[,"Rg"],
                           out.data[,"Ta"],rep(NA,nrow(out.data)),
                           out.data[,"Rh"],out.data[,"vpd"],out.data[,"ustar"])
  }
  
  # return a data frame.
  names(data.out) <- head1
  
  if (median.filter == TRUE) {
    
    #---------------- Filter NEE ----------------
    # put in median deviation filter from Brock 86 / Starkenburg 2016.
    
    NEE.filt <- rollapply(data.out$NEE,7,median,na.rm=TRUE,fill=NA)
    NEE.logi <- abs(data.out$NEE - NEE.filt) > 10
    
    data.out$NEE[NEE.logi == TRUE] <- NA
    
    # remove points that are 5 sigma away from mean?
    NEE.mu <- mean(data.out$NEE,na.rm=TRUE)
    NEE.sd <- sd(data.out$NEE,na.rm=TRUE)

    NEE.oor <- (data.out$NEE < NEE.mu-4*NEE.sd) | (data.out$NEE > NEE.mu+4*NEE.sd)

    # set out of range values to missing:
    data.out$NEE[NEE.oor == TRUE] <- NA

    #---------------- Filter LH -----------------
    # remove points that are 5 sigma away from mean?
    LH.mu <- mean(data.out$LH,na.rm=TRUE)
    LH.sd <- sd(data.out$LH,na.rm=TRUE)
    
    LH.oor <- (data.out$LH < LH.mu-4*LH.sd) | (data.out$LH > LH.mu+4*LH.sd)
    
    # set out of range values to missing:
    data.out$LH[LH.oor == TRUE] <- NA

    # just filter LH    
    data.out$LH <- rollapply(data.out$LH,filt.width,median,fill=NA)
    
    #---------------- Filter H ------------------
    # remove points that are 5 sigma away from mean?
    H.mu <- mean(data.out$H,na.rm=TRUE)
    H.sd <- sd(data.out$H,na.rm=TRUE)
    
    H.oor <- (data.out$H < H.mu-4*H.sd) | (data.out$H > H.mu+4*H.sd)
    
    # set out of range values to missing:
    data.out$H[H.oor == TRUE] <- NA
    
    # just filter H    
    data.out$H <- rollapply(data.out$H,filt.width,median,fill=NA)
    
  }
 
  attr(data.out,"lat") <- lat 
  attr(data.out,"lon") <- lon 
  attr(data.out,"tzone") <- tzone

  #------------------------------------------------------------
  # write out data file if requested.
  if (write.to.file == TRUE & expanded == FALSE) {
    if (is.null(out.path)) {
      stop("Attempting to write to file, but no path specified!")
    } else {
      data.mpi.out <- do.call(rbind,list(head1,head2,data.out))
      write.table(data.mpi.out,
                  paste0(out.path,"/",Sys.Date(),"_",neon.site,"_",year,"_forREddyProc.txt"),
                  sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
    }
  }
    
  return(data.out)
}  
  
  
  